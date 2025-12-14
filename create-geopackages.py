import os
import time
import math
import requests
import geopandas as gpd
import pandas as pd
import osmnx as ox
from shapely.geometry import shape, box
from shapely.ops import unary_union
from shapely.validation import make_valid
from shapely.geometry.base import BaseGeometry
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import warnings
import argparse
import sys
warnings.filterwarnings("ignore", message="invalid value encountered in area", category=RuntimeWarning)

# ---------------- CONFIG ----------------
equal_area_crs = "EPSG:5070"   # Albers Equal Area Conic - accurate area
web_crs = "EPSG:4326"          # WGS84 lat/lon
meter_crs = "EPSG:3857"        # Web Mercator - good for meters & OSM
today = datetime.today().strftime("%Y-%m-%d")
output_dir = "unified_geopackages_no_missing"
os.makedirs(output_dir, exist_ok=True)

TEST_MODE = False
TEST_STATES = ["06", "48", "36"]  # CA, TX, NY for testing

# Data sources
COUNTY_URL = "https://www2.census.gov/geo/tiger/TIGER2025/COUNTY/tl_2025_us_county.zip"
PLACE_BASE_URL = "https://www2.census.gov/geo/tiger/TIGER2025/PLACE/tl_2025_{statefp}_place.zip"
PARKSERVE_URL = "https://server7.tplgis.org/arcgis7/rest/services/ParkServe/ParkServe_ProdNew/MapServer/2/query"
PARKSERVE_CONDITION = "(park_designation = 'LP' OR park_designation = 'LREC')"

# OSM tags
parks_tags = {"leisure": ["park", "pitch", "sports_centre"]}
gym_tags = {"leisure": ["fitness_centre", "sports_hall"]}
bldg_tags = {"building": True}
biz_tags = {
        "amenity": True, "shop": True, "office": True
    }
failed_cities = {}

# Exclude non-contiguous
EXCLUDED_STATEFPS = {"02", "15", "60", "66", "69", "72", "78"}  # AK, HI, territories

# Tiling & OSM settings
MAX_TILE_AREA_KM2 = 20.0        # Reduced for reliability
OSM_REQUEST_TIMEOUT = 30
OSM_MAX_RETRIES = 5  # Increase from 1
OSM_RETRY_BACKOFF = 2.0  # Increase backoff
TILE_SLEEP = 0.4                # Polite delay between OSM tiles

# Indoor proxy size
INDOOR_SIDE_FT = 70.0
INDOOR_SIDE_M = INDOOR_SIDE_FT * 0.3048
INDOOR_HALF_M = INDOOR_SIDE_M / 2.0

# ---------------- Utility functions ----------------

def safe_unary_union(gdf_or_geoms):
    """Safely compute unary_union with geometry validation"""
    if isinstance(gdf_or_geoms, gpd.GeoDataFrame):
        geoms = [make_valid(g) for g in gdf_or_geoms.geometry if g is not None and not g.is_empty]
    else:
        geoms = [make_valid(g) for g in gdf_or_geoms if g is not None and not g.is_empty]
    
    if not geoms:
        return None
    
    try:
        return unary_union(geoms)
    except Exception as e:
        print(f"  Warning: unary_union failed: {e}")
        return None
    
def safe_area(geom):
    if geom is None or geom.is_empty:
        return 0.0
    try:
        geom_v = make_valid(geom)
        a = geom_v.area
        return 0.0 if math.isnan(a) else a
    except Exception:
        return 0.0

# ---------------- helpers to handle shapely vs geopandas emptiness ----------------

def is_geom_empty(obj):
    """
    Return True if obj is empty.
    Works for:
      - GeoDataFrame / GeoSeries: uses .empty (pandas)
      - Shapely geometry: uses .is_empty
      - None: True
      - other objects: fallback to False
    """
    if obj is None:
        return True
    # GeoDataFrame / GeoSeries / pandas objects
    if hasattr(obj, "empty"):
        try:
            return bool(obj.empty)
        except Exception:
            # defensive fallback
            pass
    # shapely geometry
    if isinstance(obj, BaseGeometry):
        return bool(obj.is_empty)
    return False

def to_gdf_mask(mask, target_crs):
    """
    Convert a mask (shapely geometry or GeoDataFrame/GeoSeries) to a single-row GeoDataFrame
    in target_crs. If mask is already a GeoDataFrame/GeoSeries, it will be returned (reprojected).
    """
    if mask is None:
        return gpd.GeoDataFrame(geometry=[], crs=target_crs)
    if isinstance(mask, gpd.GeoDataFrame):
        return mask.to_crs(target_crs)
    if isinstance(mask, gpd.GeoSeries):
        gdf = gpd.GeoDataFrame(geometry=mask.values, crs=mask.crs)
        return gdf.to_crs(target_crs)
    if isinstance(mask, BaseGeometry):
        gdf = gpd.GeoDataFrame(geometry=[mask], crs=target_crs)
        return gdf.to_crs(target_crs)
    # unknown type -> empty
    return gpd.GeoDataFrame(geometry=[], crs=target_crs)

# ---------------- Robust safe_clip ----------------
def safe_clip(gdf, mask):
    """
    Robust clipping that accepts:
      - gdf: GeoDataFrame (expected)
      - mask: can be GeoDataFrame, GeoSeries, shapely geometry, or None
    Falls back to overlay, then to intersects-selection if overlay fails.
    Always returns a GeoDataFrame in the gdf.crs (or web_crs if gdf has no crs).
    """
    # quick sanity
    if gdf is None:
        return gpd.GeoDataFrame(geometry=[], crs=getattr(gdf, "crs", None))
    gdf_crs = getattr(gdf, "crs", None)

    # If gdf is empty, nothing to do
    if is_geom_empty(gdf):
        return gpd.GeoDataFrame(geometry=[], crs=gdf_crs)

    # Convert mask into a GeoDataFrame in same CRS as gdf for reliable overlay/clip
    try:
        mask_gdf = to_gdf_mask(mask, gdf_crs)
    except Exception:
        # fallback: create empty mask
        mask_gdf = gpd.GeoDataFrame(geometry=[], crs=gdf_crs)

    if is_geom_empty(mask_gdf):
        # nothing to clip with -> return empty with same schema/CRS as gdf
        return gpd.GeoDataFrame(columns=gdf.columns, geometry=[], crs=gdf_crs)

    # Try gpd.clip first (fast)
    try:
        clipped = gpd.clip(gdf, mask_gdf)
        clipped = clipped.set_crs(gdf_crs, allow_override=True)
        return clipped
    except Exception as e_clip:
        # fallback to overlay
        try:
            # ensure both are same CRS
            left = gdf.to_crs(gdf_crs)
            right = mask_gdf.to_crs(gdf_crs)
            over = gpd.overlay(left, right, how="intersection")
            over = over.set_crs(gdf_crs, allow_override=True)
            return over
        except Exception as e_over:
            # final fallback: keep features whose geometry intersects any mask geometry
            try:
                # Build unary_union of mask geometries for faster intersects tests
                mask_union = unary_union([geom for geom in mask_gdf.geometry if geom is not None])
                subset = gdf[gdf.geometry.intersects(mask_union)]
                subset = subset.set_crs(gdf_crs, allow_override=True)
                return subset
            except Exception as e_final:
                # give up - return empty with original columns/CRS
                return gpd.GeoDataFrame(columns=gdf.columns, geometry=[], crs=gdf_crs)


def fetch_parkserve_features(state_fips):
    features = []
    offset = 0
    while True:
        params = {
            "where": f"Park_State_FIPS = '{state_fips}' AND {PARKSERVE_CONDITION}",
            "outFields": "*",
            "f": "geojson",
            "resultOffset": offset,
            "resultRecordCount": 1000,
        }
        try:
            r = requests.get(PARKSERVE_URL, params=params, timeout=60)
            r.raise_for_status()
            gj = r.json()
            page = gj.get("features", [])
            if not page:
                break
            features.extend(page)
            if len(page) < 1000:
                break
            offset += 1000
        except Exception as e:
            print(f"ParkServe fetch error for {state_fips}: {e}")
            break
    
    # Validate geometries before returning
    if not features:
        return []
    
    valid_features = []
    for f in features:
        if f.get("geometry"):
            try:
                geom = shape(f["geometry"])
                if geom and not geom.is_empty:
                    geom = make_valid(geom)
                    if not geom.is_empty:
                        f["geometry"] = geom.__geo_interface__
                        valid_features.append(f)
            except Exception:
                continue
    
    return valid_features


def tile_polygon_to_grid(geom, max_tile_area_km2=MAX_TILE_AREA_KM2):
    if geom is None or geom.is_empty:
        return []
    minx, miny, maxx, maxy = geom.bounds
    total_area = (maxx - minx) * (maxy - miny)
    if total_area <= (max_tile_area_km2 * 1e6):
        return [geom]
    tile_side = math.sqrt(max_tile_area_km2 * 1e6)
    nx = max(1, int(math.ceil((maxx - minx) / tile_side)))
    ny = max(1, int(math.ceil((maxy - miny) / tile_side)))
    tiles = []
    dx = (maxx - minx) / nx
    dy = (maxy - miny) / ny
    for i in range(nx):
        for j in range(ny):
            tx1 = minx + i * dx
            ty1 = miny + j * dy
            tx2 = tx1 + dx
            ty2 = ty1 + dy
            tile = box(tx1, ty1, tx2, ty2)
            inter = geom.intersection(tile)
            if not inter.is_empty:
                tiles.append(inter)
    return tiles


def osm_query_with_retries(tags, polygon_ll, state_fips=None, city_name=None):
    attempt = 0
    while attempt < OSM_MAX_RETRIES:
        try:
            gdf = ox.features_from_polygon(polygon_ll, tags)
            if gdf is None or gdf.empty:
                return gpd.GeoDataFrame(geometry=[], crs=web_crs)
            gdf = gdf.set_crs(web_crs, allow_override=True)
            return gdf
        except Exception as e:
            msg = str(e).lower()
            if "no matching features" in msg or "no features" in msg or "empty" in msg:
                return gpd.GeoDataFrame(geometry=[], crs=web_crs)
            attempt += 1
            wait = OSM_RETRY_BACKOFF ** attempt
            print(f"  OSM query failed (attempt {attempt}) for {city_name} ({state_fips}): {e}. Retrying in {wait}s...")
            time.sleep(wait)

    # After max retries, log the failure
    print(f"  OSM query failed after max retries for {city_name} ({state_fips})")
    return None


# Helper: single geometry conversion (avoid repeated GeoSeries creation)
def geom_to_crs(geom, target_crs):
    if geom is None or geom.is_empty:
        return None
    return gpd.GeoSeries([geom], crs=equal_area_crs).to_crs(target_crs).iloc[0]

# ---------------- Optimized OSM fetch ----------------
def fetch_osm_by_tiling(tags, geom_proj, state_fips=None, city_name=None, buffer_m=50):
    """Fetch OSM features by tiling, optimized to avoid redundant CRS conversions."""
    if geom_proj is None or geom_proj.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=web_crs)
    
    # Tile directly in equal-area CRS (meters)
    tiles = tile_polygon_to_grid(geom_proj, MAX_TILE_AREA_KM2)
    parts = []

    for tile in tiles:
        tile_buf = tile.buffer(buffer_m)  # buffer directly in equal-area CRS
        # Convert buffered tile to WGS84 only when querying OSM
        tile_ll = geom_to_crs(tile_buf, web_crs)
        gdf_tile = osm_query_with_retries(
            tags,
            tile_ll,
            state_fips=state_fips,
            city_name=city_name
        )
        if not gdf_tile.empty:
            gdf_tile['geometry'] = gdf_tile.geometry.apply(make_valid)
            gdf_tile = gdf_tile[~gdf_tile.geometry.is_empty]
            if not gdf_tile.empty:
                # Clip in WGS84 first
                gdf_tile = safe_clip(gdf_tile, tile_ll)
                # Convert to equal-area CRS for union / area calculations
                gdf_tile = gdf_tile.to_crs(equal_area_crs)
                parts.append(gdf_tile)
        time.sleep(TILE_SLEEP)

    if not parts:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    combined = pd.concat(parts, ignore_index=True)
    combined = combined[~combined.geometry.duplicated()]

    # Clip back to original city boundary (equal-area CRS)
    combined = safe_clip(combined, geom_proj)
    return combined

def process_city_merged(city, state_fips, parks_ps):
    """
    Process a single city:
      - ParkServe
      - OSM parks, gyms, buildings, and other businesses
      - Indoor gyms clipped to building and nearby businesses
    Returns dict for final GeoDataFrame.
    """
    city_geom = city.geometry
    city_name = city.get("NAME", city.get("NAME10", "Unknown"))

    # ---------------- ParkServe ----------------
    ps_clipped = safe_clip(parks_ps, city_geom) if not parks_ps.empty else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    # ---------------- OSM Fetch ----------------
    # Parks
    osm_parks_gdf = fetch_osm_by_tiling(
    parks_tags,
    city_geom,
    state_fips=state_fips,
    city_name=city_name
    )
    osm_parks_gdf = osm_parks_gdf[osm_parks_gdf.geometry.type.isin(["Polygon", "MultiPolygon"])]

    # Gyms
    gym_gdf = fetch_osm_by_tiling(
        gym_tags,
        city_geom,
        state_fips=state_fips,
        city_name=city_name
    )

    # Buildings
    bldg_gdf = fetch_osm_by_tiling(
        bldg_tags,
        city_geom,
        state_fips=state_fips,
        city_name=city_name
    )

    # Other businesses
    biz_gdf = fetch_osm_by_tiling(
        biz_tags,
        city_geom,
        state_fips=state_fips,
        city_name=city_name
    )

    # Validate geometries
    for gdf in [osm_parks_gdf, gym_gdf, bldg_gdf, biz_gdf]:
        if not gdf.empty:
            gdf['geometry'] = gdf.geometry.apply(make_valid)
            gdf = gdf[~gdf.geometry.is_empty]
    
    # Separate gyms into points and polygons
    if not gym_gdf.empty:
        gdf_gym_points = gym_gdf[gym_gdf.geometry.type.isin(["Point", "MultiPoint"])].copy()
        gdf_gym_polys = gym_gdf[gym_gdf.geometry.type.isin(["Polygon", "MultiPolygon"])].copy()
    else:
        gdf_gym_points = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
        gdf_gym_polys = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
    
    # Convert all to meter CRS for distance/clip operations
    gdf_gym_points = gdf_gym_points.to_crs(meter_crs)
    gdf_gym_polys = gdf_gym_polys.to_crs(meter_crs)
    gdf_bldgs = bldg_gdf.to_crs(meter_crs)
    gdf_biz = biz_gdf.to_crs(meter_crs)

    # Assign building IDs for gym/building join
    gdf_bldgs = gdf_bldgs.reset_index(drop=True)
    gyms_in_bldg = gpd.sjoin(gdf_gym_points, gdf_bldgs[["geometry"]], how="inner", predicate="within")
    gyms_in_bldg.rename(columns={"index_right": "building_id"}, inplace=True)

    biz_in_bldg = gpd.sjoin(gdf_biz, gdf_bldgs[["geometry"]], how="inner", predicate="within")

    # Handle dynamic naming of the building index column
    join_col = None
    for col in biz_in_bldg.columns:
        if col.startswith("index_") or col.endswith("_right"):
            join_col = col
            break
    
    if join_col is None:
        raise KeyError("Could not find building index column in business join result")
    
    # Group businesses by building ID
    biz_in_bldg.rename(columns={join_col: "building_id"}, inplace=True)
    biz_grouped = biz_in_bldg.groupby("building_id")

    # --- 7. Create indoor zones for point-based gyms ---
    zones = []
    half = INDOOR_HALF_M
    
    for idx, row in gyms_in_bldg.iterrows():
        pt = row.geometry
        bldg_geom = gdf_bldgs.loc[row["building_id"]].geometry
    
        # Create 61x61 ft square, clipped to building
        x, y = pt.x, pt.y
        square = box(x - half, y - half, x + half, y + half)
        zone = square.intersection(bldg_geom)
    
        # Clip by nearby businesses in same building
        if row["building_id"] in biz_grouped.groups:
            biz_points = biz_grouped.get_group(row["building_id"]).geometry
            if not biz_points.empty:
                min_dist = biz_points.distance(pt).min()
                if pd.notna(min_dist) and min_dist > 0:
                    half_dist = min(min_dist / 2.0, half)
                    zone = zone.intersection(pt.buffer(half_dist)).intersection(bldg_geom)
    
        # Preserve attributes from the original point
        zones.append({
            "geometry": zone,
            "gym_name": row.get("name", None),
        })
    
    # Create GeoDataFrame with explicit geometry column
    if zones:
        gdf_zones = gpd.GeoDataFrame(zones, crs=meter_crs, geometry='geometry')
    else:
        gdf_zones = gpd.GeoDataFrame(geometry=[], crs=meter_crs)

   # Convert CRS once
    if not gdf_zones.empty:
        indoor_gdf = gdf_zones.to_crs(equal_area_crs)
    else:
        indoor_gdf = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    if not gdf_gym_polys.empty:
        gym_polys = gdf_gym_polys.to_crs(equal_area_crs)
    else:
        gym_polys = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    # Clip everything in one consistent way
    if not indoor_gdf.empty:
        indoor_gdf = safe_clip(indoor_gdf, city_geom)
        indoor_gdf = indoor_gdf[~indoor_gdf.geometry.is_empty]

    if not gym_polys.empty:
        gym_polys = safe_clip(gym_polys, city_geom)
        gym_polys = gym_polys[~gym_polys.geometry.is_empty]

        # Standardize columns
        gym_polys = gym_polys[['geometry'] + (['name'] if 'name' in gym_polys else [])].copy()
        gym_polys.rename(columns={'name': 'gym_name'}, inplace=True)
        if 'gym_name' not in gym_polys:
            gym_polys['gym_name'] = None

    # Combine
    if not gym_polys.empty:
        indoor_gdf = pd.concat([indoor_gdf, gym_polys], ignore_index=True)

    
    # Ensure we have a valid GeoDataFrame
    if indoor_gdf.empty:
        indoor_gdf = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
    else:
        indoor_gdf = indoor_gdf.set_crs(equal_area_crs, allow_override=True)


    # ---------------- Unions and areas ----------------
    components = []

    for layer in (ps_clipped, osm_parks_gdf, indoor_gdf):
        if not layer.empty:
            components.append(safe_unary_union(layer))

    unified_geom = safe_unary_union(components) if components else None
    area_total = safe_area(unified_geom)
    to_4326_gdf = lambda gdf: gdf.to_crs(web_crs) if not gdf.empty else None
    # Convert to WGS84 for saving
    to4326 = lambda geom: geom_to_crs(geom, web_crs) if geom is not None else None

    #add attributes

    if not ps_clipped.empty:
        ps_city = ps_clipped.copy()
        ps_city["city_name"] = city_name
        ps_city["state_fips"] = state_fips
        ps_city["date"] = today
    else:
        ps_city = None

    if not osm_parks_gdf.empty:
        osm_city = gpd.clip(osm_parks_gdf, city_geom)

        if not osm_city.empty:
            osm_city = osm_city.copy()
            osm_city["city_name"] = city_name
            osm_city["state_fips"] = state_fips
            osm_city["date"] = today
    else:
        osm_city = None

    if not indoor_gdf.empty:
        indoor_city = gpd.clip(indoor_gdf, city_geom)

        if not indoor_city.empty:
            indoor_city = indoor_city.copy()
            indoor_city["city_name"] = city_name
            indoor_city["state_fips"] = state_fips
            indoor_city["date"] = today
    else:
        indoor_city = None



    print("  Processed city:", city_name)
    return {
        "data": {
            "geometry": to4326(unified_geom),
            "area_m2_parks": area_total,
            "date": today,
            "state_fips": state_fips,
            "city_name": city_name,
        },
        "layers": {
            "parkserve": to_4326_gdf(ps_city),
            "osm": to_4326_gdf(osm_city),
            "indoor": to_4326_gdf(indoor_city),
        }
}




# ---------------- Main per-state processing ----------------
def process_state(state_fips, cities_state, max_workers=None):
    """
    Process a state using parallel execution per city.
    parks_ps: ParkServe features pre-fetched
    """
    try:
        print(f"Processing state {state_fips}...")

        # ---------------- ParkServe ----------------
        features = fetch_parkserve_features(state_fips)
        parks_ps = gpd.GeoDataFrame(
            geometry=[shape(f["geometry"]) for f in features if f.get("geometry")],
            crs=web_crs
        ).to_crs(equal_area_crs) if features else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

        if not parks_ps.empty:
            parks_ps['geometry'] = parks_ps.geometry.apply(make_valid)
            parks_ps = parks_ps[~parks_ps.geometry.is_empty]

        results = []
        parkserve_layers = []
        osm_layers = []
        indoor_layers = []

        # ---------------- Parallel city processing ----------------
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_city_merged, city, state_fips, parks_ps): idx 
                       for idx, city in cities_state.iterrows()}

            for future in as_completed(futures):
                res = future.result()
                results.append(res["data"])

                if res["layers"]["parkserve"] is not None:
                    parkserve_layers.append(res["layers"]["parkserve"])
                if res["layers"]["osm"] is not None:
                    osm_layers.append(res["layers"]["osm"])
                if res["layers"]["indoor"] is not None:
                    indoor_layers.append(res["layers"]["indoor"])


        if not results:
            print(f"No cities processed for state {state_fips}")
            return state_fips

        # ---------------- Build final GeoDataFrame ----------------
        state_gdf = gpd.GeoDataFrame(results, crs="EPSG:4326")

        # Save unified layer
        gpkg_path = os.path.join(output_dir, f"state_{state_fips}_unified.gpkg")
        csv_path = os.path.join(output_dir, f"state_{state_fips}_unified.csv")

        unified = state_gdf.drop(columns=["geometry_parkserve", "geometry_osm", "geometry_indoor"], errors="ignore")
        unified = unified.set_geometry("geometry")
        unified = unified.set_crs("EPSG:4326", allow_override=True)
        unified.to_file(gpkg_path, layer="unified_park_area", driver="GPKG")

        # Save component layers
        # ---------------- Save component layers ----------------

        if parkserve_layers:
            parkserve_gdf = gpd.GeoDataFrame(
                pd.concat(parkserve_layers, ignore_index=True),
                crs="EPSG:4326"
            )
            parkserve_gdf.to_file(gpkg_path, layer="parkserve_geom", driver="GPKG")

        if osm_layers:
            osm_gdf = gpd.GeoDataFrame(
                pd.concat(osm_layers, ignore_index=True),
                crs="EPSG:4326"
            )
            osm_gdf.to_file(gpkg_path, layer="osm_geom", driver="GPKG")

        if indoor_layers:
            indoor_gdf = gpd.GeoDataFrame(
                pd.concat(indoor_layers, ignore_index=True),
                crs="EPSG:4326"
            )
            indoor_gdf.to_file(gpkg_path, layer="indoor_gyms", driver="GPKG")

        # CSV summary
        state_gdf.drop(columns=["geometry", "geometry_parkserve", "geometry_osm", "geometry_indoor"], errors="ignore").to_csv(csv_path, index=False)

        print(f"State {state_fips} completed: {len(state_gdf)} places")
        return state_fips

    except Exception as e:
        print(f"Error in state {state_fips}: {e}")
        import traceback
        traceback.print_exc()
        return None

def load_cities_for_state(state_fips):
    url = PLACE_BASE_URL.format(statefp=state_fips)
    df = gpd.read_file(url)
    df = df.to_crs(equal_area_crs)
    df['NAME'] = df['NAME'].astype(str)

    # normalize column names
    df['PLACEFP'] = df.get('PLACEFP', df.get('PLACEFP10'))
    return df

def run_missing_cities(missing_df, max_workers=6):
    """Run only the missing cities and append to the existing state gpkg."""
    states = sorted(missing_df.state_fips.unique())

    for state in states:
        print(f"\n=== Processing missing cities for state {state} ===")

        # Load all cities for that state
        cities = load_cities_for_state(state)

        # Filter only needed ones
        miss = missing_df[missing_df.state_fips == state]

        sub = cities[cities['NAME'].str.lower().isin(
            miss.city_name.str.lower()
        )].copy()

        if sub.empty:
            print(f"No matching cities found in TIGER for {state}")
            continue

        # Run your original state processor,
        # but pass only the subset of cities
        process_state(state, sub, max_workers=max_workers)

        print(f"Finished missing subset for {state}")


def main():
    args = parse_args()

    # Single state provided by SLURM array
    state_fips = args.state[0].zfill(2)
    print(f"Processing state {state_fips}")

    print("Loading counties...")
    counties_gdf = gpd.read_file(COUNTY_URL).to_crs(equal_area_crs)

    if state_fips not in counties_gdf["STATEFP"].unique():
        raise RuntimeError(f"State {state_fips} not found in county dataset")

    if state_fips in EXCLUDED_STATEFPS:
        raise RuntimeError(f"State {state_fips} is excluded")

    print("Loading city/place boundaries...")
    url = PLACE_BASE_URL.format(statefp=state_fips)
    cities_state = gpd.read_file(url).to_crs(equal_area_crs)
    cities_state["statefp"] = state_fips

    print(f"Loaded {len(cities_state)} cities for state {state_fips}")

    # Direct single-state execution
    print(f"Running process_state({state_fips}, cities_state)...")
    result = process_state(state_fips, cities_state)

    print(f"Finished state {state_fips}: {result}")
    if failed_cities:
        print("Failed cities per state:")
    for st, cities in failed_cities.items():
        print(f"State {st}: {cities}")

    # Save to CSV for retry later
    pd.DataFrame([
        {"state_fips": st, "city_name": city} 
        for st, cities in failed_cities.items() 
        for city in cities
    ]).to_csv(f"failed_cities_to_retry_{state_fips}.csv", index=False)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build unified geopackages. Optionally restrict to specific states."
    )

    parser.add_argument(
        "--state",
        nargs="+",
        help=(
            "State FIPS codes to process. Examples:\n"
            "  --state 30             (single)\n"
            "  --state 30 32 56       (multiple)\n"
            "  --state all            (process all states)\n"
        ),
        default=["all"]
    )

    return parser.parse_args()

if __name__ == "__main__":
    main()