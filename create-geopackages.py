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
OSM_PARKS_TAGS = {"leisure": ["park", "pitch"]}
OSM_SPORTS_CENTRE_TAGS = {"leisure": ["sports_centre"]}
OSM_GYM_TAGS = {"leisure": ["fitness_centre", "sports_hall"]}
bldg_tags = {"building": True}
biz_tags = {
        "amenity": True, "shop": True, "office": True
    }

# Exclude non-contiguous
EXCLUDED_STATEFPS = {"02", "15", "60", "66", "69", "72", "78"}  # AK, HI, territories

# Tiling & OSM settings
MAX_TILE_AREA_KM2 = 20.0        # Reduced for reliability
OSM_REQUEST_TIMEOUT = 30
OSM_MAX_RETRIES = 8  # Increase from 1
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
    return gpd.GeoDataFrame(geometry=[], crs=web_crs)


# Helper: single geometry conversion (avoid repeated GeoSeries creation)
def geom_to_crs(geom, target_crs):
    if geom is None or geom.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
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

def keep_name_and_sport(gdf, name_col="name", sport_col="sport"):
    """
    Reduce a GeoDataFrame to geometry + name + sport.
    Ensures name and sport exist and are strings.
    """
    if gdf is None or gdf.empty:
        return gpd.GeoDataFrame(geometry=[], crs=getattr(gdf, "crs", None))
    gdf = gdf.copy()
    
    # Handle name column
    if name_col in gdf.columns:
        gdf[name_col] = gdf[name_col].astype(str)
    else:
        gdf[name_col] = None
    
    # Handle sport column
    if sport_col in gdf.columns:
        gdf[sport_col] = gdf[sport_col].astype(str)
    else:
        gdf[sport_col] = None
    
    # Keep only geometry, name, and sport
    gdf = gdf[[name_col, sport_col, "geometry"]]
    return gdf


def process_city_merged(city, state_fips, parks_ps):
    """
    Process a single city:
      - ParkServe
      - OSM parks, gyms, buildings, and other businesses
      - Indoor gyms clipped to building and nearby businesses
    Returns dict with three layers: parks, outdoor_sports, and indoor_sports.
    """
    city_geom = city.geometry
    city_name = city.get("NAME", city.get("NAME10", "Unknown"))

    # ---------------- ParkServe ----------------
    ps_clipped = safe_clip(parks_ps, city_geom) if not parks_ps.empty else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    # ---------------- OSM Fetch ----------------
    # --- OSM parks (polygons only)
    osm_parks = fetch_osm_by_tiling(
        OSM_PARKS_TAGS,
        city_geom,
        state_fips=state_fips,
        city_name=city_name
    )
    osm_parks = osm_parks[
        osm_parks.geometry.type.isin(["Polygon", "MultiPolygon"])
    ].copy()

    # --- OSM sports centres (polygons only)
    osm_sports_centre = fetch_osm_by_tiling(
        OSM_SPORTS_CENTRE_TAGS,
        city_geom,
        state_fips=state_fips,
        city_name=city_name
    )
    osm_sports_centre = osm_sports_centre[
        osm_sports_centre.geometry.type.isin(["Polygon", "MultiPolygon"])
    ].copy()

    # --- OSM gyms (points + polygons)
    gym_gdf = fetch_osm_by_tiling(
        OSM_GYM_TAGS,
        city_geom,
        state_fips=state_fips,
        city_name=city_name
    )

    osm_gym_points = gym_gdf[
        gym_gdf.geometry.type.isin(["Point", "MultiPoint"])
    ].copy()

    osm_gym_polys = gym_gdf[
        gym_gdf.geometry.type.isin(["Polygon", "MultiPolygon"])
    ].copy()

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
    for gdf in [osm_sports_centre, osm_parks, osm_gym_points, osm_gym_polys, bldg_gdf, biz_gdf]:
        if not gdf.empty:
            gdf['geometry'] = gdf.geometry.apply(make_valid)
            gdf = gdf[~gdf.geometry.is_empty]
    
    # Convert all to meter CRS for distance/clip operations
    osm_gym_points = osm_gym_points.to_crs(meter_crs)
    osm_gym_polys = osm_gym_polys.to_crs(meter_crs)
    gdf_bldgs = bldg_gdf.to_crs(meter_crs)
    gdf_biz = biz_gdf.to_crs(meter_crs)

    # Assign building IDs for gym/building join
    gdf_bldgs = gdf_bldgs.reset_index(drop=True)
    gyms_in_bldg = gpd.sjoin(osm_gym_points, gdf_bldgs[["geometry"]], how="inner", predicate="within")
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

    # --- Create indoor zones for point-based gyms ---
    zones = []
    half = INDOOR_HALF_M
    
    for idx, row in gyms_in_bldg.iterrows():
        pt = row.geometry
        bldg_geom = gdf_bldgs.loc[row["building_id"]].geometry
    
        # Create 70x70 ft square, clipped to building
        x, y = pt.x, pt.y
        square = box(x - half, y - half, x + half, y + half)
        zone = square.intersection(bldg_geom)
    
        # Clip by nearby businesses in same building, circular
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
            "name": row.get("name", None),
        })
    
    # Create GeoDataFrame with explicit geometry column
    if zones:
        gdf_zones = gpd.GeoDataFrame(zones, crs=meter_crs, geometry='geometry')
    else:
        gdf_zones = gpd.GeoDataFrame(geometry=[], crs=meter_crs)
    
    # Convert CRS to equal_area_crs
    if not gdf_zones.empty:
        indoor_zones = gdf_zones.to_crs(equal_area_crs)
    else:
        indoor_zones = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
    
    # Convert gym polygons to equal_area_crs
    if not osm_gym_polys.empty:
        gym_polys = osm_gym_polys.to_crs(equal_area_crs)
        # Clip to city boundary
        gym_polys = safe_clip(gym_polys, city_geom)
        gym_polys = gym_polys[~gym_polys.geometry.is_empty]
        
        # Standardize columns
        if not gym_polys.empty:
            gym_polys = gym_polys[['geometry'] + (['name'] if 'name' in gym_polys else [])].copy()
            if 'name' not in gym_polys:
                gym_polys['name'] = None
    else:
        gym_polys = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    # Clip indoor zones to city boundary
    if not indoor_zones.empty:
        indoor_zones = safe_clip(indoor_zones, city_geom)
        indoor_zones = indoor_zones[~indoor_zones.geometry.is_empty]

        # First, ensure sports centres are in the correct CRS
    if not osm_sports_centre.empty:
        osm_sports_centre = osm_sports_centre.to_crs(equal_area_crs)
        osm_sports_centre = safe_clip(osm_sports_centre, city_geom)
        osm_sports_centre = osm_sports_centre[~osm_sports_centre.geometry.is_empty]

    # Separate indoor vs outdoor components of sports centres
    indoor_sports_centres = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
    outdoor_sports_centres = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    if not osm_sports_centre.empty and not gdf_bldgs.empty:
        # Convert buildings to equal_area_crs and create union
        gdf_bldgs_ea = gdf_bldgs.to_crs(equal_area_crs)
        buildings_union = gdf_bldgs_ea.union_all()
        
        indoor_parts = []
        outdoor_parts = []
        
        for idx, row in osm_sports_centre.iterrows():
            geom = row.geometry
            
            # Indoor component: intersection with buildings
            indoor_part = geom.intersection(buildings_union)
            if not indoor_part.is_empty:
                indoor_record = {'geometry': indoor_part}
                if 'sport' in row.index:
                    indoor_record['sport'] = row['sport']
                indoor_parts.append(indoor_record)
            
            # Outdoor component: difference with buildings
            outdoor_part = geom.difference(buildings_union)
            if not outdoor_part.is_empty:
                outdoor_record = {'geometry': outdoor_part}
                if 'sport' in row.index:
                    outdoor_record['sport'] = row['sport']
                outdoor_parts.append(outdoor_record)
        
        if indoor_parts:
            indoor_sports_centres = gpd.GeoDataFrame(indoor_parts, crs=equal_area_crs)
        if outdoor_parts:
            outdoor_sports_centres = gpd.GeoDataFrame(outdoor_parts, crs=equal_area_crs)
            
    elif not osm_sports_centre.empty:
        # No buildings data, treat all as outdoor
        outdoor_sports_centres = osm_sports_centre.copy()

    # --- COMBINE INDOOR SPORTS: gym polygons + indoor zones + indoor sports centres ---
    indoor_sports_list = []
    if not indoor_zones.empty:
        indoor_sports_list.append(indoor_zones)
    if not gym_polys.empty:
        indoor_sports_list.append(gym_polys)
    if not indoor_sports_centres.empty:
        indoor_sports_list.append(indoor_sports_centres)

    if indoor_sports_list:
        indoor_sports = pd.concat(indoor_sports_list, ignore_index=True)
        indoor_sports = indoor_sports.set_crs(equal_area_crs, allow_override=True)
    else:
        indoor_sports = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    # --- OUTDOOR SPORTS: outdoor sports centres ---
    outdoor_sports = outdoor_sports_centres.copy()

    # --- PARKS: union of ParkServe and OSM parks, minus outdoor sports only ---
    parks_list = []
    if not ps_clipped.empty:
        parks_list.append(ps_clipped)
    if not osm_parks.empty:
        parks_list.append(osm_parks)

    if parks_list:
        parks_combined = pd.concat(parks_list, ignore_index=True)
        parks_combined = parks_combined.set_crs(equal_area_crs, allow_override=True)
        
        # Remove areas that overlap with outdoor sports
        if not outdoor_sports.empty:
            outdoor_union = outdoor_sports.union_all()
            parks_combined['geometry'] = parks_combined.geometry.apply(
                lambda geom: geom.difference(outdoor_union) if geom and not geom.is_empty else geom
            )
            parks_combined = parks_combined[~parks_combined.geometry.is_empty]
        
        parks = parks_combined
    else:
        parks = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    

    # ---------------- Add attributes and convert to WGS84 ----------------
    def to_4326_gdf(gdf):
        if gdf is None or gdf.empty:
            return None
        return gdf.to_crs(web_crs)

    def add_attributes(gdf, city_name, state_fips):
        if gdf is None or gdf.empty:
            return None
        gdf = gdf.copy()
        gdf["city_name"] = city_name
        gdf["state_fips"] = state_fips
        return gdf

    parks = keep_name_and_sport(parks)
    outdoor_sports = keep_name_and_sport(outdoor_sports)
    indoor_sports = keep_name_and_sport(indoor_sports)

    parks = add_attributes(parks, city_name, state_fips)
    outdoor_sports = add_attributes(outdoor_sports, city_name, state_fips)
    indoor_sports = add_attributes(indoor_sports, city_name, state_fips)

    print("  Processed city:", city_name)
    return {
        "layers": {
            "parks": to_4326_gdf(parks),
            "outdoor_sports_facilities": to_4326_gdf(outdoor_sports),
            "indoor_sports_facilities": to_4326_gdf(indoor_sports)
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

        parks = []
        outdoor_sports_centres = []
        indoor_sports_facilities =  []


        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(process_city_merged, city, state_fips, parks_ps): idx
                    for idx, city in cities_state.iterrows()}

            for future in as_completed(futures):
                res = future.result()

                # old "results" collection for unified data is no longer needed
                layers = res["layers"]

                if layers["parks"] is not None:
                    parks.append(layers["parks"])
                if layers["outdoor_sports_facilities"] is not None:
                    outdoor_sports_centres.append(layers["outdoor_sports_facilities"])
                if layers["indoor_sports_facilities"] is not None:
                    indoor_sports_facilities.append(layers["indoor_sports_facilities"])


        gpkg_path = os.path.join(output_dir, f"state_{state_fips}_components.gpkg")

        layer_dict = {
            "parks": parks,
            "outdoor_sports_facilities": outdoor_sports_centres,
            "indoor_sports_facilities": indoor_sports_facilities,
        }

        for layer_name, gdfs in layer_dict.items():
            if gdfs:
                gdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), crs="EPSG:4326")
                gdf.to_file(gpkg_path, layer=layer_name, driver="GPKG")

        print(f"State {state_fips} completed: {len(cities_state)} cities processed")
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