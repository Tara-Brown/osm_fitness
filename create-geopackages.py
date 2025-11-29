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
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import warnings
warnings.filterwarnings("ignore", message="invalid value encountered in area", category=RuntimeWarning)

# ---------------- CONFIG ----------------
equal_area_crs = "EPSG:5070"   # Albers Equal Area Conic - accurate area
web_crs = "EPSG:4326"          # WGS84 lat/lon
meter_crs = "EPSG:3857"        # Web Mercator - good for meters & OSM
today = datetime.today().strftime("%Y-%m-%d")
output_dir = "unified_geopackages_"
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

# Exclude non-contiguous
EXCLUDED_STATEFPS = {"02", "15", "60", "66", "69", "72", "78"}  # AK, HI, territories

# Tiling & OSM settings
MAX_TILE_AREA_KM2 = 20.0        # Reduced for reliability
OSM_REQUEST_TIMEOUT = 30
OSM_MAX_RETRIES = 3
OSM_RETRY_BACKOFF = 2.0
TILE_SLEEP = 0.4                # Polite delay between OSM tiles

# Indoor proxy size
INDOOR_SIDE_FT = 70.0
INDOOR_SIDE_M = INDOOR_SIDE_FT * 0.3048
INDOOR_HALF_M = INDOOR_SIDE_M / 2.0

# ---------------- Utility functions ----------------
def safe_area(geom):
    if geom is None or geom.is_empty:
        return 0.0
    try:
        geom_v = make_valid(geom)
        a = geom_v.area
        return 0.0 if math.isnan(a) else a
    except Exception:
        return 0.0

def safe_clip(gdf, mask):
    """Robust clipping that falls back to overlay on failure."""
    if gdf.empty or mask.empty:
        return gpd.GeoDataFrame(geometry=[], crs=gdf.crs)
    try:
        return gpd.clip(gdf, mask)
    except Exception:
        try:
            mask_gdf = gpd.GeoDataFrame(geometry=[mask], crs=gdf.crs) if not isinstance(mask, gpd.GeoDataFrame) else mask
            return gpd.overlay(gdf, mask_gdf, how="intersection")
        except Exception:
            return gdf[gdf.intersects(mask)]

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
    return features

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

def osm_query_with_retries(tags, polygon_ll):
    attempt = 0
    while attempt < OSM_MAX_RETRIES:
        try:
            gdf = ox.geometries_from_polygon(polygon_ll, tags)
            if gdf is None or gdf.empty:
                return gpd.GeoDataFrame(geometry=[], crs=web_crs)
            gdf = gdf.set_crs(web_crs, allow_override=True)
            return gdf
        except Exception as e:
            attempt += 1
            wait = OSM_RETRY_BACKOFF ** attempt
            print(f"  OSM query failed (attempt {attempt}): {e}. Retrying in {wait}s...")
            time.sleep(wait)
    print("  OSM query failed after max retries")
    return gpd.GeoDataFrame(geometry=[], crs=web_crs)

def fetch_osm_by_tiling(tags, geom_proj, buffer_m=50):
    if geom_proj is None or geom_proj.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=web_crs)
    geom_m = gpd.GeoSeries([geom_proj], crs=equal_area_crs).to_crs(meter_crs).iloc[0]
    tiles = tile_polygon_to_grid(geom_m, MAX_TILE_AREA_KM2)
    parts = []
    for tile in tiles:
        tile_buf = tile.buffer(buffer_m)
        tile_ll = gpd.GeoSeries([tile_buf], crs=meter_crs).to_crs(web_crs).iloc[0]
        gdf_tile = osm_query_with_retries(tags, tile_ll)
        if not gdf_tile.empty:
            gdf_tile = safe_clip(gdf_tile, tile_ll)
            parts.append(gdf_tile)
        time.sleep(TILE_SLEEP)  # Be nice to OSM
    if not parts:
        return gpd.GeoDataFrame(geometry=[], crs=web_crs)
    combined = pd.concat(parts, ignore_index=True)
    combined = combined.set_crs(web_crs, allow_override=True)
    # Fast deduplication in meter CRS
    combined_m = combined.to_crs(meter_crs)
    combined_m = combined_m[~combined_m.geometry.duplicated()]
    return combined_m.to_crs(web_crs)

# ---------------- Indoor handling (70 ft squares, clipped to building) ----------------
def load_osm_indoor_for_city(city_geom_proj):
    if city_geom_proj is None or city_geom_proj.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    city_geom_m = gpd.GeoSeries([city_geom_proj], crs=equal_area_crs).to_crs(meter_crs).iloc[0]

    # Fetch data
    gyms_ll = fetch_osm_by_tiling(gym_tags, city_geom_proj)
    bldgs_ll = fetch_osm_by_tiling(bldg_tags, city_geom_proj)

    if gyms_ll.empty or bldgs_ll.empty:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)  # No buildings → no indoor

    gyms_m = gyms_ll.to_crs(meter_crs)
    bldgs_m = bldgs_ll.to_crs(meter_crs).reset_index(drop=True)
    bldgs_m["bldg_idx"] = bldgs_m.index

    # Keep polygon gyms as-is
    gym_poly = gyms_m[gyms_m.geometry.type.isin(["Polygon", "MultiPolygon"])].copy()
    gym_pts = gyms_m[gyms_m.geometry.type == "Point"].copy()

    indoor_geoms = []

    # Handle point gyms inside buildings
    if not gym_pts.empty:
        joined = gpd.sjoin(gym_pts, bldgs_m[["geometry", "bldg_idx"]], how="left", predicate="within")
        for _, row in joined.iterrows():
            if pd.isna(row["bldg_idx"]):
                continue  # Not inside any building
            pt = row.geometry
            bldg_geom = bldgs_m.loc[row["bldg_idx"]].geometry
            square = box(pt.x - INDOOR_HALF_M, pt.y - INDOOR_HALF_M,
                         pt.x + INDOOR_HALF_M, pt.y + INDOOR_HALF_M)
            indoor_zone = square.intersection(bldg_geom)
            if not indoor_zone.is_empty:
                indoor_geoms.append({
                    "geometry": indoor_zone,
                    "gym_name": row.get("name"),
                    "source": "osm_point_in_building"
                })

    # Add polygon gyms
    if not gym_poly.empty:
        for _, row in gym_poly.iterrows():
            indoor_geoms.append({
                "geometry": row.geometry,
                "gym_name": row.get("name"),
                "source": "osm_polygon"
            })

    if not indoor_geoms:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

    indoor_gdf = gpd.GeoDataFrame(indoor_geoms, crs=meter_crs)
    indoor_gdf = safe_clip(indoor_gdf, city_geom_m)
    indoor_gdf = indoor_gdf[~indoor_gdf.geometry.is_empty]
    return indoor_gdf.to_crs(equal_area_crs)

# ---------------- OSM Parks ----------------
def get_osm_parks_for_city(city_geom_proj):
    parks_ll = fetch_osm_by_tiling(parks_tags, city_geom_proj)
    if parks_ll.empty:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
    parks_poly = parks_ll[parks_ll.geometry.type.isin(["Polygon", "MultiPolygon"])]
    if parks_poly.empty:
        return gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
    parks_proj = parks_poly.to_crs(equal_area_crs)
    city_gdf = gpd.GeoDataFrame(geometry=[city_geom_proj], crs=equal_area_crs)
    clipped = safe_clip(parks_proj, city_gdf)
    return clipped[~clipped.geometry.is_empty]

# ---------------- Main per-state processing ----------------
def process_state(state_fips, cities_gdf):
    try:
        print(f"Processing state {state_fips}...")
        features = fetch_parkserve_features(state_fips)
        parks_ps = gpd.GeoDataFrame(
            geometry=[shape(f["geometry"]) for f in features if f.get("geometry")],
            crs=web_crs
        ).to_crs(equal_area_crs) if features else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

        cities_state = cities_gdf[cities_gdf["statefp"] == state_fips].copy()
        results = []

        for _, city in cities_state.iterrows():
            city_geom = city.geometry
            city_name = city.get("NAME", city.get("NAME10", "Unknown"))

            # ParkServe
            ps_clipped = safe_clip(parks_ps, city_geom) if not parks_ps.empty else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
            ps_union = unary_union(ps_clipped.geometry) if not ps_clipped.empty else None
            area_ps = safe_area(ps_union)

            # OSM Parks
            osm_gdf = get_osm_parks_for_city(city_geom)
            osm_union = unary_union(osm_gdf.geometry) if not osm_gdf.empty else None
            area_osm = safe_area(osm_union)

            # Indoor gyms
            indoor_gdf = load_osm_indoor_for_city(city_geom)
            indoor_union = unary_union(indoor_gdf.geometry) if not indoor_gdf.empty else None
            area_indoor = safe_area(indoor_union)

            # Unified
            components = [g for g in [ps_union, osm_union, indoor_union] if g is not None]
            unified_geom = unary_union(components) if components else None
            area_total = safe_area(unified_geom)

            results.append({
                "state_fips": state_fips,
                "placefp": city.get("PLACEFP"),
                "city_name": city_name,
                "geometry": unified_geom,
                "geometry_parkserve": ps_union,
                "geometry_osm": osm_union,
                "geometry_indoor": indoor_union,
                "area_m2_parks": area_total,
                "area_m2_parkserve": area_ps,
                "area_m2_osm": area_osm,
                "area_m2_indoor": area_indoor,
                "date": today,
            })

        if not results:
            print(f"No cities processed for state {state_fips}")
            return state_fips

        state_gdf = gpd.GeoDataFrame(results, crs=equal_area_crs)

        # Save
        gpkg_path = os.path.join(output_dir, f"state_{state_fips}_unified.gpkg")
        csv_path = os.path.join(output_dir, f"state_{state_fips}_unified.csv")

        # Unified layer
        unified = state_gdf.drop(columns=["geometry_parkserve", "geometry_osm", "geometry_indoor"], errors="ignore")
        unified = unified.set_geometry("geometry")
        unified.to_file(gpkg_path, layer="unified_park_area", driver="GPKG")

        # Source layers
        for name, col in [("parkserve_geom", "geometry_parkserve"),
                          ("osm_geom", "geometry_osm"),
                          ("indoor_gyms", "geometry_indoor")]:
            layer = state_gdf.drop(columns=[c for c in ["geometry", "geometry_parkserve", "geometry_osm", "geometry_indoor"] if c != col], errors="ignore")
            layer = layer.rename(columns={col: "geometry"}).set_geometry("geometry")
            layer = layer[~layer.geometry.is_empty]
            if layer.empty:
                layer = gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
            layer.to_file(gpkg_path, layer=name, driver="GPKG")

        # CSV
        state_gdf.drop(columns=["geometry", "geometry_parkserve", "geometry_osm", "geometry_indoor"], errors="ignore").to_csv(csv_path, index=False)

        print(f"State {state_fips} completed: {len(state_gdf)} places")
        return state_fips

    except Exception as e:
        print(f"Error in state {state_fips}: {e}")
        return None

def main():
    print("Loading counties and places...")
    counties_gdf = gpd.read_file(COUNTY_URL).to_crs(equal_area_crs)
    state_fips_list = sorted([s for s in counties_gdf["STATEFP"].unique() 
                              if s not in EXCLUDED_STATEFPS])

    if TEST_MODE:
        print(f"TEST MODE: running on {TEST_STATES}")
        state_fips_list = [s for s in state_fips_list if s in TEST_STATES]

    city_gdfs = []
    for fips in state_fips_list:
        url = PLACE_BASE_URL.format(statefp=fips.zfill(2))
        try:
            gdf = gpd.read_file(url)                  # ← read once
            gdf["statefp"] = fips                      # ← add state FIPS
            gdf = gdf.to_crs(equal_area_crs)           # ← reproject once
            city_gdfs.append(gdf)
            print(f"Loaded places for state {fips}")
        except Exception as e:
            print(f"Failed to load places for state {fips}: {e}")

    if not city_gdfs:
        raise RuntimeError("No place data loaded")

    cities_gdf = pd.concat(city_gdfs, ignore_index=True)
    print(f"All place boundaries loaded: {len(cities_gdf)} cities total")

    # ---------------- Parallel execution ----------------
    max_workers = 8  # Hellgate gives you 8 CPUs → use them all
    print(f"Starting processing of {len(state_fips_list)} states with {max_workers} workers...")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_state, fips, cities_gdf) 
                   for fips in state_fips_list]
        for future in as_completed(futures):
            result = future.result()
            if result:
                print(f"Completed state {result}")

    print("All done! Output in:", output_dir)


if __name__ == "__main__":
    main()