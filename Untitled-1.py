import os
import time
import requests
import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import osmnx as ox

# ---------- CONFIG ----------
COUNTY_URL = "https://www2.census.gov/geo/tiger/TIGER2025/COUNTY/tl_2025_us_county.zip"
PLACE_BASE_URL = "https://www2.census.gov/geo/tiger/TIGER2025/PLACE/tl_2025_{statefp}_place.zip"
PARKSERVE_URL = "https://server7.tplgis.org/arcgis7/rest/services/ParkServe/ParkServe_ProdNew/MapServer/2/query"
PARKSERVE_CONDITION = "(park_designation = 'LP' OR park_designation = 'LREC')"
OSM_TAGS = {"leisure": ["park", "pitch"]}
WGS84 = "EPSG:4326"
EQUAL_AREA_CRS = "EPSG:5070"

MAX_WORKERS_PER_STATE = 6
SLEEP_BETWEEN_REQUESTS = 0.5
SLEEP_BETWEEN_STATES = 6.0
REQUEST_RETRIES = 3
RETRY_BACKOFF = 2.0

OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)
today = datetime.utcnow().strftime("%Y-%m-%d")

# OSMnx settings
ox.settings.use_cache = True
ox.settings.log_console = False
ox.settings.timeout = 180

# ---------- Helpers ----------
def read_tiger_counties():
    gdf = gpd.read_file(COUNTY_URL).to_crs(WGS84)
    gdf["statefp"] = gdf["STATEFP"]
    gdf["countyfp"] = gdf["COUNTYFP"]
    gdf["geoid"] = gdf["GEOID"]
    gdf["county_name"] = gdf["NAME"]
    return gdf

def read_tiger_places_for_state(statefp):
    url = PLACE_BASE_URL.format(statefp=str(statefp).zfill(2))
    gdf = gpd.read_file(url).to_crs(WGS84)
    gdf["statefp"] = str(statefp).zfill(2)
    return gdf

def parkserve_fetch_by_county(county_name, state_fips):
    features = []
    offset = 0
    while True:
        where_clause = f"Park_State_FIPS = '{state_fips}' AND Park_County = '{county_name} County' AND {PARKSERVE_CONDITION}"
        params = {"where": where_clause, "outFields": "*", "f": "geojson", "resultOffset": offset, "resultRecordCount": 1000}
        for attempt in range(1, REQUEST_RETRIES + 1):
            try:
                r = requests.get(PARKSERVE_URL, params=params, timeout=60)
                r.raise_for_status()
                page = r.json().get("features", [])
                break
            except Exception:
                time.sleep(RETRY_BACKOFF ** attempt)
        if not page:
            break
        features.extend(page)
        if len(page) < 1000:
            break
        offset += 1000
    if not features:
        return gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)
    parks = gpd.GeoDataFrame.from_features(features, crs=WGS84)
    if "geometry" in parks:
        parks["geometry"] = parks.geometry.buffer(0)
        parks = parks[~parks.geometry.is_empty]
    return parks

def osm_fetch_by_county(county_geom):
    for attempt in range(1, REQUEST_RETRIES + 1):
        try:
            gdf = ox.features.features_from_polygon(county_geom, tags=OSM_TAGS)
            if gdf is None:
                return gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)
            gdf["geometry"] = gdf.geometry.buffer(0)
            gdf = gdf[~gdf.geometry.is_empty]
            return gdf
        except Exception:
            time.sleep(RETRY_BACKOFF ** attempt)
    return gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)

def compute_areas_only(osm_gdf, park_gdf):
    """Compute areas in equal-area CRS, no geometry saved"""
    osm_area = 0.0
    park_area = 0.0
    combined_area = 0.0
    overlap_area = 0.0

    if not osm_gdf.empty:
        osm_proj = osm_gdf.to_crs(EQUAL_AREA_CRS)
        osm_union = unary_union(osm_proj.geometry)
        osm_area = osm_union.area
    else:
        osm_union = None

    if not park_gdf.empty:
        park_proj = park_gdf.to_crs(EQUAL_AREA_CRS)
        park_union = unary_union(park_proj.geometry)
        park_area = park_union.area
    else:
        park_union = None

    if osm_union and park_union:
        combined_union = unary_union([osm_union, park_union])
        combined_area = combined_union.area
        overlap_area = osm_area + park_area - combined_area
    elif osm_union:
        combined_area = osm_area
    elif park_union:
        combined_area = park_area

    return osm_area, park_area, combined_area, overlap_area

def process_county(county_row, cities_gdf_state):
    try:
        cities_in_county = cities_gdf_state[cities_gdf_state.intersects(county_row.geometry)]
        if cities_in_county.empty:
            return pd.DataFrame([{
                "geoid": county_row["geoid"],
                "statefp": county_row["statefp"],
                "countyfp": county_row["countyfp"],
                "county_name": county_row["county_name"],
                "osm_area_m2": 0.0,
                "parkserve_area_m2": 0.0,
                "combined_area_m2": 0.0,
                "overlap_area_m2": 0.0,
                "date": today
            }])

        osm_gdf = osm_fetch_by_county(county_row.geometry)
        park_gdf = parkserve_fetch_by_county(county_row["county_name"], county_row["statefp"])

        if not osm_gdf.empty:
            osm_gdf = gpd.clip(osm_gdf, cities_in_county)
        if not park_gdf.empty:
            park_gdf = gpd.clip(park_gdf, cities_in_county)

        osm_area, park_area, combined_area, overlap_area = compute_areas_only(osm_gdf, park_gdf)

        return pd.DataFrame([{
            "geoid": county_row["geoid"],
            "statefp": county_row["statefp"],
            "countyfp": county_row["countyfp"],
            "county_name": county_row["county_name"],
            "osm_area_m2": osm_area,
            "parkserve_area_m2": park_area,
            "combined_area_m2": combined_area,
            "overlap_area_m2": overlap_area,
            "date": today
        }])
    except Exception as e:
        print(f"ERROR processing county {county_row['county_name']}: {e}")
        return pd.DataFrame([{
            "geoid": county_row["geoid"],
            "statefp": county_row["statefp"],
            "countyfp": county_row["countyfp"],
            "county_name": county_row["county_name"],
            "osm_area_m2": float("nan"),
            "parkserve_area_m2": float("nan"),
            "combined_area_m2": float("nan"),
            "overlap_area_m2": float("nan"),
            "date": today
        }])

# ---------- State runner ----------
def run_state(statefp, counties_gdf, cities_gdf_all, max_workers=MAX_WORKERS_PER_STATE):
    print(f"\n=== PROCESSING STATE {statefp} ===")
    state_counties = counties_gdf[counties_gdf["statefp"] == statefp].reset_index(drop=True)
    cities_state = cities_gdf_all[cities_gdf_all["STATEFP"] == statefp]

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as exe:
        futures = {exe.submit(process_county, row, cities_state): idx for idx, row in state_counties.iterrows()}
        for i, fut in enumerate(as_completed(futures)):
            results.append(fut.result())
            if i % 10 == 0:
                time.sleep(SLEEP_BETWEEN_REQUESTS)

    if not results:
        return None
    return pd.concat(results, ignore_index=True)

# ---------- Nationwide runner ----------
def run_nationwide(counties_gdf, cities_gdf_all):
    all_states = sorted(counties_gdf["statefp"].unique())
    per_state_results = []
    for statefp in all_states:
        try:
            state_df = run_state(statefp, counties_gdf, cities_gdf_all)
            if state_df is not None:
                per_state_results.append(state_df)
        except Exception as e:
            print(f"State {statefp} failed: {e}")
        time.sleep(SLEEP_BETWEEN_STATES)

    final_df = pd.concat(per_state_results, ignore_index=True)
    final_csv = os.path.join(OUTPUT_DIR, "county_parks_all.csv")
    final_df.to_csv(final_csv, index=False)
    print(f"\nFINISHED: wrote {len(final_df)} counties -> {final_csv}")
    return final_df

# ---------- Main ----------
if __name__ == "__main__":
    print("Loading TIGER county and place data...")
    counties = read_tiger_counties()
    state_list = sorted(counties["statefp"].unique())
    place_gdfs = []
    for st in state_list:
        try:
            g = read_tiger_places_for_state(st)
            place_gdfs.append(g)
            print(f"Loaded {len(g)} places for state {st}")
        except Exception as e:
            print(f"Failed to load places for {st}: {e}")
    places_all = pd.concat(place_gdfs, ignore_index=True)
    counties = counties.to_crs(WGS84)

    final_df = run_nationwide(counties, places_all)
