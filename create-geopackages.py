import os
import requests
import geopandas as gpd
from shapely.geometry import shape
from shapely.ops import unary_union
import osmnx as ox
import pandas as pd
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import math
from shapely.validation import make_valid
import warnings
warnings.filterwarnings("ignore", message="invalid value encountered in area", category=RuntimeWarning)


# --- Config ---
equal_area_crs = "EPSG:5070"  # Equal-area projection for area accuracy
today = datetime.today().strftime("%Y-%m-%d")
output_dir = "unified_geopackages_"
os.makedirs(output_dir, exist_ok=True)

# --- TEST MODE ---
TEST_MODE = False       # Change to False for full run
TEST_STATES = ["30", "53"]  # Example: Montana (30), Washington (53)

# --- Data Sources ---
COUNTY_URL = "https://www2.census.gov/geo/tiger/TIGER2025/COUNTY/tl_2025_us_county.zip"
PLACE_BASE_URL = "https://www2.census.gov/geo/tiger/TIGER2025/PLACE/tl_2025_{statefp}_place.zip"

PARKSERVE_URL = "https://server7.tplgis.org/arcgis7/rest/services/ParkServe/ParkServe_ProdNew/MapServer/2/query"
PARKSERVE_CONDITION = "(park_designation = 'LP' OR park_designation = 'LREC')"

tags = {"leisure": ["park", "pitch", "sports_centre"]}


# --- Helper Functions ---
def fetch_parkserve_features(state_fips):
    """Fetch ParkServe features for a state FIPS."""
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
    return features

def safe_area(geom):
    """Return valid area (m²) or 0.0 if geometry is invalid/empty."""
    if geom is None or geom.is_empty:
        return 0.0
    try:
        geom_valid = make_valid(geom)
        a = geom_valid.area
        return 0.0 if math.isnan(a) else a
    except Exception:
        return 0.0


def get_osm_parks(city_geom):
    """Download OSM parks clipped to a city boundary."""
    try:
        # Convert to lat/lon for OSM query
        geom_ll = gpd.GeoSeries([city_geom], crs=equal_area_crs).to_crs("EPSG:4326").iloc[0]

        # Download OSM parks within the city
        gdf_osm = ox.features_from_polygon(geom_ll, tags=tags)

        # Filter to park-like polygons
        gdf_osm = gdf_osm[gdf_osm.geometry.type.isin(["Polygon", "MultiPolygon"])]

        # Clip and ensure CRS matches output
        gdf_osm = gdf_osm.clip(geom_ll)
        gdf_osm = gdf_osm.set_crs("EPSG:4326")

        return gdf_osm

    except Exception as e:
        print(f"OSM error: {e}")
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")



def process_state(state_fips):
    """Process one state: ParkServe + OSM → unified GeoPackage (3 layers)."""
    try:
        print(f"Processing state {state_fips}...")

        # --- Load ParkServe parks ---
        features = fetch_parkserve_features(state_fips)
        parks = (
            gpd.GeoDataFrame(
                geometry=[shape(f["geometry"]) for f in features if f.get("geometry")],
                crs="EPSG:4326",
            ).to_crs(equal_area_crs)
            if features else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
        )
        parks["geometry"] = parks.buffer(0)
        parks = parks[~parks.geometry.is_empty]

        # --- Cities for this state ---
        cities_state = cities_gdf[cities_gdf["statefp"] == state_fips]
        results = []

        # --- Process each city ---
        for _, city in cities_state.iterrows():
            city_geom = city.geometry
            city_name = city["NAME"]

            # --- ParkServe parks clipped to city ---
            ps_city = (
                gpd.overlay(parks, gpd.GeoDataFrame(geometry=[city_geom], crs=equal_area_crs), how="intersection")
                if not parks.empty
                else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
            )

            if not ps_city.empty:
                parkserve_union = unary_union(ps_city.geometry)
                area_m2_parkserve = safe_area(parkserve_union)
            else:
                parkserve_union = None
                area_m2_parkserve = 0.0

            # --- OSM parks clipped to city ---
            osm_raw = get_osm_parks(city_geom)
            osm_city = osm_raw.to_crs(equal_area_crs) if not osm_raw.empty else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)

            if not osm_city.empty:
                osm_union = unary_union(osm_city.geometry)
                area_m2_osm = safe_area(osm_union)
            else:
                osm_union = None
                area_m2_osm = 0.0

            # --- Combined (unified) parks ---
            combined = pd.concat([ps_city, osm_city], ignore_index=True) if not ps_city.empty or not osm_city.empty else gpd.GeoDataFrame(geometry=[], crs=equal_area_crs)
            if not combined.empty:
                unified_geom = unary_union(combined.geometry)
                area_m2_parks = safe_area(unified_geom)
            else:
                unified_geom = None
                area_m2_parks = 0.0

            # --- Append results for this city ---
            results.append({
                "state_fips": state_fips,
                "placefp": city["PLACEFP"],
                "city_name": city_name,
                "geometry": unified_geom,            # unified
                "geometry_parkserve": parkserve_union,  # ParkServe only
                "geometry_osm": osm_union,              # OSM only
                "area_m2_parks": area_m2_parks,
                "area_m2_parkserve": area_m2_parkserve,
                "area_m2_osm": area_m2_osm,
                "date": today,
            })

        # --- Convert to GeoDataFrame ---
        state_gdf = gpd.GeoDataFrame(results, crs=equal_area_crs)

        # --- Save outputs ---
        gpkg_path = os.path.join(output_dir, f"state_{state_fips}_unified.gpkg")
        csv_path = os.path.join(output_dir, f"state_{state_fips}_unified.csv")

        # 1️⃣ Unified layer
        unified_layer = state_gdf.drop(columns=["geometry_parkserve", "geometry_osm"])
        unified_layer.set_crs(equal_area_crs, allow_override=True).to_file(
            gpkg_path, driver="GPKG", layer="unified_park_area"
        )

        # 2️⃣ ParkServe-only layer
        ps_layer = state_gdf.drop(columns=["geometry", "geometry_osm"]).rename(columns={"geometry_parkserve": "geometry"})
        ps_layer = ps_layer.set_geometry("geometry").set_crs(equal_area_crs, allow_override=True)
        ps_layer.to_file(gpkg_path, driver="GPKG", layer="parkserve_geom")

        # 3️⃣ OSM-only layer
        osm_layer = state_gdf.drop(columns=["geometry", "geometry_parkserve"]).rename(columns={"geometry_osm": "geometry"})
        osm_layer = osm_layer.set_geometry("geometry").set_crs(equal_area_crs, allow_override=True)
        osm_layer.to_file(gpkg_path, driver="GPKG", layer="osm_geom")

        # 4️⃣ CSV (no geometry)
        state_gdf.drop(columns=["geometry", "geometry_parkserve", "geometry_osm"]).to_csv(csv_path, index=False)

        print(f" Finished state {state_fips}: {len(state_gdf)} cities saved.")
        return state_fips

    except Exception as e:
        print(f" Error in state {state_fips}: {e}")
        return None



# --- Load counties and cities once globally ---
print("Loading counties and cities...")
counties_gdf = gpd.read_file(COUNTY_URL).to_crs(equal_area_crs)
state_fips_list = sorted(counties_gdf["STATEFP"].unique().tolist())

# Apply test mode subset
if TEST_MODE:
    print(f" TEST MODE ACTIVE — processing {TEST_STATES}")
    state_fips_list = TEST_STATES

# Load city shapefiles for all (or test) states
city_gdfs = []
for fips in state_fips_list:
    url = PLACE_BASE_URL.format(statefp=fips.zfill(2))
    try:
        gdf = gpd.read_file(url)
        gdf["statefp"] = fips
        gdf = gdf.to_crs(equal_area_crs)
        city_gdfs.append(gdf)
    except Exception as e:
        print(f"Failed to load {url}: {e}")
cities_gdf = pd.concat(city_gdfs, ignore_index=True)


# --- Parallel Execution ---
if __name__ == "__main__":
    max_workers = 6 # adjust to your CPU (4–8 ideal)
    print(f"\n🚀 Starting parallel processing with {max_workers} workers...")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_state, fips) for fips in state_fips_list]
        for f in as_completed(futures):
            result = f.result()
            if result:
                print(f"🎯 Completed {result}")

    print("\n✅ All states processed and saved.")

