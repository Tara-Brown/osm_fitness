#!/usr/bin/env python3
"""
Resumable nationwide city park area pipeline
- Per-state ParkServe cache: cache/parkserve_{statefp}.geojson
- Per-city optional OSM cache: cache/osm/{statefp}_{placefp}.gpkg
- Per-state CSV checkpoint: output/city_parks_{statefp}_{date}.csv
"""

import os
import time
import json
import requests
import geopandas as gpd
import pandas as pd
import osmnx as ox
from shapely.geometry import shape
from shapely.ops import unary_union
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---------- CONFIG ----------
TODAY = datetime.utcnow().strftime("%Y-%m-%d")
OUTPUT_DIR = "output"
CACHE_DIR = "cache"
PARKSERVE_CACHE_DIR = os.path.join(CACHE_DIR, "parkserve")
OSM_CACHE_DIR = os.path.join(CACHE_DIR, "osm")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PARKSERVE_CACHE_DIR, exist_ok=True)
os.makedirs(OSM_CACHE_DIR, exist_ok=True)

# Data sources
PLACE_BASE_URL = "https://www2.census.gov/geo/tiger/TIGER2025/PLACE/tl_2025_{statefp}_place.zip"
PARKSERVE_URL = "https://server7.tplgis.org/arcgis7/rest/services/ParkServe/ParkServe_ProdNew/MapServer/2/query"
PARKSERVE_CONDITION = "(park_designation = 'LP' OR park_designation = 'LREC')"
OSM_TAGS = {"leisure": ["park", "pitch"]}

# CRS
WGS84 = "EPSG:4326"
EQUAL_AREA_CRS = "EPSG:5070"

# Parallelism & politeness
MAX_WORKERS_PER_STATE = 6        # per-state thread count for cities
MAX_STATES_AT_ONCE = 1           # process states sequentially by default for safety
SLEEP_BETWEEN_REQUESTS = 0.3
SLEEP_BETWEEN_STATES = 3.0

# Retry/backoff
REQUEST_RETRIES = 3
RETRY_BACKOFF = 2.0

# Caching toggles
OSM_CACHE = True                 # saves per-city OSM GPKG after fetch
PARKSERVE_CACHE = True           # saves per-state ParkServe GeoJSON

# OSMnx settings
ox.settings.use_cache = True
ox.settings.log_console = False
ox.settings.timeout = 180

# ---------- HELPERS ----------
def read_tiger_places_for_state(statefp: str) -> gpd.GeoDataFrame:
    """Download/load TIGER PLACE for a state; normalize statefp column and set WGS84 CRS."""
    url = PLACE_BASE_URL.format(statefp=str(statefp).zfill(2))
    gdf = gpd.read_file(url).to_crs(WGS84)
    gdf["statefp"] = str(statefp).zfill(2)
    return gdf

def parkserve_fetch_by_state(state_fips: str) -> gpd.GeoDataFrame:
    """
    Fetch ParkServe features for a state; uses numeric FIPS in WHERE clause.
    Caches to cache/parkserve_{statefp}.geojson if PARKSERVE_CACHE True.
    """
    cache_path = os.path.join(PARKSERVE_CACHE_DIR, f"parkserve_{state_fips}.geojson")
    if PARKSERVE_CACHE and os.path.exists(cache_path):
        try:
            g = gpd.read_file(cache_path)
            print(f"  Loaded cached ParkServe for {state_fips} ({len(g)} features)")
            return g
        except Exception:
            print(f"  Warning: failed to read cached ParkServe {cache_path}; refetching.")

    features = []
    offset = 0
    where_clause = f"Park_State_FIPS = {int(state_fips)} AND {PARKSERVE_CONDITION}"

    while True:
        params = {
            "where": where_clause,
            "outFields": "*",
            "f": "geojson",
            "resultOffset": offset,
            "resultRecordCount": 1000,
        }
        page = None
        for attempt in range(1, REQUEST_RETRIES + 1):
            try:
                r = requests.get(PARKSERVE_URL, params=params, timeout=60)
                r.raise_for_status()
                page = r.json().get("features", [])
                break
            except Exception as exc:
                wait = (RETRY_BACKOFF ** (attempt - 1)) + 0.1 * attempt
                print(f"    ParkServe request attempt {attempt} failed for state {state_fips}: {exc} — sleeping {wait:.1f}s")
                time.sleep(wait)
        if page is None:
            print(f"    ParkServe requests failed after retries for state {state_fips}; proceeding with what we have.")
            break
        if not page:
            break
        features.extend(page)
        if len(page) < 1000:
            break
        offset += 1000

    if not features:
        return gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)

    # Convert features to geometries (robust)
    geoms = []
    for f in features:
        geom = f.get("geometry")
        if not geom:
            continue
        try:
            geoms.append(shape(geom))
        except Exception:
            continue

    gdf = gpd.GeoDataFrame(geometry=geoms, crs=WGS84)
    # attempt to fix invalid geometries
    try:
        gdf["geometry"] = gdf.geometry.buffer(0)
    except Exception:
        pass
    gdf = gdf[~gdf.geometry.is_empty]

    if PARKSERVE_CACHE:
        try:
            gdf.to_file(cache_path, driver="GeoJSON")
            print(f"  Cached ParkServe to {cache_path} ({len(gdf)} features)")
        except Exception as exc:
            print(f"  Warning: failed to write ParkServe cache: {exc}")

    return gdf

def osm_cache_path(statefp: str, placefp: str) -> str:
    safe = f"{statefp}_{placefp}"
    return os.path.join(OSM_CACHE_DIR, f"osm_{safe}.gpkg")

def osm_fetch_by_city(city_geom):
    """Fetch OSM parks for a single city polygon."""
    for attempt in range(1, REQUEST_RETRIES + 1):
        try:
            gdf = ox.geometries.geometries_from_polygon(city_geom, tags=OSM_TAGS)
            if gdf is None:
                return gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)
            gdf["geometry"] = gdf.geometry.buffer(0)
            gdf = gdf[~gdf.geometry.is_empty]
            return gdf
        except Exception:
            sleep_time = RETRY_BACKOFF ** attempt
            print(f"OSM attempt {attempt} failed: {e}; sleeping {sleep_time}s")
            time.sleep(sleep_time)
    return gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)


def compute_areas(osm_gdf: gpd.GeoDataFrame, park_gdf: gpd.GeoDataFrame, city_geom):
    """
    Clip both layers to city polygon in WGS84, reproject clipped layers to equal-area,
    then compute individual, union, and intersection areas (m^2).
    """
    if (osm_gdf is None or osm_gdf.empty) and (park_gdf is None or park_gdf.empty):
        return 0.0, 0.0, 0.0, 0.0

    # Ensure WGS84 on inputs
    if osm_gdf is None:
        osm_gdf = gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)
    else:
        if osm_gdf.crs is None:
            osm_gdf = osm_gdf.set_crs(WGS84, allow_override=True)
        else:
            osm_gdf = osm_gdf.to_crs(WGS84)

    if park_gdf is None:
        park_gdf = gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)
    else:
        if park_gdf.crs is None:
            park_gdf = park_gdf.set_crs(WGS84, allow_override=True)
        else:
            park_gdf = park_gdf.to_crs(WGS84)

    # City geometry as GeoDataFrame in WGS84
    city_g = gpd.GeoDataFrame(geometry=[city_geom], crs=WGS84)

    # Clip in WGS84 (same CRS)
    try:
        osm_clip = gpd.clip(osm_gdf, city_g)
    except Exception:
        osm_clip = gpd.clip(osm_gdf.to_crs(WGS84), city_g)
    try:
        park_clip = gpd.clip(park_gdf, city_g)
    except Exception:
        park_clip = gpd.clip(park_gdf.to_crs(WGS84), city_g)

    # Fast exit if both clipped are empty
    if osm_clip.empty and park_clip.empty:
        return 0.0, 0.0, 0.0, 0.0

    # Reproject clipped to equal-area for area calc
    if not osm_clip.empty:
        osm_clip = osm_clip.to_crs(EQUAL_AREA_CRS)
        try:
            osm_union = unary_union(osm_clip.geometry)
            osm_area = osm_union.area
        except Exception:
            osm_area = osm_clip.geometry.area.sum()
            osm_union = None
    else:
        osm_area = 0.0
        osm_union = None

    if not park_clip.empty:
        park_clip = park_clip.to_crs(EQUAL_AREA_CRS)
        try:
            park_union = unary_union(park_clip.geometry)
            park_area = park_union.area
        except Exception:
            park_area = park_clip.geometry.area.sum()
            park_union = None
    else:
        park_area = 0.0
        park_union = None

    # Compute overlap & combined
    overlap_area = 0.0
    combined_area = 0.0
    if osm_union is not None and park_union is not None:
        try:
            overlap_area = osm_union.intersection(park_union).area
            combined_area = osm_area + park_area - overlap_area
        except Exception:
            # fallback to overlay if unary ops fail
            try:
                combined = gpd.overlay(osm_clip, park_clip, how="union")
                overlap = gpd.overlay(osm_clip, park_clip, how="intersection")
                combined_area = combined.geometry.area.sum()
                overlap_area = overlap.geometry.area.sum()
            except Exception:
                combined_area = osm_area + park_area
                overlap_area = 0.0
    else:
        combined_area = osm_area if osm_union is not None else park_area
        overlap_area = 0.0

    return float(osm_area), float(park_area), float(combined_area), float(overlap_area)

# ---------- Per-city processing ----------
def process_city(city_row, parkserve_state):
    """Compute areas for a single city geometry. parkserve_state is a GeoDataFrame (state-level)."""
    try:
        statefp = city_row["statefp"]
        placefp = str(city_row.get("PLACEFP", "")).zfill(5) if "PLACEFP" in city_row else ""
        # 1) fetch OSM parks (with caching)
        osm_gdf = osm_fetch_by_city(city_row.geometry, statefp=statefp, placefp=placefp)

        # 2) restrict parkserve_state to the city's bbox before compute (clip will further refine)
        if parkserve_state is None or parkserve_state.empty:
            park_gdf = gpd.GeoDataFrame(columns=["geometry"], crs=WGS84)
        else:
            # clip parkserve_state by city polygon to reduce size
            try:
                park_gdf = gpd.clip(parkserve_state, gpd.GeoDataFrame(geometry=[city_row.geometry], crs=WGS84))
            except Exception:
                park_gdf = parkserve_state  # fallback (will be clipped downstream)

        osm_area, park_area, combined_area, overlap_area = compute_areas(osm_gdf, park_gdf, city_row.geometry)

        return {
            "statefp": statefp,
            "placefp": placefp,
            "city_name": city_row.get("NAME", ""),
            "osm_area_m2": osm_area,
            "parkserve_area_m2": park_area,
            "combined_area_m2": combined_area,
            "overlap_area_m2": overlap_area,
            "date": TODAY,
        }
    except Exception as exc:
        print(f"    City processing error {city_row.get('NAME', 'unknown')}: {exc}")
        return {
            "statefp": city_row.get("statefp", ""),
            "placefp": city_row.get("PLACEFP", ""),
            "city_name": city_row.get("NAME", ""),
            "osm_area_m2": float("nan"),
            "parkserve_area_m2": float("nan"),
            "combined_area_m2": float("nan"),
            "overlap_area_m2": float("nan"),
            "date": TODAY,
        }

# ---------- State runner with resume & caches ----------
def run_state(statefp, places_all_gdf):
    out_csv = os.path.join(OUTPUT_DIR, f"city_parks_{statefp}_{TODAY}.csv")
    if os.path.exists(out_csv):
        print(f"  Skipping state {statefp} (already processed): {out_csv}")
        return pd.read_csv(out_csv)

    print(f"\n=== Processing state {statefp} ===")
    # load/prepare places for this state
    cities_state = places_all_gdf[places_all_gdf["statefp"] == statefp].reset_index(drop=True)
    if cities_state.empty:
        print(f"  No places found for state {statefp}")
        return pd.DataFrame()

    # Load or fetch ParkServe for this state (cached)
    parkserve_cache = os.path.join(PARKSERVE_CACHE_DIR, f"parkserve_{statefp}.geojson")
    if PARKSERVE_CACHE and os.path.exists(parkserve_cache):
        try:
            parkserve_state = gpd.read_file(parkserve_cache)
            print(f"  Loaded ParkServe cache for {statefp} ({len(parkserve_state)} features)")
        except Exception:
            parkserve_state = parkserve_fetch_by_state(statefp)
    else:
        parkserve_state = parkserve_fetch_by_state(statefp)
        # save done in fetch function if enabled

    # trim columns to minimize memory
    if parkserve_state is not None and not parkserve_state.empty:
        parkserve_state = parkserve_state[["geometry"]].copy()

    results = []
    with ThreadPoolExecutor(max_workers=min(MAX_WORKERS_PER_STATE, max(1, os.cpu_count() or 2))) as exe:
        futures = {exe.submit(process_city, row, parkserve_state): idx for idx, row in cities_state.iterrows()}
        for i, fut in enumerate(as_completed(futures)):
            res = fut.result()
            results.append(res)
            if i % 10 == 0:
                time.sleep(SLEEP_BETWEEN_REQUESTS)

    state_df = pd.DataFrame(results)
    # save per-state CSV checkpoint
    state_df.to_csv(out_csv, index=False)
    print(f"  Saved state CSV: {out_csv} ({len(state_df)} rows)")
    return state_df

# ---------- Nationwide runner ----------
def run_nationwide():
    print("=== Starting nationwide run (resumable) ===")
    # build state list from TIGER place files (download on-demand)
    state_list = [f"{i:02d}" for i in range(1, 57)]
    # load all places in a memory-friendly way (one state at a time to build a combined df)
    places_all = []
    for st in state_list:
        try:
            g = read_tiger_places_for_state(st)
            places_all.append(g)
            print(f"Loaded {len(g)} places for state {st}")
        except Exception as exc:
            print(f"  Failed to load places for {st}: {exc}")

    if not places_all:
        raise RuntimeError("No place shapefiles loaded; aborting.")

    places_all_gdf = pd.concat(places_all, ignore_index=True)

    # make sure statefp exists as strings zero-padded
    if "statefp" not in places_all_gdf.columns:
        raise RuntimeError("places_all missing 'statefp' column")

    # process states sequentially (safe); parallelize per state if desired
    processed_states = []
    for statefp in sorted(places_all_gdf["statefp"].unique()):
        try:
            df_state = run_state(statefp, places_all_gdf)
            if df_state is not None and not df_state.empty:
                processed_states.append(df_state)
        except Exception as exc:
            print(f"State {statefp} failed: {exc}")
        time.sleep(SLEEP_BETWEEN_STATES)

    if processed_states:
        combined = pd.concat(processed_states, ignore_index=True)
        combined_path = os.path.join(OUTPUT_DIR, f"city_parks_nationwide_{TODAY}.csv")
        combined.to_csv(combined_path, index=False)
        print(f"\n✅ Nationwide complete. Saved {combined_path} ({len(combined)} rows)")
    else:
        print("\nNo state outputs were produced.")

# ---------- MAIN ----------
if __name__ == "__main__":
    run_nationwide()
