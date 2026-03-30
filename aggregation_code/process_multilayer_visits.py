#!/usr/bin/env python3
"""
Process cellphone mobility TSV files and match visits to multiple GeoPackage layers.
Uses city-level filtering for efficient spatial joins.

"""
import argparse
import sys
from pathlib import Path
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from datetime import datetime
import warnings
import re
import json
from collections import defaultdict
import pyarrow
import shapely

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Column names in TSV files  
TIME_COL = 'Unix Timestamp of Visit'
LON_COL = 'Lon of Visit'
LAT_COL = 'Lat of Visit'
DEVICE_COL = 'Hashed Device ID'

# GeoPackage layers to process
NEW_LAYERS = [
    'parks',
    'sports_areas',
]

# Processing parameters
CHUNK_SIZE = 100000
CRS = 'EPSG:4326'

STATES_DIR = Path("/mnt/beegfs/projects/cellphone_mobility/states")

# Different files use different abbreviations, unfortunately
# State mappings
STATE_ABBREV = {
    'AL': 'Alabama', 'AK': 'Alaska', 'AZ': 'Arizona', 'AR': 'Arkansas',
    'CA': 'California', 'CO': 'Colorado', 'CT': 'Connecticut', 'DE': 'Delaware',
    'FL': 'Florida', 'GA': 'Georgia', 'HI': 'Hawaii', 'ID': 'Idaho',
    'IL': 'Illinois', 'IN': 'Indiana', 'IA': 'Iowa', 'KS': 'Kansas',
    'KY': 'Kentucky', 'LA': 'Louisiana', 'ME': 'Maine', 'MD': 'Maryland',
    'MA': 'Massachusetts', 'MI': 'Michigan', 'MN': 'Minnesota', 'MS': 'Mississippi',
    'MO': 'Missouri', 'MT': 'Montana', 'NE': 'Nebraska', 'NV': 'Nevada',
    'NH': 'New Hampshire', 'NJ': 'New Jersey', 'NM': 'New Mexico', 'NY': 'New York',
    'NC': 'North Carolina', 'ND': 'North Dakota', 'OH': 'Ohio', 'OK': 'Oklahoma',
    'OR': 'Oregon', 'PA': 'Pennsylvania', 'RI': 'Rhode Island', 'SC': 'South Carolina',
    'SD': 'South Dakota', 'TN': 'Tennesee', 'TX': 'Texas', 'UT': 'Utah',
    'VT': 'Vermont', 'VA': 'Virginia', 'WA': 'Washington', 'WV': 'West Virginia',
    'WI': 'Wisconsin', 'WY': 'Wyoming', 'DC': 'District of Columbia'
}

STATE_FIPS = {
    'AL': '01', 'AK': '02', 'AZ': '04', 'AR': '05', 'CA': '06', 'CO': '08',
    'CT': '09', 'DE': '10', 'FL': '12', 'GA': '13', 'HI': '15', 'ID': '16',
    'IL': '17', 'IN': '18', 'IA': '19', 'KS': '20', 'KY': '21', 'LA': '22',
    'ME': '23', 'MD': '24', 'MA': '25', 'MI': '26', 'MN': '27', 'MS': '28',
    'MO': '29', 'MT': '30', 'NE': '31', 'NV': '32', 'NH': '33', 'NJ': '34',
    'NM': '35', 'NY': '36', 'NC': '37', 'ND': '38', 'OH': '39', 'OK': '40',
    'OR': '41', 'PA': '42', 'RI': '44', 'SC': '45', 'SD': '46', 'TN': '47',
    'TX': '48', 'UT': '49', 'VT': '50', 'VA': '51', 'WA': '53', 'WV': '54',
    'WI': '55', 'WY': '56', 'DC': '11'
}

# =============================================================================
# UTILITY FUNCTIONS (from original script)
# =============================================================================

import re
from pathlib import Path

LOOKUP_DIR = Path('/mnt/beegfs/projects/cellphone_mobility/testing_files/lookups')
GPKG_DIR = Path("/mnt/beegfs/hellgate/home/tb208541/osm_fitness/unified_geopackages_polygon_level")

def get_paths_from_filename(filename):
    """
    Given a TSV filename, return:
      - state_abbrev (e.g. 'OR')
      - state_name   (e.g. 'Oregon')
      - gpkg_path    (Path to state geopackage)
      - lookup_path  (Path to ubermedia lookup parquet)
    """
    city_name, state_name = parse_tsv_filename(filename)
    
    # Get abbrev + FIPS
    state_abbrev = next(k for k, v in STATE_ABBREV.items() if v == state_name)
    fips = STATE_FIPS[state_abbrev]
    
    # --- Geopackage ---
    gpkg_path = GPKG_DIR / f"state_{fips}_components.gpkg"
    
    # --- Ubermedia lookup ---
    # Build all candidate strings this state might appear as in filenames
    candidates = {
        state_name.lower(),                    # oregon
        state_name.replace(' ', '_').lower(),  # new_mexico
        state_abbrev.upper(),                  # OR
        f"_{state_abbrev.upper()}",            # _OR  (leading underscore variant)
    }
    
    lookup_path = None
    for f in LOOKUP_DIR.glob("ubermedia_lookup_*.parquet"):
        # Extract the part after 'ubermedia_lookup_'
        suffix = f.stem.replace("ubermedia_lookup_", "").strip("_")
        if suffix.lower() in {c.lower().strip("_") for c in candidates}:
            lookup_path = f
            break
    
    if lookup_path is None:
        raise FileNotFoundError(
            f"No ubermedia lookup found for {state_name} ({state_abbrev}) in {LOOKUP_DIR}\n"
            f"  Tried candidates: {candidates}"
        )
    
    return state_abbrev, state_name, gpkg_path, lookup_path

def normalize_city_name(name):
    """Remove double-underscore sections, strip numbers, replace underscores with spaces, lowercase."""
    # Remove everything within double underscores
    name = re.sub(r'__.*?__', '', name)
    # Remove trailing/leading underscores left behind
    name = name.strip('_')
    # Remove trailing numbers (e.g. Los_Angeles_5 -> Los_Angeles)
    name = re.sub(r'_\d+$', '', name)
    # Replace underscores with spaces
    name = name.replace('_', ' ')
    # Lowercase and strip
    return name.strip().lower()

def parse_tsv_filename(filename):
    if isinstance(filename, Path):
        stem = filename.stem
    else:
        stem = Path(filename).stem
    
    stem = re.sub(r'_pin_report(_\d{3})?$', '', stem)
    parts = stem.split('_')
    
    if len(parts) < 3:
        raise ValueError(f"Cannot parse filename: {filename}")
    
    parts = parts[1:]  # Skip leading numeric ID

    # Build lookup: normalized state name -> abbrev
    # Handles both 'California' and 'New_Mexico' style
    full_name_to_abbrev = {v.replace(' ', '_').lower(): k for k, v in STATE_ABBREV.items()}
    full_name_to_abbrev.update({v.lower(): k for k, v in STATE_ABBREV.items()})

    state_abbrev = None
    state_parts_count = 0

    # Try matching 1, 2, or 3 trailing parts as a state name (handles "New_Mexico", "Rhode_Island", etc.)
    for n in range(3, 0, -1):
        candidate = '_'.join(parts[-n:]).lower()
        if candidate in full_name_to_abbrev:
            state_abbrev = full_name_to_abbrev[candidate]
            state_parts_count = n
            break
        # Legacy: single uppercase abbreviation (WA, CA, etc.)
        if n == 1 and parts[-1].upper() in STATE_ABBREV:
            state_abbrev = parts[-1].upper()
            state_parts_count = 1
            break

    if not state_abbrev:
        raise ValueError(f"Cannot identify state in filename: {filename}")

    state = STATE_ABBREV[state_abbrev]
    city_name = '_'.join(parts[:-state_parts_count])  # e.g. 'Albuquerque_1_b'
    
    return city_name, state


def find_tsv_parts(tsv_path):
    """Find all parts of a multi-part TSV file."""
    parent = tsv_path.parent
    base_name = re.sub(r'_pin_report(_\d{3})?\.tsv$', '', tsv_path.name)
    
    if tsv_path.exists():
        parts = [tsv_path]
    else:
        parts = []
    
    # Look for numbered parts
    pattern = f"{base_name}_pin_report_*.tsv"
    numbered_parts = sorted(parent.glob(pattern))
    
    if numbered_parts:
        return numbered_parts
    
    return parts if parts else [tsv_path]


def find_geojson(city_name, state):
    """
    Find the geojson file for a city by matching filename.
    """
    state_dir_name = state.replace(' ', '_')
    #remove duplicates
    city_name = re.sub(r'__.*?__', '_', city_name) 
    city_name = re.sub(r'__','_',city_name)
    # Look in state directories
    for state_dir in STATES_DIR.glob(f"{state_dir_name}_*"):
        # Try 1: Direct match with underscores (from TSV parsing)
        geojson_path = state_dir / f"{city_name}.geojson"
        if geojson_path.exists():
            return geojson_path
        
        # Try 2: Replace underscores with spaces (for filenames like 'South Hill.geojson')
        name_with_spaces = city_name.replace('_', ' ')
        alt_path = state_dir / f"{name_with_spaces}.geojson"
        if alt_path.exists():
            return alt_path
        
        # Try 3: Replace spaces with underscores (reverse)
        name_with_underscores = city_name.replace(' ', '_')
        alt_path2 = state_dir / f"{name_with_underscores}.geojson"
        if alt_path2.exists():
            return alt_path2
        
        # Try 4: Name with numbers (some city parts do not have an initial numbering)
        alt_path3 = state_dir / f"{city_name}_1.geojson"
        if alt_path3.exists():
            return alt_path3
        
        # Try 5: Strip trailing _a or _b variant suffix
        name_no_variant = re.sub(r'_[ab]$', '', city_name)
        if name_no_variant != city_name:
            for variant in [name_no_variant, name_no_variant.replace('_', ' ')]:
                alt_path = state_dir / f"{variant}.geojson"
                if alt_path.exists():
                    return alt_path
        # Try 6: spaces in name part but underscore before number retained
        name_space_then_number = re.sub(r'^(.+)_(\d+)$', lambda m: m.group(1).replace('_', ' ') + '_' + m.group(2), name_no_variant)
        # San_Diego_6 -> San Diego_6
        alt_path = state_dir / f"{name_space_then_number}.geojson"
        if alt_path.exists():
            return alt_path
    # If no direct match found, search inside geojson files
    normalized_city = normalize_city_name(city_name)
    for state_dir in STATES_DIR.glob(f"{state_dir_name}_*"):
        for geojson_path in state_dir.glob("*.geojson"):
            try:
                with open(geojson_path) as f:
                    data = json.load(f)
                
                for feature in data.get('features', []):
                    props = feature.get('properties', {})
                    geojson_city = props.get('city') or props.get('name') or props.get('NAME')
                    if geojson_city and normalize_city_name(geojson_city) == normalized_city:
                        return geojson_path
            except Exception as e:
                continue
    
    print(f"  ⚠️  Geojson not found for city '{city_name}' in {state}")
    print(f"     Searched in: {STATES_DIR}/{state_dir_name}_*")
    return None

def load_geojson_cities(geojson_path):
    """Load list of cities from geojson file."""
    try:
        with open(geojson_path) as f:
            data = json.load(f)
        
        cities = []
        for feature in data.get('features', []):
            props = feature.get('properties', {})
            
            # Check for all_cities field (comma-separated list)
            if 'all_cities' in props:
                all_cities_str = props['all_cities']
                city_list = [normalize_city_name(c.strip()) for c in all_cities_str.split(',')]
                cities.extend(city_list)
        
        return cities
    except Exception as e:
        print(f"  Error loading geojson: {e}")
        return []


def load_all_city_polygons(state, cities, layers, gpkg_path=None):
    """
    Load polygons for all city/layer combinations.
    Returns: {(city, layer): GeoDataFrame}
    feature_id is set to gdf.index (fid) for stable cross-reload identity.
    """
    if gpkg_path is None:
        fips = STATE_FIPS.get(state)
        if not fips:
            print(f"  No FIPS code for {state}")
            return {}
        gpkg_path = GPKG_DIR / f"state_{fips}_components.gpkg"

    if not gpkg_path.exists():
        print(f"  Geopackage not found: {gpkg_path}")
        return {}

    city_layer_map = {}

    is_wyoming = False
    if state:
        is_wyoming = state.strip().lower() in ("wyoming", "wy")
    elif gpkg_path:
        is_wyoming = "state_56" in gpkg_path.stem

    for layer in layers:
        try:
            gdf = gpd.read_file(gpkg_path, layer=layer)

            if 'city_name' not in gdf.columns:
                print(f"  ⚠️  Layer {layer} has no 'city_name' column, skipping")
                continue

            # Use fid (gdf.index) as feature_id — stable across reloads
            gdf['feature_id'] = gdf.index

            if is_wyoming:
                all_cities = gdf['city_name'].dropna().unique()
                print(f"  Wyoming detected — loading all {len(all_cities)} cities from {layer}")
                for city_name in all_cities:
                    city_gdf = gdf[gdf['city_name'] == city_name].copy()
                    if len(city_gdf) > 0:
                        city_layer_map[(city_name, layer)] = city_gdf
                        print(f"    Found {len(city_gdf)} features for {city_name} in {layer}")
                continue

            for city in cities:
                city_base = re.sub(r'_\d+$', '', city)
                city_no_variant = re.sub(r'_[ab]$', '', city)
                city_base2 = re.sub(r'_\d+$', '', city_no_variant)

                city_gdf = gdf[
                    (gdf['city_name'].apply(normalize_city_name) == normalize_city_name(city)) |
                    (gdf['city_name'].apply(normalize_city_name) == normalize_city_name(city_base)) |
                    (gdf['city_name'].apply(normalize_city_name) == normalize_city_name(city_no_variant)) |
                    (gdf['city_name'].apply(normalize_city_name) == normalize_city_name(city_base2))
                ].copy()

                if len(city_gdf) > 0:
                    city_layer_map[(city, layer)] = city_gdf
                    print(f"    Found {len(city_gdf)} features for {city} in {layer}")

        except Exception as e:
            print(f"  ⚠️  Error loading layer {layer}: {e}")
            continue

    return city_layer_map


def filter_city_from_gdf(gdf, city_code):
    """
    Replicate the city filter from load_all_city_polygons for use in
    parquet rewriting. Assumes gdf.index is fid and feature_id = gdf.index.
    Returns (filtered_gdf, mask).
    """
    city_base = re.sub(r'_\d+$', '', city_code)
    city_no_variant = re.sub(r'_[ab]$', '', city_code)
    city_base2 = re.sub(r'_\d+$', '', city_no_variant)

    targets = {
        normalize_city_name(city_code),
        normalize_city_name(city_base),
        normalize_city_name(city_no_variant),
        normalize_city_name(city_base2),
    }

    mask = gdf['city_name'].apply(normalize_city_name).isin(targets)
    filtered = gdf[mask].copy()
    filtered['feature_id'] = filtered.index  # fid
    return filtered, mask


def process_tsv_chunk(chunk):
    """Convert TSV chunk to GeoDataFrame."""
    # Parse timestamp - Unix timestamp format
    chunk['datetime'] = pd.to_datetime(chunk[TIME_COL], unit='s', errors='coerce')
    chunk['date'] = chunk['datetime'].dt.date
    
    # Create geometry
    geometry = shapely.points(chunk[LON_COL].to_numpy(), chunk[LAT_COL].to_numpy())
    gdf = gpd.GeoDataFrame(chunk, geometry=geometry, crs=CRS)
    
    return gdf



def chunk_tsv_by_day(tsv_path, target_rows=100000):
    """
    Generator that yields complete days from TSV file.
    Reads ~target_rows, then continues to end of that day.
    Much faster than grouping every chunk.
    
    Args:
        tsv_path: Path to TSV file
        target_rows: Approximate rows per chunk (will read more to complete day)
    
    Yields:
        GeoDataFrame with complete day(s) of data
    """
    try:
        reader = pd.read_csv(tsv_path, sep='\t', chunksize=target_rows, iterator=True)
        
        buffer = []
        last_day = None
        
        for chunk in reader:
            # Parse timestamps
            chunk['datetime'] = pd.to_datetime(chunk[TIME_COL], unit='s', errors='coerce')
            chunk['date'] = chunk['datetime'].dt.date
            
            # Add this chunk to buffer
            buffer.append(chunk)
            
            # Check the last date in this chunk
            valid_dates = chunk['date'].dropna()
            if len(valid_dates) == 0:
                continue
            
            chunk_last_day = valid_dates.iloc[-1]
            
            # If this is our first chunk, just note the last day and continue
            if last_day is None:
                last_day = chunk_last_day
                continue
            
            # If the last day changed, we need to split the buffer
            if chunk_last_day != last_day:
                # Concatenate all buffered chunks
                combined = pd.concat(buffer, ignore_index=True)
                
                # Split: complete days vs. incomplete final day
                complete_days = combined[combined['date'] < chunk_last_day]
                incomplete_day = combined[combined['date'] == chunk_last_day]
                
                # Yield complete days if any
                if len(complete_days) > 0:
                    geometry = [Point(xy) for xy in zip(complete_days[LON_COL], complete_days[LAT_COL])]
                    gdf = gpd.GeoDataFrame(complete_days, geometry=geometry, crs=CRS)
                    yield gdf
                
                # Keep incomplete day for next iteration
                buffer = [incomplete_day] if len(incomplete_day) > 0 else []
                last_day = chunk_last_day
        
        # Yield final buffer
        if buffer:
            complete_day = pd.concat(buffer, ignore_index=True)
            if len(complete_day) > 0:
                geometry = shapely.points(
                    complete_day[LON_COL].to_numpy(),  # complete_day not complete_days
                    complete_day[LAT_COL].to_numpy()
                )
                gdf = gpd.GeoDataFrame(complete_day, geometry=geometry, crs=CRS)
                yield gdf
                
    except Exception as e:
        print(f"    Error in chunk_tsv_by_day: {e}")
        raise


def aggregate_daily_visits(points_gdf, layer_name):
    """Aggregate unique device visits by day."""
    if len(points_gdf) == 0:
        return pd.DataFrame()
    
    # Get unique devices per date/feature_id with their demographics - grabs the first instance of each
    unique_devices = points_gdf.groupby(['date', 'feature_id', DEVICE_COL]).first().reset_index()

    # Now aggregate across unique devices
    daily = unique_devices.groupby(['date', 'feature_id']).agg(
    unique_visits=(DEVICE_COL, 'nunique'),
    sum_median_income=('median_income', lambda x: x.sum(skipna=True)),
    sum_average_income=('average_income', lambda x: x.sum(skipna=True)),
    sum_pct_education=('pct_education', lambda x: x.sum(skipna=True)),
    n_demographic_matched=('median_income', 'count')  # count non-NaN
).reset_index()
    
    daily['layer'] = layer_name

    return daily


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_tsv_file(tsv_path, output_dir, by_polygon=False, gpkg_path=None):
    """Process a single TSV file with city-level filtering."""
    
    state_abbrev, state_name, gpkg_path, lookup_path = get_paths_from_filename(tsv_path)

    print(f"  State: {state_name} ({state_abbrev})")
    print(f"  GPKG:  {gpkg_path}")
    print(f"  Lookup: {lookup_path}")

    DEMOGRAPHICS_DF = pd.read_parquet(lookup_path).set_index('Hashed Ubermedia Id')
    DEMO_COLS = ['median_income', 'average_income', 'pct_education']
    for col in DEMO_COLS:
        if col in DEMOGRAPHICS_DF.columns:
            DEMOGRAPHICS_DF[col] = DEMOGRAPHICS_DF[col].where(DEMOGRAPHICS_DF[col] >= 0, other=pd.NA)

    print(f"\n{'='*80}")
    print(f"Processing: {tsv_path.name}")
    print(f"{'='*80}")
    
    try:
        city_name, state = parse_tsv_filename(tsv_path)
        print(f"  City Name: {city_name}")
        print(f"  State: {state}")
    except Exception as e:
        print(f"  ❌ Error parsing filename: {e}")
        return {}
    
    tsv_parts = find_tsv_parts(tsv_path)
    print(f"  Found {len(tsv_parts)} TSV part(s)")
    
    # --- Geojson loading with fallback scanner ---
    geojson_path = find_geojson(city_name, state)
    
    def _build_polygons_from_geojson(gj_path):
        """Helper: load cities from a geojson, build city_layer_map and all_polygons."""
        cities = load_geojson_cities(gj_path)
        print(f"    DEBUG {gj_path.name}: cities={cities}") 
        if not cities:
            return None, None, None, None
        clm = load_all_city_polygons(state, cities, NEW_LAYERS, gpkg_path)
        if not clm:
            return None, None, None, None
        poly_attrs = pd.concat([
            polygons[['feature_id'] + [col for col in polygons.columns
                                        if col not in ['geometry', 'feature_id']]]
            .assign(source_city=city, layer=layer)
            .drop_duplicates('feature_id')
            for (city, layer), polygons in clm.items()
        ], ignore_index=True)
        all_poly = gpd.GeoDataFrame(
            pd.concat([
                polygons.assign(layer=layer, source_city=city)
                for (city, layer), polygons in clm.items()
            ], ignore_index=True),
            crs=CRS
        )
        all_poly['source_city'] = all_poly['source_city'].astype('category')
        all_poly.sindex
        return cities, clm, poly_attrs, all_poly

    def _sample_points_from_tsv(tsv_part, n=500):
        """Read the first chunk and return a small GeoDataFrame sample."""
        try:
            sample_df = pd.read_csv(tsv_part, sep='\t', nrows=n)
            sample_df['datetime'] = pd.to_datetime(sample_df[TIME_COL], unit='s', errors='coerce')
            sample_df['date'] = sample_df['datetime'].dt.date
            geom = shapely.points(sample_df[LON_COL].to_numpy(), sample_df[LAT_COL].to_numpy())
            return gpd.GeoDataFrame(sample_df, geometry=geom, crs=CRS)
        except Exception as e:
            print(f"    ⚠️  Could not sample TSV for probe: {e}")
            return None

    def _any_points_match(sample_gdf, all_poly):
        """Return True if at least one point intersects any polygon."""
        if sample_gdf is None or all_poly is None or len(all_poly) == 0:
            return False
        try:
            result = gpd.sjoin(sample_gdf, all_poly, how='inner', predicate='intersects')
            return len(result) > 0
        except Exception:
            return False

    # Try the matched geojson first, then fall back to scanning all geojsons in state
    cities = city_layer_map = polygon_attrs_lookup = all_polygons = None
    used_geojson = None

    if geojson_path:
        print(f"  ✅ Found geojson: {geojson_path.name}")
        cities, city_layer_map, polygon_attrs_lookup, all_polygons = \
            _build_polygons_from_geojson(geojson_path)

        print(f"  DEBUG cities: {cities}")
        print(f"  DEBUG all_polygons rows: {0 if all_polygons is None else len(all_polygons)}")
                
        if all_polygons is not None and len(tsv_parts) > 0:
            sample = _sample_points_from_tsv(tsv_parts[0])
            if _any_points_match(sample, all_polygons):
                used_geojson = geojson_path
                print(f"  ✅ Probe confirmed spatial match: {geojson_path.name}")
            else:
                print(f"  ⚠️  Probe found no spatial matches for {geojson_path.name}, will scan state...")
                cities = city_layer_map = polygon_attrs_lookup = all_polygons = None
        elif all_polygons is not None:
            used_geojson = geojson_path  # no TSV parts to probe, trust the match
    else:
        print(f"  ⚠️  No geojson found via filename match, scanning state geojsons...")

    # Fallback: scan all geojsons in the state directory
    if used_geojson is None:
        state_dir_name = state.replace(' ', '_')
        candidate_geojsons = []
        for state_dir in STATES_DIR.glob(f"{state_dir_name}_*"):
            candidate_geojsons.extend(sorted(state_dir.glob("*.geojson")))
        
        # Skip the one we already tried
        if geojson_path:
            candidate_geojsons = [p for p in candidate_geojsons if p != geojson_path]
        
        print(f"  🔍 Scanning {len(candidate_geojsons)} candidate geojson(s) in {state}...")
        sample = _sample_points_from_tsv(tsv_parts[0]) if tsv_parts else None

        for candidate in candidate_geojsons:
            c_cities, c_clm, c_attrs, c_polys = _build_polygons_from_geojson(candidate)
            if c_polys is None:
                continue
            if _any_points_match(sample, c_polys):
                print(f"  ✅ Fallback match found: {candidate.name}")
                cities, city_layer_map = c_cities, c_clm
                polygon_attrs_lookup, all_polygons = c_attrs, c_polys
                used_geojson = candidate
                break
            else:
                print(f"    ↳ No match: {candidate.name}")

    if used_geojson is None or city_layer_map is None:
        print(f"  ❌ No spatial match found in any geojson for {state}, skipping.")
        return {}

    print(f"  Using geojson: {used_geojson.name}")
    city_code = used_geojson.stem
    variant_match = re.search(r'_([ab])$', city_name)
    if variant_match and variant_match.group(1) == 'b' and not city_code.endswith('_b'):
        city_code = f"{city_code}_b"

    layers_to_process = [
        layer for layer in NEW_LAYERS
        if not (output_dir / f"{city_code}_{layer}_{state}_daily.parquet").exists()
    ]

    if not layers_to_process:
        print(f" ⏭️ All output files already exist for {city_code}, skipping.")
        return {'skipped': True, 'city_code': city_code}
    
    print(f"  ✅ Loaded polygons for {len(city_layer_map)} city/layer combinations")

    # === Rest of processing (unchanged) ===
    daily_aggregations = {key: [] for key in city_layer_map.keys()}
    stats = {
        'total_chunks': 0,
        'total_points': 0,
        'points_per_layer': defaultdict(int),
        'state': state,
        'city_name': city_name
    }

    for part_num, tsv_part in enumerate(tsv_parts, 1):
        if not tsv_part.exists():
            print(f"  ⚠️  Part {part_num} does not exist: {tsv_part.name}")
            continue
        
        print(f"  📄 Processing part {part_num}/{len(tsv_parts)}: {tsv_part.name}")
        
        try:
            chunk_count = 0
            for gdf_chunk in chunk_tsv_by_day(tsv_part, target_rows=CHUNK_SIZE):
                chunk_count += 1
                stats['total_chunks'] += 1
                stats['total_points'] += len(gdf_chunk)

                gdf_chunk = gdf_chunk.join(DEMOGRAPHICS_DF, on=DEVICE_COL, how='left')

                all_matches = gpd.sjoin(
                    gdf_chunk, all_polygons, how='inner', predicate='intersects'
                )
                
                all_matches['layer'] = all_matches['layer'].astype(str)

                for layer in all_matches['layer'].unique():
                    stats['points_per_layer'][layer] += (all_matches['layer'] == layer).sum()

                for layer in layers_to_process:
                    layer_matches = all_matches[all_matches['layer'] == layer]
                    if len(layer_matches) == 0:
                        continue
                    
                    for source_city, city_layer_matches in layer_matches.groupby('source_city', observed=True):
                        daily_chunk = aggregate_daily_visits(city_layer_matches, layer)
                        if len(daily_chunk) > 0:
                            daily_aggregations[(source_city, layer)].append(daily_chunk)
                
                del gdf_chunk
                
                if chunk_count % 10 == 0:
                    print(f"    Processed {chunk_count} day-chunks...")
        
        except Exception as e:
            print(f"  ⚠️  Error reading {tsv_part.name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n  ✅ Processed {stats['total_chunks']} chunks, {stats['total_points']:,} total points")
    
    if stats['points_per_layer']:
        print(f"\n  📍 Points matched per layer:")
        for layer, count in sorted(stats['points_per_layer'].items()):
            print(f"    {layer}: {count:,} points")
    else:
        print(f"\n  ⚠️  No points matched any layer!")
        return stats
    
    output_dir.mkdir(parents=True, exist_ok=True)
    saved_files = []
    layer_aggregations = defaultdict(list)

    for (city, layer), chunk_list in daily_aggregations.items():
        if not chunk_list:
            continue
        
        combined = pd.concat(chunk_list, ignore_index=True)
        
        final = combined.groupby(['date', 'layer', 'feature_id']).agg(
            unique_visits=('unique_visits', 'sum'),
            sum_median_income=('sum_median_income', 'sum'),
            sum_average_income=('sum_average_income', 'sum'),
            sum_pct_education=('sum_pct_education', 'sum'),
            n_demographic_matched=('n_demographic_matched', 'sum')
        ).reset_index()

        final = final.merge(
            polygon_attrs_lookup, on=['feature_id', 'layer'], how='left', suffixes=('', '_poly')
        )
        final['date'] = pd.to_datetime(final['date'])
        final = final.sort_values('date')
        layer_aggregations[layer].append(final)

    for layer, city_dfs in layer_aggregations.items():
        if not city_dfs:
            continue
        layer_combined = pd.concat(city_dfs, ignore_index=True)
        output_file = output_dir / f"{city_code}_{layer}_{state}_daily.parquet"
        layer_combined.to_parquet(output_file, index=False, engine='pyarrow')
        saved_files.append(output_file)
        print(f"  ✅ Saved {layer} (all cities): {output_file.name} ({len(layer_combined)} rows)")

    stats['saved_files'] = saved_files
    return stats

def main():
    parser = argparse.ArgumentParser(description='Process cellphone mobility data with city-level filtering')
    
    parser.add_argument('tsv_file', type=Path, help='TSV file to process')
    parser.add_argument('gpkg_file', type=Path, help='GeoPackage file')
    parser.add_argument('output_dir', type=Path, help='Output directory')
    parser.add_argument('--by-polygon', action='store_true', help='Aggregate by polygon')
    
    args = parser.parse_args()
    
    if not args.tsv_file.exists():
        print(f"❌ TSV file not found: {args.tsv_file}")
        sys.exit(1)
    
    if not args.gpkg_file.exists():
        print(f"❌ GeoPackage file not found: {args.gpkg_file}")
        sys.exit(1)
    
    # Process individual city file
    stats = process_tsv_file(
        args.tsv_file,
        args.output_dir,
        by_polygon=args.by_polygon
    )
    
    if stats.get('saved_files'):
        print(f"\n✅ Processing complete!")
        print(f"   Total points: {stats['total_points']:,}")
        print(f"   Files saved: {len(stats['saved_files'])}")
        sys.exit(0)
    else:
        print(f"\n⚠️  Processing complete but no output files created")
        sys.exit(1)


def combine_all_cities(output_dir, by_polygon=False):
    """Combine all city CSV files into layer-level files."""
    
    suffix = '_by_polygon' if by_polygon else ''
    pattern = f"*{suffix}_daily.parquet"

    all_files = list(output_dir.glob(pattern))
    
    if not all_files:
        print(f"No files found matching {pattern}")
        return
    
    print(f"\n{'='*80}")
    print(f"Combining {len(all_files)} city files")
    print(f"{'='*80}")
    
    # Group files by layer
    layer_files = defaultdict(list)
    for f in all_files:
        # Filename format: {city_code}_{layer}_{state}_daily.parquet
        # parts[-1] = 'daily', parts[-2] = state name(s), parts[-3] = layer
        parts = f.stem.split('_')
        layer = parts[-3] if len(parts) >= 4 else 'unknown'
        layer_files[layer].append(f)
    
    # Combine each layer
    combined_files = []
    for layer, files in layer_files.items():
        print(f"\n  Combining {len(files)} files for layer: {layer}")
        
        dfs = []
        for f in files:
            df = pd.read_parquet(f)
            dfs.append(df)
        
        combined = pd.concat(dfs, ignore_index=True)
        
        # Sort by date and city
        combined['date'] = pd.to_datetime(combined['date'])
        combined = combined.sort_values(['date', 'city_name'])
        
        output_file = output_dir / f"ALL_CITIES_{layer}{suffix}_daily.parquet"
        combined.to_parquet(output_file, index=False, engine='pyarrow')
        
        print(f"  ✅ Saved: {output_file.name} ({len(combined)} rows)")
    
    return combined_files


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == '--combine':
        # Usage: python script.py --combine output_dir [--by-polygon]
        output_dir = Path(sys.argv[2])
        by_polygon = '--by-polygon' in sys.argv
        combine_all_cities(output_dir, by_polygon)
    else:
        main()

