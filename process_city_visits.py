#!/usr/bin/env python3
"""
Process city visit data from TSV files, join with geojson and geopackage layers,
and generate monthly visit HTML visualizations.
"""

import pandas as pd
import geopandas as gpd
import plotly.express as px
from pathlib import Path
import json
import re
from shapely.validation import make_valid
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

TSV_DIR = Path("test_cities")
STATES_DIR = Path("states")
GPKG_DIR = Path("unified_geopackages_updated")
OUT_DIR = Path("html_outputs")
CHUNK_SIZE = 500_000
TIME_COL = "Unix Timestamp of Visit"
LON_COL = "Lon of Visit"
LAT_COL = "Lat of Visit"

# US State FIPS codes
STATE_FIPS = {
    'Alabama': '01', 'Alaska': '02', 'Arizona': '04', 'Arkansas': '05',
    'California': '06', 'Colorado': '08', 'Connecticut': '09', 'Delaware': '10',
    'Florida': '12', 'Georgia': '13', 'Hawaii': '15', 'Idaho': '16',
    'Illinois': '17', 'Indiana': '18', 'Iowa': '19', 'Kansas': '20',
    'Kentucky': '21', 'Louisiana': '22', 'Maine': '23', 'Maryland': '24',
    'Massachusetts': '25', 'Michigan': '26', 'Minnesota': '27', 'Mississippi': '28',
    'Missouri': '29', 'Montana': '30', 'Nebraska': '31', 'Nevada': '32',
    'New Hampshire': '33', 'New Jersey': '34', 'New Mexico': '35', 'New York': '36',
    'North Carolina': '37', 'North Dakota': '38', 'Ohio': '39', 'Oklahoma': '40',
    'Oregon': '41', 'Pennsylvania': '42', 'Rhode Island': '44', 'South Carolina': '45',
    'South Dakota': '46', 'Tennessee': '47', 'Texas': '48', 'Utah': '49',
    'Vermont': '50', 'Virginia': '51', 'Washington': '53', 'West Virginia': '54',
    'Wisconsin': '55', 'Wyoming': '56', 'District of Columbia': '11'
}

# State abbreviation to full name mapping
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
    'SD': 'South Dakota', 'TN': 'Tennessee', 'TX': 'Texas', 'UT': 'Utah',
    'VT': 'Vermont', 'VA': 'Virginia', 'WA': 'Washington', 'WV': 'West Virginia',
    'WI': 'Wisconsin', 'WY': 'Wyoming', 'DC': 'District of Columbia'
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_tsv_filename(filename):
    """
    Parse TSV filename to extract city code, city name, and state.
    
    Examples: 
        10057518_Tulsa_3_Oklahoma_pin_report_000.tsv -> ('Tulsa_3', 'Oklahoma')
        10056845_Casper_Wyoming_pin_report_000.tsv -> ('Casper', 'Wyoming')
        10057155_Roswell_a_New_Mexico_pin_report.tsv -> ('Roswell', 'New Mexico')
    
    Logic:
    - Skip initial numbers (first part)
    - Collect city name parts (including numbered suffixes like _3) until we hit _a, _b (exclude these), or identify the state name
    - State name is the part(s) between the city and _pin_report that match a known state
    """
    stem = filename.stem  # Remove .tsv
    
    # Remove _pin_report_XXX suffix if present
    stem = re.sub(r'_pin_report(_\d+)?$', '', stem)
    
    # Split by underscores
    parts = stem.split('_')
    
    # First part should be the numeric code (skip it)
    if not parts or not parts[0].isdigit():
        raise ValueError(f"Cannot parse numeric code from: {filename}")
    
    # Start from index 1 (skip the initial number)
    idx = 1
    city_parts = []
    state_parts = []
    
    max_iterations = len(parts) * 2  # Safety limit
    iteration = 0
    
    # Collect parts for city_code
    # Include numbered suffixes but stop when we hit _a, _b (and exclude them), or when we find the state name
    while idx < len(parts) and iteration < max_iterations:
        iteration += 1
        part = parts[idx]
        
        # Check for _a or _b suffixes FIRST (exclude them from city_code)
        if part in ('a', 'b'):
            idx += 1
            # The next part(s) should be the state
            if idx < len(parts):
                # Check for two-word state
                if idx + 1 < len(parts):
                    potential_state_double = f"{parts[idx]} {parts[idx + 1]}"
                    if potential_state_double in STATE_FIPS:
                        state_parts = [parts[idx], parts[idx + 1]]
                        break
                # Check for single-word state
                if parts[idx] in STATE_FIPS:
                    state_parts = [parts[idx]]
                    break
                # Check for state abbreviation
                if parts[idx].upper() in STATE_ABBREV:
                    state_parts = [STATE_ABBREV[parts[idx].upper()]]
                    break
            break
        
        # Only start checking for states after we have at least one city part
        # This handles cases like "Idaho_Falls_ID" where "Idaho" is both city and state
        elif len(city_parts) > 0:
            # Check if this could be the start of a state name
            # Try matching single-word and two-word states
            potential_state_single = part
            potential_state_double = None
            if idx + 1 < len(parts):
                potential_state_double = f"{part} {parts[idx + 1]}"
            
            # IMPORTANT: Check abbreviations FIRST (most specific)
            # This prevents "Falls_ID" from treating "Falls" as part of state name
            if idx + 1 < len(parts) and parts[idx + 1].upper() in STATE_ABBREV:
                # Next part is a state abbreviation, so current part is last of city
                city_parts.append(part)
                idx += 1
                state_parts = [STATE_ABBREV[parts[idx].upper()]]
                break
            # Check if current part itself is state abbreviation
            elif potential_state_single.upper() in STATE_ABBREV:
                # Expand abbreviation to full name
                state_parts = [STATE_ABBREV[potential_state_single.upper()]]
                break
            # Check two-word states (they're more specific)
            elif potential_state_double and potential_state_double in STATE_FIPS:
                state_parts = [part, parts[idx + 1]]
                break
            # Check single-word states
            elif potential_state_single in STATE_FIPS:
                state_parts = [part]
                break
            else:
                # Not a state, add to city parts
                city_parts.append(part)
                idx += 1
        else:
            # First iteration - just add to city parts
            city_parts.append(part)
            idx += 1
    
    if iteration >= max_iterations:
        raise ValueError(f"Parser exceeded maximum iterations for: {filename}. Parts: {parts}, idx: {idx}, city_parts: {city_parts}")
    
    if not city_parts:
        raise ValueError(f"Cannot parse city name from: {filename}")
    
    if not state_parts:
        raise ValueError(f"Cannot parse state name from: {filename}")
    
    city_code = '_'.join(city_parts)
    state = ' '.join(state_parts)
    
    return city_code, state


def find_tsv_parts(tsv_path):
    """
    Find all parts of a TSV file (e.g., _000, _001, _002).
    
    Returns list of Path objects sorted by part number.
    """
    # Check if this is already a _XXX file
    if re.search(r'_pin_report_\d{3}\.tsv$', tsv_path.name):
        # Find all related parts by replacing the number with a wildcard
        base = re.sub(r'_pin_report_\d{3}\.tsv$', '_pin_report_', tsv_path.name)
        parts = sorted(tsv_path.parent.glob(f"{base}*.tsv"))
        return parts
    elif '_pin_report.tsv' in tsv_path.name:
        # This is a single file, but check if numbered parts exist
        base = tsv_path.name.replace('.tsv', '_')
        parts = sorted(tsv_path.parent.glob(f"{base}*.tsv"))
        if parts:
            return parts
        else:
            return [tsv_path]
    else:
        # No _pin_report suffix, just return the file
        return [tsv_path]


def find_geojson(city_code, state):
    """
    Find the geojson file for a given city code and state.
    
    Format: states/<state_name>_<numbers>/<city_code>.geojson
    Note: State names in directories use underscores instead of spaces
    Also tries replacing underscores with spaces in city names (Idaho_Falls -> Idaho Falls)
    """
    # Replace spaces with underscores for directory name
    state_dir_name = state.replace(' ', '_')
    
    # Look for state directories
    state_dirs = list(STATES_DIR.glob(f"{state_dir_name}_*"))
    
    if not state_dirs:
        print(f"  ❌ No state directory found for {state}")
        return None
    
    # Try different variations of the city code
    city_variations = [
        city_code,  # Original: Idaho_Falls
        city_code.replace('_', ' '),  # With spaces: Idaho Falls
        city_code.replace('_', '')  # No separator: IdahoFalls
    ]
    
    # Search in all matching state directories
    for state_dir in state_dirs:
        for city_variant in city_variations:
            geojson_path = state_dir / f"{city_variant}.geojson"
            if geojson_path.exists():
                return geojson_path
    
    print(f"  ❌ No geojson found for {city_code} in {state}")
    print(f"      Tried: {', '.join(city_variations)}")
    return None


def load_geojson_cities(geojson_path):
    """
    Load cities from geojson file.
    """
    with open(geojson_path, 'r') as f:
        data = json.load(f)
    
    if not data.get('features'):
        return []
    
    # Extract cities from all_cities property
    all_cities = data['features'][0]['properties'].get('all_cities', '')
    cities = [c.strip() for c in all_cities.split(',') if c.strip()]
    
    return cities


def normalize_city_name(city_name):
    """
    Normalize city name for matching:
    - Remove numbered suffixes like _2, _3
    - Replace underscores with spaces (Idaho_Falls -> idaho falls)
    - Convert to lowercase for case-insensitive matching
    - Strip whitespace
    """
    # Remove _<number> suffix
    normalized = re.sub(r'_\d+$', '', city_name)
    # Replace underscores with spaces
    normalized = normalized.replace('_', ' ')
    return normalized.lower().strip()


def load_geopackage_layers(state, layer_name, city=None):
    """
    Load a specific layer from the geopackage for a state.
    Optionally filter to a specific city (case-insensitive, ignores numbered suffixes).
    """
    fips = STATE_FIPS.get(state)
    if not fips:
        print(f"  ❌ No FIPS code found for state: {state}")
        return None
    
    gpkg_path = GPKG_DIR / f"state_{fips}_unified.gpkg"
    if not gpkg_path.exists():
        print(f"  ❌ Geopackage not found: {gpkg_path}")
        return None
    
    try:
        gdf = gpd.read_file(gpkg_path, layer=layer_name)
        gdf = gdf.set_crs(4326, allow_override=True)
        
        # Validate polygons for non-unified_park_area layers
        if layer_name != 'unified_park_area':
            gdf['geometry'] = gdf['geometry'].apply(make_valid)
        
        # Filter to city if specified
        if city:
            normalized_city = normalize_city_name(city)
            
            if 'city_name' in gdf.columns:
                # Create normalized version for comparison
                gdf['_normalized_city'] = gdf['city_name'].apply(normalize_city_name)
                gdf = gdf[gdf['_normalized_city'] == normalized_city]
                gdf = gdf.drop(columns=['_normalized_city'])
            elif 'city' in gdf.columns:
                gdf['_normalized_city'] = gdf['city'].apply(normalize_city_name)
                gdf = gdf[gdf['_normalized_city'] == normalized_city]
                gdf = gdf.drop(columns=['_normalized_city'])
        
        return gdf
    except Exception as e:
        print(f"  ⚠️  Error loading layer {layer_name}: {e}")
        return None


def process_tsv_chunk(chunk_df):
    """
    Convert a TSV chunk to a GeoDataFrame with points.
    """
    # Ensure required columns exist
    if TIME_COL not in chunk_df.columns or LON_COL not in chunk_df.columns or LAT_COL not in chunk_df.columns:
        raise ValueError(f"Missing required columns in TSV")
    
    # Convert timestamp to datetime
    chunk_df[TIME_COL] = pd.to_datetime(chunk_df[TIME_COL], unit='s')
    
    # Create GeoDataFrame
    gdf_points = gpd.GeoDataFrame(
        chunk_df,
        geometry=gpd.points_from_xy(chunk_df[LON_COL], chunk_df[LAT_COL]),
        crs="EPSG:4326"
    )
    
    return gdf_points


def aggregate_monthly(points_gdf, layer_name):
    """
    Aggregate points by month.
    """
    if len(points_gdf) == 0:
        return pd.DataFrame(columns=['month', 'visits', 'layer'])
    
    monthly = points_gdf.groupby(
        pd.Grouper(key=TIME_COL, freq='ME')  # Month end frequency
    ).size().reset_index()
    monthly.columns = ['month', 'visits']
    monthly['layer'] = layer_name
    
    return monthly


def aggregate_chunk_monthly(points_gdf, layer_name):
    """
    Aggregate a chunk of points by month immediately.
    Returns monthly counts rather than storing all points.
    """
    if len(points_gdf) == 0:
        return pd.DataFrame(columns=['month', 'visits', 'layer'])
    
    monthly = points_gdf.groupby(
        pd.Grouper(key=TIME_COL, freq='ME')  # Month end frequency
    ).size().reset_index()
    monthly.columns = ['month', 'visits']
    monthly['layer'] = layer_name
    
    return monthly


def load_all_city_polygons(state, cities, layers):
    """
    Load all polygons for all cities and layers.
    Returns a dictionary mapping (city, layer) -> polygons.
    Ensures each polygon knows which city it belongs to.
    Note: Uses normalized city names (lowercase, no numbered suffixes) for matching.
    """
    city_layer_map = {}
    
    for city in cities:
        for layer in layers:
            polygons = load_geopackage_layers(state, layer, city)
            if polygons is not None and len(polygons) > 0:
                # Ensure we have a 'city' column for verification
                # Use the ORIGINAL city name (with suffix) for tracking
                if 'city_name' in polygons.columns and 'city' not in polygons.columns:
                    polygons = polygons.copy()
                    polygons['city'] = city  # Use original city name with suffix
                elif 'city' not in polygons.columns:
                    polygons = polygons.copy()
                    polygons['city'] = city
                key = (city, layer)
                city_layer_map[key] = polygons
    
    return city_layer_map


# ============================================================================
# MAIN PROCESSING
# ============================================================================

def main():
    """
    Main processing function.
    """
    # Create output directory
    OUT_DIR.mkdir(exist_ok=True)
    
    # Get all TSV files
    tsv_files = sorted(TSV_DIR.glob("*.tsv"))
    
    # Group TSV files by their base name (to handle multi-part files)
    processed_bases = set()
    
    for tsv_file in tsv_files:
        # Skip if we've already processed this base
        # Remove both _pin_report_XXX and _pin_report patterns
        base_name = re.sub(r'_pin_report(_\d{3})?\.tsv$', '', tsv_file.name)
        if base_name in processed_bases:
            continue
        
        print(f"\n{'='*80}")
        print(f"Processing: {tsv_file.name}")
        print(f"{'='*80}")
        
        try:
            # Parse filename
            city_code, state = parse_tsv_filename(tsv_file)
            print(f"  City Code: {city_code}")
            print(f"  State: {state}")
            
            # Find all parts of this TSV
            tsv_parts = find_tsv_parts(tsv_file)
            print(f"  Found {len(tsv_parts)} TSV part(s)")
            
            # Mark as processed
            processed_bases.add(base_name)
            
            # Find corresponding geojson
            geojson_path = find_geojson(city_code, state)
            if not geojson_path:
                continue
            
            print(f"  ✅ Found geojson: {geojson_path}")
            
            # Load cities from geojson
            cities = load_geojson_cities(geojson_path)
            print(f"  Found {len(cities)} cities: {', '.join(cities)}")
            
            # Get geopackage layers
            fips = STATE_FIPS.get(state)
            if not fips:
                print(f"  ❌ No FIPS code for {state}")
                continue
            
            gpkg_path = GPKG_DIR / f"state_{fips}_unified.gpkg"
            if not gpkg_path.exists():
                print(f"  ❌ Geopackage not found: {gpkg_path}")
                continue
            
            # Get available layers
            layers = ['unified_park_area', 'parkserve_geom', 'osm_geom', 'indoor_gyms']
            
            # ================================================================
            # Load all city/layer polygons once
            # ================================================================
            print(f"\n  📐 Loading all city/layer polygons...")
            city_layer_map = load_all_city_polygons(state, cities, layers)
            
            if not city_layer_map:
                print(f"  ❌ No polygons found for any city/layer combination")
                continue
            
            print(f"  ✅ Loaded polygons for {len(city_layer_map)} city/layer combinations")
            
            # ================================================================
            # Initialize storage for monthly aggregations (not raw points!)
            # ================================================================
            city_layer_monthly = {key: [] for key in city_layer_map.keys()}
            
            # ================================================================
            # Process TSV files chunk by chunk, aggregating immediately
            # ================================================================
            print(f"\n  📊 Processing TSV files in chunks...")
            total_chunks = 0
            chunk_point_counts = {}  # Track points found per city/layer
            
            for tsv_part in tsv_parts:
                print(f"    Reading {tsv_part.name}...")
                
                # Validate file before processing
                if not tsv_part.exists():
                    print(f"    ⚠️  File does not exist: {tsv_part}")
                    continue
                
                # Try to get actual file size (handles sparse files better)
                try:
                    file_size = tsv_part.stat().st_size
                    print(f"      File stat size: {file_size / (1024*1024):.2f} MB")
                except Exception as e:
                    print(f"    ⚠️  Could not stat file: {e}")
                    file_size = 0
                
                # Try to peek at first few lines to validate format
                # This is more reliable than size check for sparse files
                try:
                    with open(tsv_part, 'r') as f:
                        first_line = f.readline().strip()
                        if not first_line:
                            print(f"    ⚠️  File appears to be empty (no content)")
                            print(f"    ⚠️  This might be a sparse file or filesystem issue")
                            print(f"    ℹ️  Skipping this file - you may want to check it manually:")
                            print(f"        ls -lh {tsv_part}")
                            print(f"        file {tsv_part}")
                            print(f"        head -5 {tsv_part}")
                            continue
                        
                        print(f"      First line preview: {first_line[:100]}...")
                        
                        # Check if it looks like a TSV (has tabs)
                        if '\t' not in first_line:
                            print(f"    ⚠️  File doesn't appear to be tab-separated")
                            continue
                        
                        # Check for required columns
                        columns = first_line.split('\t')
                        print(f"      Found {len(columns)} columns")
                        if TIME_COL not in columns or LON_COL not in columns or LAT_COL not in columns:
                            print(f"    ⚠️  Missing required columns. Expected: {TIME_COL}, {LON_COL}, {LAT_COL}")
                            print(f"    ⚠️  Found columns: {columns[:10]}")
                            continue
                        
                        # Check if there's at least one data row
                        second_line = f.readline().strip()
                        if not second_line:
                            print(f"    ⚠️  File has header but no data rows")
                            continue
                        
                        print(f"      ✅ File validation passed")
                        
                except FileNotFoundError:
                    print(f"    ⚠️  File not found when trying to open")
                    continue
                except PermissionError:
                    print(f"    ⚠️  Permission denied when trying to open file")
                    continue
                except Exception as e:
                    print(f"    ⚠️  Error validating file: {e}")
                    continue
                
                try:
                    for chunk_num, chunk in enumerate(pd.read_csv(tsv_part, sep='\t', chunksize=CHUNK_SIZE)):
                        total_chunks += 1
                        gdf_chunk = process_tsv_chunk(chunk)
                        
                        # For each city/layer combination, filter and aggregate this chunk
                        for (city, layer), polygons in city_layer_map.items():
                            # Spatial join for this specific city/layer
                            # The polygons already have 'city' column to ensure correct assignment
                            points_in_layer = gpd.sjoin(
                                gdf_chunk,
                                polygons,
                                how='inner',
                                predicate='within'
                            )
                            
                            if len(points_in_layer) > 0:
                                # Debug: Track how many points we're finding
                                key = (city, layer)
                                if key not in chunk_point_counts:
                                    chunk_point_counts[key] = 0
                                
                                before_filter = len(points_in_layer)
                                
                                # Verify we're only getting points for THIS city
                                # Use normalized comparison since polygons might have different case/format
                                if 'city' in points_in_layer.columns:
                                    normalized_target = normalize_city_name(city)
                                    points_in_layer = points_in_layer[
                                        points_in_layer['city'].apply(normalize_city_name) == normalized_target
                                    ]
                                
                                after_filter = len(points_in_layer)
                                chunk_point_counts[key] += after_filter
                                
                                if before_filter != after_filter:
                                    print(f"      DEBUG: Filtered {before_filter - after_filter} points not in {city} ({layer})")
                                
                                if len(points_in_layer) > 0:
                                    # Aggregate immediately - don't store raw points!
                                    monthly = aggregate_chunk_monthly(points_in_layer, layer)
                                    city_layer_monthly[(city, layer)].append(monthly)
                        
                        # Free memory
                        del gdf_chunk
                        
                        if (chunk_num + 1) % 10 == 0:
                            print(f"      Processed {chunk_num + 1} chunks...")
                        
                except Exception as e:
                    print(f"    ⚠️  Error reading {tsv_part.name}: {e}")
                    continue
            
            print(f"  ✅ Processed {total_chunks} total chunks")
            
            # Debug: Print summary of points found
            if chunk_point_counts:
                print(f"\n  📍 Points found per city/layer:")
                for (city, layer), count in sorted(chunk_point_counts.items()):
                    print(f"    {city} / {layer}: {count:,} points")
            else:
                print(f"\n  ⚠️  No points found in any city/layer combination!")
            
            # ================================================================
            # Sum monthly aggregations for each city
            # ================================================================
            cities_with_data = []
            
            for city in cities:
                print(f"\n  📍 Processing city: {city}")
                city_visits = []
                
                # Process each layer for this city
                for layer in layers:
                    key = (city, layer)
                    
                    if key not in city_layer_map:
                        print(f"    Layer {layer}: No polygons")
                        continue
                    
                    if not city_layer_monthly[key]:
                        continue
                    
                    # Combine all monthly aggregations and sum by month
                    all_monthly = pd.concat(city_layer_monthly[key], ignore_index=True)
                    final_monthly = all_monthly.groupby(['month', 'layer'])['visits'].sum().reset_index()
                    
                    total_visits = final_monthly['visits'].sum()
                    print(f"      {layer}: {total_visits:,} visits")
                    
                    city_visits.append(final_monthly)
                
                # Generate HTML visualization
                if city_visits:
                    city_visits_df = pd.concat(city_visits, ignore_index=True)
                    
                    fig = px.line(
                        city_visits_df,
                        x='month',
                        y='visits',
                        color='layer',
                        line_dash='layer',
                        title=f"{city} - Monthly Visits"
                    )
                    
                    html_path = OUT_DIR / f"{city_code}_{city.replace(' ', '_')}.html"
                    fig.write_html(html_path)
                    print(f"      ✅ Saved: {html_path.name}")
                    cities_with_data.append(city)
            
            # ================================================================
            # FALLBACK: If no cities had data, search other geojsons in the state
            # ================================================================
            if not cities_with_data and not chunk_point_counts:
                print(f"\n  ⚠️  No matches found with {geojson_path.name}")
                print(f"  🔍 Searching for matching geojson in other {state} cities...")
                
                # Find all geojsons in this state
                state_dir_name = state.replace(' ', '_')
                all_geojsons = []
                for state_dir in STATES_DIR.glob(f"{state_dir_name}_*"):
                    all_geojsons.extend(state_dir.glob("*.geojson"))
                
                print(f"  Found {len(all_geojsons)} total geojson(s) in {state}")
                
                # Try each geojson to find matches
                for alt_geojson in all_geojsons:
                    if alt_geojson == geojson_path:
                        continue  # Skip the one we already tried
                    
                    print(f"\n  Trying: {alt_geojson.name}")
                    
                    # Load cities from this geojson
                    alt_cities = load_geojson_cities(alt_geojson)
                    if not alt_cities:
                        continue
                    
                    # Load polygons for these cities
                    alt_city_layer_map = load_all_city_polygons(state, alt_cities, layers)
                    if not alt_city_layer_map:
                        continue
                    
                    # Quick test: check if this geojson matches
                    test_found = False
                    test_city_name = None
                    tsv_first = tsv_parts[0]
                    try:
                        test_chunk = pd.read_csv(tsv_first, sep='\t', nrows=10000)
                        test_gdf = process_tsv_chunk(test_chunk)
                        
                        for (test_city, test_layer), test_polygons in alt_city_layer_map.items():
                            test_join = gpd.sjoin(test_gdf, test_polygons, how='inner', predicate='within')
                            if len(test_join) > 0:
                                print(f"    ✅ Found {len(test_join)} matches in {test_city}!")
                                test_found = True
                                test_city_name = test_city
                                break
                        
                        if test_found:
                            # Found the correct geojson! Now actually process it
                            correct_city_code = alt_geojson.stem
                            print(f"\n  🎯 CORRECT GEOJSON: {alt_geojson.name}")
                            print(f"  📝 Correct city code for this data: {correct_city_code}")
                            print(f"\n  ♻️  Reprocessing with correct geojson...")
                            
                            # Use the already-loaded polygons and process the chunks we already read
                            alt_city_layer_monthly = {key: [] for key in alt_city_layer_map.keys()}
                            
                            # Re-process all TSV parts with correct polygons
                            for tsv_part in tsv_parts:
                                print(f"    Processing {tsv_part.name}...")
                                try:
                                    for chunk_num, chunk in enumerate(pd.read_csv(tsv_part, sep='\t', chunksize=CHUNK_SIZE)):
                                        gdf_chunk = process_tsv_chunk(chunk)
                                        
                                        for (alt_city, alt_layer), alt_polygons in alt_city_layer_map.items():
                                            points_in_layer = gpd.sjoin(
                                                gdf_chunk,
                                                alt_polygons,
                                                how='inner',
                                                predicate='within'
                                            )
                                            
                                            if len(points_in_layer) > 0:
                                                if 'city' in points_in_layer.columns:
                                                    normalized_target = normalize_city_name(alt_city)
                                                    points_in_layer = points_in_layer[
                                                        points_in_layer['city'].apply(normalize_city_name) == normalized_target
                                                    ]
                                                
                                                if len(points_in_layer) > 0:
                                                    monthly = aggregate_chunk_monthly(points_in_layer, alt_layer)
                                                    alt_city_layer_monthly[(alt_city, alt_layer)].append(monthly)
                                        
                                        del gdf_chunk
                                except Exception as e:
                                    print(f"      Error: {e}")
                                    continue
                            
                            # Generate HTMLs for each city
                            for alt_city in alt_cities:
                                alt_city_visits = []
                                
                                for alt_layer in layers:
                                    key = (alt_city, alt_layer)
                                    if key not in alt_city_layer_map or not alt_city_layer_monthly[key]:
                                        continue
                                    
                                    all_monthly = pd.concat(alt_city_layer_monthly[key], ignore_index=True)
                                    final_monthly = all_monthly.groupby(['month', 'layer'])['visits'].sum().reset_index()
                                    
                                    total_visits = final_monthly['visits'].sum()
                                    if total_visits > 0:
                                        print(f"    {alt_city} / {alt_layer}: {total_visits:,} visits")
                                        alt_city_visits.append(final_monthly)
                                
                                if alt_city_visits:
                                    city_visits_df = pd.concat(alt_city_visits, ignore_index=True)
                                    
                                    fig = px.line(
                                        city_visits_df,
                                        x='month',
                                        y='visits',
                                        color='layer',
                                        line_dash='layer',
                                        title=f"{alt_city} - Monthly Visits"
                                    )
                                    
                                    html_path = OUT_DIR / f"{correct_city_code}_{alt_city.replace(' ', '_')}.html"
                                    fig.write_html(html_path)
                                    print(f"    ✅ Saved: {html_path.name}")
                            
                            break  # Found and processed, stop searching
                            
                    except Exception as e:
                        print(f"    Error testing: {e}")
                        continue
            
            # Free all city/layer data for this TSV
            del city_layer_monthly
        
        except Exception as e:
            print(f"  ❌ Error processing {tsv_file.name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*80}")
    print("✅ Processing complete!")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
