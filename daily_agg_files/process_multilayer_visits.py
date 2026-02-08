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
    'outdoor_sports_facilities',
    'indoor_sports_facilities'
]

# Processing parameters
CHUNK_SIZE = 100000
CRS = 'EPSG:4326'

# Directories
STATES_DIR = Path("/mnt/beegfs/projects/cellphone_mobility/states")
GPKG_DIR = Path("/mnt/beegfs/projects/cellphone_mobility/unified_geopackages_updated")

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
    'SD': 'South Dakota', 'TN': 'Tennessee', 'TX': 'Texas', 'UT': 'Utah',
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

def normalize_city_name(name):
    """Normalize city name for comparison."""
    return name.lower().strip().replace(' ', '_').replace('-', '_')


def parse_tsv_filename(filename):
    """
    Parse TSV filename to extract city name and state.
    
    Examples: 
        10056774_Arlington_WA_pin_report.tsv -> ('Arlington', 'WA')
        10056812_Seattle_1_WA_pin_report_005.tsv -> ('Seattle_1', 'WA')
    
    Returns: (city_name, state)
    """
    if isinstance(filename, Path):
        stem = filename.stem
    else:
        stem = Path(filename).stem
    
    # Remove _pin_report and any trailing numbers
    stem = re.sub(r'_pin_report(_\d{3})?$', '', stem)
    
    # Split by underscore
    parts = stem.split('_')
    
    if len(parts) < 3:
        raise ValueError(f"Cannot parse filename: {filename}")
    
    # First part is the numeric ID (skip it)
    # State is typically the last component (WA, CA, etc)
    state_abbrev = parts[-1]
    
    # City name is everything between the ID and state
    city_parts = parts[1:-1]
    city_name = '_'.join(city_parts)
    
    # Convert state abbreviation to full name
    if state_abbrev in STATE_ABBREV:
        state = STATE_ABBREV[state_abbrev]
    else:
        state = state_abbrev
    
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
                city_list = [c.strip() for c in all_cities_str.split(',')]
                cities.extend(city_list)
        
        return cities
    except Exception as e:
        print(f"  Error loading geojson: {e}")
        return []


def load_all_city_polygons(state, cities, layers, gpkg_path=None):
    """
    Load polygons for all city/layer combinations.
    Returns: {(city, layer): GeoDataFrame}
    """
    if gpkg_path is None:
        fips = STATE_FIPS.get(state)
        if not fips:
            print(f"  No FIPS code for {state}")
            return {}
        
        gpkg_path = GPKG_DIR / f"state_{fips}_unified.gpkg"
    
    if not gpkg_path.exists():
        print(f"  Geopackage not found: {gpkg_path}")
        return {}
    
    city_layer_map = {}
    
    for layer in layers:
        try:
            # Load layer
            gdf = gpd.read_file(gpkg_path, layer=layer)
            
            # **EXPLODE MULTIPOLYGONS HERE**
            original_count = len(gdf)
            gdf = gdf.explode(index_parts=False).reset_index(drop=True)
            exploded_count = len(gdf)
            
            if exploded_count > original_count:
                print(f"    Exploded {layer}: {original_count} → {exploded_count} polygons")
            
            # Recreate feature_id after exploding
            gdf['feature_id'] = range(len(gdf))
            
            # Check for city column (try both 'city' and 'city_name')
            city_col = None
            if 'city_name' in gdf.columns:
                city_col = 'city_name'
            else:
                print(f"  ⚠️  Layer {layer} has no 'city_name' column, skipping")
                continue
            
            # Filter to requested cities
            for city in cities:
                normalized_city = normalize_city_name(city)
                
                # Also try without number suffix (e.g., Seattle_2 -> Seattle)
                city_base = re.sub(r'_\d+$', '', city)  # Remove trailing _2, _1, etc.
                normalized_city_base = normalize_city_name(city_base)
                
                # Match by full name OR base name (without number)
                city_gdf = gdf[
                    (gdf[city_col].apply(normalize_city_name) == normalized_city) |
                    (gdf[city_col].apply(normalize_city_name) == normalized_city_base)
                ].copy()
                
                if len(city_gdf) > 0:
                    # feature_id already set above after exploding
                    city_layer_map[(city, layer)] = city_gdf
                    print(f"    Found {len(city_gdf)} features for {city} in {layer}")
        
        except Exception as e:
            print(f"  ⚠️  Error loading layer {layer}: {e}")
            continue
    
    return city_layer_map


def process_tsv_chunk(chunk):
    """Convert TSV chunk to GeoDataFrame."""
    # Parse timestamp - Unix timestamp format
    chunk['datetime'] = pd.to_datetime(chunk[TIME_COL], unit='s', errors='coerce')
    chunk['date'] = chunk['datetime'].dt.date
    
    # Create geometry
    geometry = [Point(xy) for xy in zip(chunk[LON_COL], chunk[LAT_COL])]
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
                geometry = [Point(xy) for xy in zip(complete_day[LON_COL], complete_day[LAT_COL])]
                gdf = gpd.GeoDataFrame(complete_day, geometry=geometry, crs=CRS)
                yield gdf
                
    except Exception as e:
        print(f"    Error in chunk_tsv_by_day: {e}")
        raise


def aggregate_daily_visits(points_gdf, layer_name, by_polygon=False):
    """Aggregate unique device visits by day."""
    if len(points_gdf) == 0:
        return pd.DataFrame()
    
    # Make sure we have the date column
    if 'date' not in points_gdf.columns:
        if 'datetime' in points_gdf.columns:
            points_gdf['date'] = points_gdf['datetime'].dt.date
        elif TIME_COL in points_gdf.columns:
            points_gdf['datetime'] = pd.to_datetime(points_gdf[TIME_COL], unit='s', errors='coerce')
            points_gdf['date'] = points_gdf['datetime'].dt.date
        else:
            return pd.DataFrame()
    
    if by_polygon and 'feature_id' in points_gdf.columns:
        # First compute per-device stats within each polygon-date
        device_stats = points_gdf.groupby(['date', 'feature_id', DEVICE_COL]).agg(
            ping_count=(DEVICE_COL, 'size'),
            time_range=('datetime', lambda x: (x.max() - x.min()).total_seconds() / 60 if len(x) > 1 else 0),
            avg_time_between_pings=('datetime', lambda x: x.diff().mean().total_seconds() / 60 if len(x) > 1 else None)
        ).reset_index()
        
        # Then aggregate to polygon-date level
        daily = device_stats.groupby(['date', 'feature_id']).agg(
            unique_visits=(DEVICE_COL, 'nunique'),
            total_pings=('ping_count', 'sum'),
            avg_pings_per_device=('ping_count', 'mean'),
            avg_time_range_minutes=('time_range', 'mean'),
            avg_time_between_pings_minutes=('avg_time_between_pings', 'mean')
        ).reset_index()
        
        daily['layer'] = layer_name
        daily.rename(columns={'feature_id': 'polygon_id'}, inplace=True)
    else:
        daily = points_gdf.groupby('date').agg(
            unique_visits=(DEVICE_COL, 'nunique')
        ).reset_index()
        daily['layer'] = layer_name
    
    return daily


# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_tsv_file(tsv_path, output_dir, by_polygon=False, gpkg_path=None):
    """Process a single TSV file with city-level filtering."""
    
    print(f"\n{'='*80}")
    print(f"Processing: {tsv_path.name}")
    print(f"{'='*80}")
    
    # Parse filename
    try:
        city_name, state = parse_tsv_filename(tsv_path)
        print(f"  City Name: {city_name}")
        print(f"  State: {state}")
    except Exception as e:
        print(f"  ❌ Error parsing filename: {e}")
        return {}
    
    # Find TSV parts
    tsv_parts = find_tsv_parts(tsv_path)
    print(f"  Found {len(tsv_parts)} TSV part(s)")
    
    # Find geojson
    geojson_path = find_geojson(city_name, state)
    if not geojson_path:
        print(f"  ❌ No geojson found, skipping")
        return {}
    
    print(f"  ✅ Found geojson: {geojson_path.name}")
    
    # Load cities from geojson
    cities = load_geojson_cities(geojson_path)
    if not cities:
        print(f"  ❌ No cities found in geojson")
        return {}
    
    print(f"  Found {len(cities)} cities: {', '.join(cities)}")
    
    # Extract city code from geojson filename
    city_code = geojson_path.stem
    print(f"  City Code: {city_code}")
    
    # Load city/layer polygons
    print(f"\n  📐 Loading city/layer polygons...")
    city_layer_map = load_all_city_polygons(state, cities, NEW_LAYERS, gpkg_path)
    
    if not city_layer_map:
        print(f"  ❌ No polygons found for any city/layer combination")
        return {}
    
    print(f"  ✅ Loaded polygons for {len(city_layer_map)} city/layer combinations")
    
    # Initialize storage
    daily_aggregations = {key: [] for key in city_layer_map.keys()}
    stats = {
        'total_chunks': 0,
        'total_points': 0,
        'points_per_layer': defaultdict(int),
        'city_code': city_code,
        'state': state,
        'city_name': city_name
    }
    
    # Process TSV parts
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
                
                # gdf_chunk is already processed (has geometry, date, datetime)
                
                # Match to each city/layer
                for (city, layer), polygons in city_layer_map.items():
                    try:
                        # Save date column before join (geopackage has its own 'date' column)
                        chunk_dates = gdf_chunk['date'].copy() if 'date' in gdf_chunk.columns else None
                        chunk_datetime = gdf_chunk['datetime'].copy() if 'datetime' in gdf_chunk.columns else None
                        
                        points_in_layer = gpd.sjoin(
                            gdf_chunk,
                            polygons,
                            how='inner',
                            predicate='within'
                        )
                        
                        # Restore our date columns (overwrite geopackage's date)
                        if chunk_dates is not None:
                            points_in_layer['date'] = chunk_dates.loc[points_in_layer.index]
                        if chunk_datetime is not None:
                            points_in_layer['datetime'] = chunk_datetime.loc[points_in_layer.index]
                        
                        if len(points_in_layer) > 0:
                            # Additional city filter if city column exists
                            if 'city' in points_in_layer.columns:
                                normalized_target = normalize_city_name(city)
                                points_in_layer = points_in_layer[
                                    points_in_layer['city'].apply(normalize_city_name) == normalized_target
                                ]
                            
                            if len(points_in_layer) > 0:
                                stats['points_per_layer'][layer] += len(points_in_layer)
                                
                                daily_chunk = aggregate_daily_visits(
                                    points_in_layer,
                                    layer,
                                    by_polygon=by_polygon
                                )
                                
                                if len(daily_chunk) > 0:
                                    daily_aggregations[(city, layer)].append(daily_chunk)
                    
                    except Exception as e:
                        print(f"    ⚠️  Error matching to {city}/{layer}: {e}")
                        continue
                
                del gdf_chunk
                
                if chunk_count % 10 == 0:
                    print(f"    Processed {chunk_count} day-chunks...")
        
        except Exception as e:
            print(f"  ⚠️  Error reading {tsv_part.name}: {e}")
            continue
    
    print(f"\n  ✅ Processed {stats['total_chunks']} chunks, {stats['total_points']:,} total points")
    
    if stats['points_per_layer']:
        print(f"\n  📍 Points matched per layer:")
        for layer, count in sorted(stats['points_per_layer'].items()):
            print(f"    {layer}: {count:,} points")
    else:
        print(f"\n  ⚠️  No points matched any layer!")
        return stats
    
    # Save results
    output_dir.mkdir(parents=True, exist_ok=True)
    saved_files = []
    
    for (city, layer), chunk_list in daily_aggregations.items():
        if not chunk_list:
            continue
        
        combined = pd.concat(chunk_list, ignore_index=True)
        
        #For multiple chunks, the means are somewhat off. For current purposes, I consider it sufficient (chunks are around the same size)
        if by_polygon:
            final = combined.groupby(['date', 'layer', 'polygon_id']).agg(
                unique_visits=('unique_visits', 'sum'),
                total_pings=('total_pings', 'sum'),
                avg_pings_per_device=('avg_pings_per_device', 'mean'),
                avg_time_range_minutes=('avg_time_range_minutes', 'mean'),
                avg_time_between_pings_minutes=('avg_time_between_pings_minutes', 'mean')
            ).reset_index()
        else:
            final = combined.groupby(['date', 'layer']).agg(
                unique_visits=('unique_visits', 'sum')
            ).reset_index()
        
        final['date'] = pd.to_datetime(final['date'])
        final = final.sort_values('date')
        final['city'] = city
        final['city_code'] = city_code
        
        suffix = '_by_polygon' if by_polygon else ''
        output_file = output_dir / f"{city_code}_{city.replace(' ', '_')}_{layer}{suffix}_daily.csv"
        final.to_csv(output_file, index=False)
        saved_files.append(output_file)
        
        print(f"  ✅ Saved {city}/{layer}: {output_file.name} ({len(final)} rows)")
    
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
    
    stats = process_tsv_file(
        args.tsv_file,
        args.output_dir,
        by_polygon=args.by_polygon,
        gpkg_path=args.gpkg_file
    )
    
    if stats.get('saved_files'):
        print(f"\n✅ Processing complete!")
        print(f"   Total points: {stats['total_points']:,}")
        print(f"   Files saved: {len(stats['saved_files'])}")
        sys.exit(0)
    else:
        print(f"\n⚠️  Processing complete but no output files created")
        sys.exit(1)


if __name__ == "__main__":
    main()
