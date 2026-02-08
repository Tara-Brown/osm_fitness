#!/usr/bin/env python3
"""
Process a single TSV file (wrapper around main script for parallel processing)
Usage: python process_single_tsv.py <tsv_file_path>
"""

import sys
from pathlib import Path
import re

if len(sys.argv) != 2:
    print("Usage: python process_single_tsv.py <tsv_file_path>")
    sys.exit(1)

tsv_file = Path(sys.argv[1])

# Handle case where file doesn't exist but parts might (e.g., given "file.tsv" but "file_000.tsv" exists)
if not tsv_file.exists():
    # Try to find parts
    parent = tsv_file.parent
    stem = tsv_file.stem
    
    # Check if this is a base name without _pin_report or without _XXX
    # Try finding actual files
    possible_patterns = [
        f"{stem}_000.tsv",  # Has _pin_report but no number
        f"{stem}_pin_report_000.tsv",  # Missing _pin_report entirely
    ]
    
    found = None
    for pattern in possible_patterns:
        candidate = parent / pattern
        if candidate.exists():
            found = candidate
            break
    
    # If still not found, try globbing
    if not found:
        # Try to find any file matching the base pattern
        matches = list(parent.glob(f"{stem.split('_pin_report')[0]}_pin_report*.tsv"))
        if matches:
            found = sorted(matches)[0]  # Take the first one
    
    if found:
        print(f"Note: {tsv_file.name} not found, using {found.name} instead")
        tsv_file = found
    else:
        print(f"Error: File not found: {tsv_file}")
        print(f"  Tried patterns: {possible_patterns}")
        sys.exit(1)

print(f"Processing single TSV: {tsv_file}")
print(f"File exists: {tsv_file.exists()}")
print(f"File size: {tsv_file.stat().st_size / (1024*1024):.2f} MB")
sys.stdout.flush()  # Force output to appear immediately

# Import all the functions and constants from the main script
print("Importing modules...")
sys.stdout.flush()
from process_city_visits import (
    STATE_FIPS, STATE_ABBREV, TSV_DIR, STATES_DIR, GPKG_DIR, OUT_DIR,
    CHUNK_SIZE, TIME_COL, LON_COL, LAT_COL,
    normalize_city_name, parse_tsv_filename, find_tsv_parts,
    find_geojson, load_geojson_cities, load_geopackage_layers,
    process_tsv_chunk, aggregate_monthly, aggregate_chunk_monthly,
    load_all_city_polygons
)
import pandas as pd
import geopandas as gpd
import plotly.express as px

# Create output directory
OUT_DIR.mkdir(exist_ok=True)

print(f"\n{'='*80}")
print(f"Processing: {tsv_file.name}")
print(f"{'='*80}")
sys.stdout.flush()

try:
    # Parse filename
    print("Parsing filename...")
    sys.stdout.flush()
    city_code, state = parse_tsv_filename(tsv_file)
    print(f"  City Code: {city_code}")
    print(f"  State: {state}")
    sys.stdout.flush()
    
    # Find ALL parts of this TSV (to handle multi-part files)
    tsv_parts = find_tsv_parts(tsv_file)
    print(f"  Found {len(tsv_parts)} TSV part(s)")
    
    # Find corresponding geojson
    geojson_path = find_geojson(city_code, state)
    if not geojson_path:
        sys.exit(1)
    
    print(f"  ✅ Found geojson: {geojson_path}")
    
    # Load cities from geojson
    cities = load_geojson_cities(geojson_path)
    print(f"  Found {len(cities)} cities: {', '.join(cities)}")
    
    # Get geopackage layers
    fips = STATE_FIPS.get(state)
    if not fips:
        print(f"  ❌ No FIPS code for {state}")
        sys.exit(1)
    
    gpkg_path = GPKG_DIR / f"state_{fips}_unified.gpkg"
    if not gpkg_path.exists():
        print(f"  ❌ Geopackage not found: {gpkg_path}")
        sys.exit(1)
    
    # Get available layers
    layers = ['unified_park_area', 'parkserve_geom', 'osm_geom', 'indoor_gyms']
    
    # Load all city/layer polygons once
    print(f"\n  📐 Loading all city/layer polygons...")
    city_layer_map = load_all_city_polygons(state, cities, layers)
    
    if not city_layer_map:
        print(f"  ❌ No polygons found for any city/layer combination")
        sys.exit(1)
    
    print(f"  ✅ Loaded polygons for {len(city_layer_map)} city/layer combinations")
    
    # Initialize storage for monthly aggregations
    city_layer_monthly = {key: [] for key in city_layer_map.keys()}
    
    # Process TSV files chunk by chunk, aggregating immediately
    print(f"\n  📊 Processing TSV files in chunks...")
    total_chunks = 0
    chunk_point_counts = {}
    
    for tsv_part in tsv_parts:
        print(f"    Reading {tsv_part.name}...")
        
        # Validate file
        if not tsv_part.exists():
            print(f"    ⚠️  File does not exist: {tsv_part}")
            continue
        
        try:
            file_size = tsv_part.stat().st_size
            print(f"      File stat size: {file_size / (1024*1024):.2f} MB")
        except Exception as e:
            print(f"    ⚠️  Could not stat file: {e}")
            file_size = 0
        
        try:
            with open(tsv_part, 'r') as f:
                first_line = f.readline().strip()
                if not first_line:
                    print(f"    ⚠️  File appears to be empty (no content)")
                    continue
                
                print(f"      First line preview: {first_line[:100]}...")
                
                if '\t' not in first_line:
                    print(f"    ⚠️  File doesn't appear to be tab-separated")
                    continue
                
                columns = first_line.split('\t')
                print(f"      Found {len(columns)} columns")
                if TIME_COL not in columns or LON_COL not in columns or LAT_COL not in columns:
                    print(f"    ⚠️  Missing required columns. Expected: {TIME_COL}, {LON_COL}, {LAT_COL}")
                    print(f"    ⚠️  Found columns: {columns[:10]}")
                    continue
                
                second_line = f.readline().strip()
                if not second_line:
                    print(f"    ⚠️  File has header but no data rows")
                    continue
                
                print(f"      ✅ File validation passed")
                
        except Exception as e:
            print(f"    ⚠️  Error validating file: {e}")
            continue
        
        try:
            for chunk_num, chunk in enumerate(pd.read_csv(tsv_part, sep='\t', chunksize=CHUNK_SIZE)):
                total_chunks += 1
                gdf_chunk = process_tsv_chunk(chunk)
                
                # For each city/layer combination, filter and aggregate this chunk
                for (city, layer), polygons in city_layer_map.items():
                    points_in_layer = gpd.sjoin(
                        gdf_chunk,
                        polygons,
                        how='inner',
                        predicate='within'
                    )
                    
                    if len(points_in_layer) > 0:
                        key = (city, layer)
                        if key not in chunk_point_counts:
                            chunk_point_counts[key] = 0
                        
                        before_filter = len(points_in_layer)
                        
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
                            monthly = aggregate_chunk_monthly(points_in_layer, layer)
                            city_layer_monthly[(city, layer)].append(monthly)
                
                del gdf_chunk
                
                if (chunk_num + 1) % 10 == 0:
                    print(f"      Processed {chunk_num + 1} chunks...")
                
        except Exception as e:
            print(f"    ⚠️  Error reading {tsv_part.name}: {e}")
            continue
    
    print(f"  ✅ Processed {total_chunks} total chunks")
    
    if chunk_point_counts:
        print(f"\n  📍 Points found per city/layer:")
        for (city, layer), count in sorted(chunk_point_counts.items()):
            print(f"    {city} / {layer}: {count:,} points")
    else:
        print(f"\n  ⚠️  No points found in any city/layer combination!")
    
    # Sum monthly aggregations for each city
    cities_with_data = []
    
    for city in cities:
        print(f"\n  📍 Processing city: {city}")
        city_visits = []
        
        for layer in layers:
            key = (city, layer)
            
            if key not in city_layer_map:
                continue
            
            if not city_layer_monthly[key]:
                continue
            
            all_monthly = pd.concat(city_layer_monthly[key], ignore_index=True)
            final_monthly = all_monthly.groupby(['month', 'layer'])['visits'].sum().reset_index()
            
            total_visits = final_monthly['visits'].sum()
            print(f"      {layer}: {total_visits:,} visits")
            
            city_visits.append(final_monthly)
        
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
    
    # FALLBACK: Search other geojsons if no matches
    if not cities_with_data and not chunk_point_counts:
        print(f"\n  ⚠️  No matches found with {geojson_path.name}")
        print(f"  🔍 Searching for matching geojson in other {state} cities...")
        
        state_dir_name = state.replace(' ', '_')
        all_geojsons = []
        for state_dir in STATES_DIR.glob(f"{state_dir_name}_*"):
            all_geojsons.extend(state_dir.glob("*.geojson"))
        
        print(f"  Found {len(all_geojsons)} total geojson(s) in {state}")
        
        for alt_geojson in all_geojsons:
            if alt_geojson == geojson_path:
                continue
            
            print(f"\n  Trying: {alt_geojson.name}")
            
            alt_cities = load_geojson_cities(alt_geojson)
            if not alt_cities:
                continue
            
            alt_city_layer_map = load_all_city_polygons(state, alt_cities, layers)
            if not alt_city_layer_map:
                continue
            
            test_found = False
            tsv_first = tsv_parts[0]
            try:
                test_chunk = pd.read_csv(tsv_first, sep='\t', nrows=10000)
                test_gdf = process_tsv_chunk(test_chunk)
                
                for (test_city, test_layer), test_polygons in alt_city_layer_map.items():
                    test_join = gpd.sjoin(test_gdf, test_polygons, how='inner', predicate='within')
                    if len(test_join) > 0:
                        print(f"    ✅ Found {len(test_join)} matches in {test_city}!")
                        test_found = True
                        break
                
                if test_found:
                    correct_city_code = alt_geojson.stem
                    print(f"\n  🎯 CORRECT GEOJSON: {alt_geojson.name}")
                    print(f"  📝 Correct city code for this data: {correct_city_code}")
                    print(f"\n  ♻️  Reprocessing with correct geojson...")
                    
                    # Use the already-loaded polygons and process all TSV parts
                    alt_city_layer_monthly = {key: [] for key in alt_city_layer_map.keys()}
                    
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
    
    print(f"\n{'='*80}")
    print("✅ Processing complete!")
    print(f"{'='*80}")

except Exception as e:
    print(f"  ❌ Error processing {tsv_file.name}: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
