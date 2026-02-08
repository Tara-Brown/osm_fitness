#!/usr/bin/env python3
"""
Create a lookup table linking cities to smoke levels by date.

This script:
1. Loads the 10km grid shapefile with grid IDs
2. Loads smoke data CSV (grid_id, date, smokePM_pred)
3. Loads city polygons from geopackage
4. Spatially joins cities to grids (weighted by overlap area)
5. Creates a city-date-smoke lookup CSV

Output: city_smoke_lookup.csv with columns:
    - city_code
    - city_name
    - date
    - avg_smoke_pm (area-weighted average across all overlapping grids)
    - max_smoke_pm (maximum smoke across overlapping grids)
    - min_smoke_pm (minimum smoke across overlapping grids)
    - n_grids (number of grids that overlap this city)
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import geopandas as gpd
import numpy as np
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')


def parse_smoke_date(date_str):
    """Convert YYYYMMDD string to YYYY-MM-DD format."""
    try:
        return datetime.strptime(str(date_str), '%Y%m%d').strftime('%Y-%m-%d')
    except:
        return str(date_str)


def load_smoke_data(smoke_csv):
    """Load and prepare smoke data."""
    print(f"Loading smoke data from {smoke_csv}...")
    
    df = pd.read_csv(smoke_csv)
    
    # Ensure column names are correct
    expected_cols = ['grid_id_10km', 'date', 'smokePM_pred']
    if not all(col in df.columns for col in expected_cols):
        print(f"❌ Error: Expected columns {expected_cols}, found {df.columns.tolist()}")
        sys.exit(1)
    
    # Convert date to standard format
    df['date'] = df['date'].apply(parse_smoke_date)
    
    print(f"  ✅ Loaded {len(df):,} smoke records")
    print(f"  Date range: {df['date'].min()} to {df['date'].max()}")
    print(f"  Unique grids: {df['grid_id_10km'].nunique():,}")
    
    return df


def load_grid_shapefile(grid_shp):
    """Load 10km grid shapefile."""
    print(f"\nLoading grid shapefile from {grid_shp}...")
    
    grid_gdf = gpd.read_file(grid_shp)
    
    # Ensure CRS is set
    if grid_gdf.crs is None:
        print("  ⚠️  No CRS found, assuming EPSG:4326")
        grid_gdf = grid_gdf.set_crs('EPSG:4326')
    
    print(f"  ✅ Loaded {len(grid_gdf):,} grid cells")
    print(f"  CRS: {grid_gdf.crs}")
    print(f"  Bounds: {grid_gdf.total_bounds}")
    
    return grid_gdf


def load_city_polygons(gpkg_path, layer_name='parks'):
    """
    Load city polygons from geopackage.
    Uses the first layer to get unique city polygons.
    """
    print(f"\nLoading city polygons from {gpkg_path}...")
    
    # Load the layer
    cities_gdf = gpd.read_file(gpkg_path, layer=layer_name)
    
    print(f"  ✅ Loaded {len(cities_gdf)} features from {layer_name} layer")
    
    # Check for city identifier columns
    city_cols = ['city', 'city_name']
    city_col = None
    for col in city_cols:
        if col in cities_gdf.columns:
            city_col = col
            break
    
    if city_col is None:
        print(f"  ⚠️  Warning: No city name column found in {cities_gdf.columns.tolist()}")
        print(f"  Will use index as city identifier")
        cities_gdf['city_name'] = cities_gdf.index.astype(str)
    else:
        cities_gdf['city_name'] = cities_gdf[city_col]
    
    print(f"  CRS: {cities_gdf.crs}")
    
    # Ensure same CRS as grids
    if cities_gdf.crs != 'EPSG:4326':
        print(f"  Reprojecting from {cities_gdf.crs} to EPSG:4326...")
        cities_gdf= cities_gdf.to_crs('EPSG:4326')
    
    return cities_gdf


def calculate_overlap_weights(cities_gdf, grid_gdf):
    """
    Calculate area-weighted overlaps between cities and grid cells.
    Returns a DataFrame with city-grid relationships and overlap weights.
    """
    print(f"\nCalculating city-grid overlaps...")
    print(f"  This may take a few minutes for large datasets...")
    
    # Ensure both are in same projected CRS for area calculations
    # Use equal area projection
    proj_crs = 'EPSG:5070'  # Albers Equal Area for continental US
    
    cities_proj = cities_gdf.to_crs(proj_crs)
    grid_proj = grid_gdf.to_crs(proj_crs)
    
    # Calculate city areas
    cities_proj['city_area'] = cities_proj.geometry.area
    
    # Spatial join
    overlaps = gpd.sjoin(cities_proj, grid_proj, how='inner', predicate='intersects')
    
    print(f"  ✅ Found {len(overlaps)} city-grid intersections")
    
    # Calculate intersection areas
    print(f"  Calculating intersection areas...")
    overlap_list = []
    
    for idx, row in overlaps.iterrows():
        city_geom = row.geometry
        grid_id = row['grid_id_10km']
        
        # Get the grid geometry
        grid_geom = grid_proj[grid_proj['grid_id_10km'] == grid_id].geometry.iloc[0]
        
        # Calculate intersection
        intersection = city_geom.intersection(grid_geom)
        intersection_area = intersection.area
        
        # Weight is intersection area / total city area
        weight = intersection_area / row['city_area']
        
        overlap_list.append({
            'city_name': row['city_name'],
            'grid_id_10km': grid_id,
            'overlap_weight': weight
        })
    
    overlap_df = pd.DataFrame(overlap_list)
    
    print(f"  ✅ Calculated weights for {len(overlap_df)} overlaps")
    
    return overlap_df


def create_city_smoke_lookup(smoke_df, overlap_df):
    """
    Create final lookup table with area-weighted smoke values by city and date.
    """
    print(f"\nCreating city-smoke lookup table...")
    
    # Merge smoke data with overlaps
    merged = overlap_df.merge(smoke_df, on='grid_id_10km', how='inner')
    
    print(f"  ✅ Merged {len(merged):,} city-grid-date records")
    
    # Calculate weighted average smoke per city-date
    print(f"  Calculating weighted averages...")
    
    # Weighted average: sum(smoke * weight) / sum(weight)
    grouped = merged.groupby(['city_code', 'city_name', 'date']).agg({
        'smokePM_pred': ['mean', 'min', 'max'],  # Basic stats
        'overlap_weight': 'sum',  # Total weight
        'grid_id_10km': 'count'  # Number of grids
    }).reset_index()
    
    # Flatten column names
    grouped.columns = ['city_name', 'date', 
                      'avg_smoke_pm', 'min_smoke_pm', 'max_smoke_pm',
                      'total_weight', 'n_grids']
    
    # Calculate truly weighted average
    merged['weighted_smoke'] = merged['smokePM_pred'] * merged['overlap_weight']
    weighted_avg = merged.groupby(['city_name', 'date'])['weighted_smoke'].sum().reset_index()
    weighted_avg = weighted_avg.merge(
        grouped[[ 'city_name', 'date', 'total_weight']], 
        on=['city_name', 'date']
    )
    weighted_avg['weighted_avg_smoke_pm'] = weighted_avg['weighted_smoke'] / weighted_avg['total_weight']
    
    # Merge back
    lookup = grouped.merge(
        weighted_avg[['city_name', 'date', 'weighted_avg_smoke_pm']],
        on=['city_name', 'date']
    )
    
    # Reorder columns
    lookup = lookup[[
         'city_name', 'date',
        'weighted_avg_smoke_pm', 'avg_smoke_pm', 'min_smoke_pm', 'max_smoke_pm',
        'n_grids'
    ]]
    
    lookup = lookup.sort_values(['city_name', 'date'])
    
    print(f"  ✅ Created lookup with {len(lookup):,} city-date records")
    print(f"  Cities: {lookup['city_name'].nunique()}")
    print(f"  Date range: {lookup['date'].min()} to {lookup['date'].max()}")
    
    return lookup


def main():
    parser = argparse.ArgumentParser(
        description='Create city-smoke lookup table from grid data and city polygons',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python create_city_smoke_lookup.py \\
        --smoke-csv smoke_data.csv \\
        --grid-shp grid_10km.shp \\
        --city-gpkg cities.gpkg \\
        --output city_smoke_lookup.csv
        """
    )
    
    parser.add_argument('--smoke-csv', type=Path, required=True,
                       help='CSV file with smoke data (grid_id_10km, date, smokePM_pred)')
    parser.add_argument('--grid-shp', type=Path, required=True,
                       help='Shapefile with 10km grid (must have grid_id column)')
    parser.add_argument('--city-gpkg', type=Path, required=True,
                       help='GeoPackage with city polygons')
    parser.add_argument('--city-layer', type=str, default='parks',
                       help='Layer name in geopackage to use for cities (default: parks)')
    parser.add_argument('--output', type=Path, required=True,
                       help='Output CSV file')
    parser.add_argument('--save-overlaps', type=Path,
                       help='Optional: Save city-grid overlap weights to CSV')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.smoke_csv.exists():
        print(f"❌ Smoke CSV not found: {args.smoke_csv}")
        sys.exit(1)
    
    if not args.grid_shp.exists():
        print(f"❌ Grid shapefile not found: {args.grid_shp}")
        sys.exit(1)
    
    if not args.city_gpkg.exists():
        print(f"❌ City geopackage not found: {args.city_gpkg}")
        sys.exit(1)
    
    # Load data
    smoke_df = load_smoke_data(args.smoke_csv)
    grid_gdf = load_grid_shapefile(args.grid_shp)
    cities_gdf = load_city_polygons(args.city_gpkg, args.city_layer)
    
    # Calculate overlaps
    overlap_df = calculate_overlap_weights(cities_gdf, grid_gdf)
    
    if args.save_overlaps:
        print(f"\nSaving overlap weights to {args.save_overlaps}...")
        overlap_df.to_csv(args.save_overlaps, index=False)
        print(f"  ✅ Saved overlap weights")
    
    # Create lookup table
    lookup_df = create_city_smoke_lookup(smoke_df, overlap_df)
    
    # Save output
    print(f"\nSaving lookup table to {args.output}...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    lookup_df.to_csv(args.output, index=False)
    
    print(f"\n{'='*80}")
    print(f"✅ SUCCESS!")
    print(f"{'='*80}")
    print(f"Output saved to: {args.output}")
    print(f"Rows: {len(lookup_df):,}")
    print(f"Cities: {lookup_df['city_code'].nunique()}")
    print(f"Dates: {lookup_df['date'].nunique()}")
    print(f"\nSample of output:")
    print(lookup_df.head(10).to_string())


if __name__ == "__main__":
    main()
