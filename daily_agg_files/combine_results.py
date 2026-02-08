#!/usr/bin/env python3
"""
Combine and analyze multilayer processing results.

This script provides utilities to:
1. Combine individual city CSV files into layer-specific or master files
2. Generate summary statistics
3. Create comparison visualizations

Usage:
    # Combine all results into layer-specific files
    python combine_results.py --input-dir ./output --output-file combined_by_layer.csv --group-by layer
    
    # Create master file with all data
    python combine_results.py --input-dir ./output --output-file master_all_cities.csv
    
    # Generate summary statistics
    python combine_results.py --input-dir ./output --summarize --output-file summary_stats.csv
"""

import argparse
from pathlib import Path
import pandas as pd
import sys

def find_output_files(input_dir, pattern='*_daily.csv'):
    """Find all output CSV files."""
    return list(Path(input_dir).glob(pattern))


def combine_all_files(input_dir, output_file):
    """
    Combine all individual city/layer CSV files into one master file.
    
    Args:
        input_dir: Directory containing individual CSV files
        output_file: Path for combined output file
    """
    files = find_output_files(input_dir)
    
    if not files:
        print(f"No CSV files found in {input_dir}")
        return
    
    print(f"Found {len(files)} CSV files to combine")
    
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file.name}: {e}")
            continue
    
    if not dfs:
        print("No data to combine")
        return
    
    # Combine all dataframes
    combined = pd.concat(dfs, ignore_index=True)
    
    # Sort by city, layer, and date
    combined['date'] = pd.to_datetime(combined['date'])
    combined = combined.sort_values(['city', 'layer', 'date'])
    
    # Save
    combined.to_csv(output_file, index=False)
    print(f"\n✅ Combined {len(files)} files into {output_file}")
    print(f"   Total rows: {len(combined):,}")
    print(f"   Cities: {combined['city'].nunique()}")
    print(f"   Layers: {combined['layer'].nunique()}")
    print(f"   Date range: {combined['date'].min()} to {combined['date'].max()}")


def combine_by_layer(input_dir, output_dir):
    """
    Create separate combined files for each layer.
    
    Args:
        input_dir: Directory containing individual CSV files
        output_dir: Directory for layer-specific output files
    """
    files = find_output_files(input_dir)
    
    if not files:
        print(f"No CSV files found in {input_dir}")
        return
    
    print(f"Found {len(files)} CSV files")
    
    # Group files by layer
    layer_files = {}
    for file in files:
        # Extract layer name from filename
        # Format: {code}_{city}_{layer}_daily.csv
        parts = file.stem.split('_')
        
        # Find the layer (it's the part before '_daily' or '_by_polygon')
        if '_by_polygon_daily' in file.stem:
            # Handle polygon-level files
            layer_idx = file.stem.rindex('_by_polygon_daily')
            layer_part = file.stem[:layer_idx]
            layer = '_'.join(layer_part.split('_')[2:])  # Skip code and city
        else:
            layer_idx = file.stem.rindex('_daily')
            layer_part = file.stem[:layer_idx]
            layer = '_'.join(layer_part.split('_')[2:])  # Skip code and city
        
        if layer not in layer_files:
            layer_files[layer] = []
        layer_files[layer].append(file)
    
    # Combine files for each layer
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for layer, layer_file_list in layer_files.items():
        print(f"\nProcessing layer: {layer}")
        print(f"  Files: {len(layer_file_list)}")
        
        dfs = []
        for file in layer_file_list:
            try:
                df = pd.read_csv(file)
                dfs.append(df)
            except Exception as e:
                print(f"  Error reading {file.name}: {e}")
                continue
        
        if dfs:
            combined = pd.concat(dfs, ignore_index=True)
            combined['date'] = pd.to_datetime(combined['date'])
            combined = combined.sort_values(['city', 'date'])
            
            output_file = output_dir / f"combined_{layer}.csv"
            combined.to_csv(output_file, index=False)
            print(f"  ✅ Saved: {output_file.name} ({len(combined):,} rows)")


def create_summary_statistics(input_dir, output_file):
    """
    Create summary statistics from all output files.
    
    Args:
        input_dir: Directory containing individual CSV files
        output_file: Path for summary statistics CSV
    """
    files = find_output_files(input_dir)
    
    if not files:
        print(f"No CSV files found in {input_dir}")
        return
    
    print(f"Analyzing {len(files)} CSV files...")
    
    summary_data = []
    
    for file in files:
        try:
            df = pd.read_csv(file)
            df['date'] = pd.to_datetime(df['date'])
            
            # Extract metadata
            city = df['city'].iloc[0] if 'city' in df.columns else 'Unknown'
            city_code = df['city_code'].iloc[0] if 'city_code' in df.columns else 'Unknown'
            layer = df['layer'].iloc[0] if 'layer' in df.columns else 'Unknown'
            
            # Calculate statistics
            stats = {
                'city_code': city_code,
                'city': city,
                'layer': layer,
                'total_days': len(df),
                'date_start': df['date'].min(),
                'date_end': df['date'].max(),
                'total_unique_visits': df['unique_visits'].sum(),
                'avg_daily_visits': df['unique_visits'].mean(),
                'median_daily_visits': df['unique_visits'].median(),
                'max_daily_visits': df['unique_visits'].max(),
                'min_daily_visits': df['unique_visits'].min(),
                'std_daily_visits': df['unique_visits'].std(),
            }
            
            # If polygon-level data
            if 'polygon_id' in df.columns:
                stats['num_polygons'] = df['polygon_id'].nunique()
                stats['avg_visits_per_polygon'] = df.groupby('polygon_id')['unique_visits'].sum().mean()
            
            summary_data.append(stats)
            
        except Exception as e:
            print(f"Error processing {file.name}: {e}")
            continue
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values(['city', 'layer'])
    
    # Save
    summary_df.to_csv(output_file, index=False)
    print(f"\n✅ Summary statistics saved to {output_file}")
    print(f"   Total city-layer combinations: {len(summary_df)}")
    
    # Print top performers
    print(f"\n📊 Top 10 City-Layer Combinations by Total Visits:")
    top_10 = summary_df.nlargest(10, 'total_unique_visits')[['city', 'layer', 'total_unique_visits']]
    print(top_10.to_string(index=False))


def compare_layers_by_city(input_dir, output_file):
    """
    Create a comparison showing each layer's performance by city.
    
    Args:
        input_dir: Directory containing individual CSV files
        output_file: Path for comparison CSV
    """
    files = find_output_files(input_dir)
    
    if not files:
        print(f"No CSV files found in {input_dir}")
        return
    
    print(f"Creating layer comparison from {len(files)} files...")
    
    # Load all data
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file.name}: {e}")
            continue
    
    if not dfs:
        print("No data to analyze")
        return
    
    combined = pd.concat(dfs, ignore_index=True)
    
    # Aggregate by city and layer
    comparison = combined.groupby(['city', 'layer']).agg(
        total_visits=('unique_visits', 'sum'),
        avg_daily_visits=('unique_visits', 'mean'),
        num_days=('date', 'nunique')
    ).reset_index()
    
    # Pivot for easier comparison
    pivot = comparison.pivot(index='city', columns='layer', values='total_visits')
    pivot = pivot.fillna(0)
    pivot['total_all_layers'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('total_all_layers', ascending=False)
    
    # Save
    pivot.to_csv(output_file)
    print(f"\n✅ Layer comparison saved to {output_file}")
    print(f"\nTop 10 cities by total visits across all layers:")
    print(pivot.head(10).to_string())


def main():
    parser = argparse.ArgumentParser(
        description='Combine and analyze multilayer processing results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input-dir', type=Path, required=True,
                       help='Directory containing individual city CSV files')
    parser.add_argument('--output-file', type=Path,
                       help='Output file path (for single-file operations)')
    parser.add_argument('--output-dir', type=Path,
                       help='Output directory (for multi-file operations)')
    
    # Operation modes
    parser.add_argument('--combine-all', action='store_true',
                       help='Combine all files into one master file')
    parser.add_argument('--combine-by-layer', action='store_true',
                       help='Create separate files for each layer')
    parser.add_argument('--summarize', action='store_true',
                       help='Generate summary statistics')
    parser.add_argument('--compare-layers', action='store_true',
                       help='Compare layer performance by city')
    
    args = parser.parse_args()
    
    # Default: combine all if no specific operation specified
    if not any([args.combine_all, args.combine_by_layer, args.summarize, args.compare_layers]):
        args.combine_all = True
    
    # Execute requested operations
    if args.combine_all:
        if not args.output_file:
            args.output_file = Path(args.input_dir) / 'combined_all_cities.csv'
        combine_all_files(args.input_dir, args.output_file)
    
    if args.combine_by_layer:
        if not args.output_dir:
            args.output_dir = Path(args.input_dir) / 'combined_by_layer'
        combine_by_layer(args.input_dir, args.output_dir)
    
    if args.summarize:
        if not args.output_file:
            args.output_file = Path(args.input_dir) / 'summary_statistics.csv'
        create_summary_statistics(args.input_dir, args.output_file)
    
    if args.compare_layers:
        if not args.output_file:
            args.output_file = Path(args.input_dir) / 'layer_comparison.csv'
        compare_layers_by_city(args.input_dir, args.output_file)
    
    print("\n✅ All operations complete!")


if __name__ == '__main__':
    main()
