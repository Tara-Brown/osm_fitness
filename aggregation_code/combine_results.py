#!/usr/bin/env python3
"""
Combine and analyze multilayer processing results.

This script provides utilities to:
1. Combine all city parquet files into a single master file
2. Generate summary statistics including 95th and 98th percentiles
3. Create comparison visualizations

Usage:
    # Combine all results into a single master file (default)
    python combine_results.py --input-dir ./output --output-file combined_all_cities.parquet

    # Generate summary statistics with percentiles
    python combine_results.py --input-dir ./output --summarize --output-file summary_stats.parquet

    # Compare layer performance by city
    python combine_results.py --input-dir ./output --compare-layers --output-file layer_comparison.parquet
"""

import argparse
from pathlib import Path
import pandas as pd
import sys


def find_output_files(input_dir, pattern='daily_device_instances*'):
    """Find all output parquet files, excluding already-combined ALL_CITIES files."""
    all_files = list(Path(input_dir).glob(pattern))
    return [f for f in all_files if not f.name.startswith('ALL_CITIES_')]


def print_percentiles(series, label="unique_visits"):
    """Print the 95th and 98th percentiles for a given Series."""
    p95 = series.quantile(0.95)
    p98 = series.quantile(0.98)
    print(f"\n📊 Percentiles for {label}:")
    print(f"   95th percentile: {p95:,.2f}")
    print(f"   98th percentile: {p98:,.2f}")


def combine_all_files(input_dir, output_file):
    """
    Combine all individual city/layer parquet files into one master file.

    Args:
        input_dir: Directory containing individual parquet files
        output_file: Path for combined output file
    """
    files = find_output_files(input_dir)

    if not files:
        print(f"No parquet files found in {input_dir}")
        return

    print(f"Found {len(files)} parquet files to combine")

    dfs = []
    for file in files:
        try:
            df = pd.read_parquet(file)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file.name}: {e}")
            continue

    if not dfs:
        print("No data to combine")
        return

    combined = pd.concat(dfs, ignore_index=True)
    combined['date'] = pd.to_datetime(combined['date'])

    sort_cols = [c for c in ['city_name', 'layer', 'date'] if c in combined.columns]
    combined = combined.sort_values(sort_cols)

    output_file = Path(output_file).with_suffix('.parquet')
    combined.to_parquet(output_file, index=False, engine='pyarrow')

    print(f"\n✅ Combined {len(files)} files into {output_file}")
    print(f"   Total rows: {len(combined):,}")
    if 'city_name' in combined.columns:
        print(f"   Cities: {combined['city_name'].nunique()}")
    if 'layer' in combined.columns:
        print(f"   Layers: {combined['layer'].nunique()}")
    print(f"   Date range: {combined['date'].min()} to {combined['date'].max()}")

    if 'unique_visits' in combined.columns:
        print_percentiles(combined['unique_visits'], label="unique_visits (all cities/layers)")


def create_summary_statistics(input_dir, output_file):
    """
    Create summary statistics from all output files, including percentiles.

    Args:
        input_dir: Directory containing individual parquet files
        output_file: Path for summary statistics parquet
    """
    files = find_output_files(input_dir)

    if not files:
        print(f"No parquet files found in {input_dir}")
        return

    print(f"Analyzing {len(files)} parquet files...")

    summary_data = []
    all_visits = []

    for file in files:
        try:
            df = pd.read_parquet(file)
            df['date'] = pd.to_datetime(df['date'])

            if 'unique_visits' in df.columns:
                all_visits.extend(df['unique_visits'].dropna().tolist())

            city = df['city_name'].iloc[0] if 'city_name' in df.columns else file.stem
            layer = df['layer'].iloc[0] if 'layer' in df.columns else 'unknown'

            stats = {
                'city': city,
                'layer': layer,
                'source_file': file.name,
                'total_rows': len(df),
                'date_start': df['date'].min(),
                'date_end': df['date'].max(),
            }

            if 'unique_visits' in df.columns:
                stats.update({
                    'total_unique_visits': df['unique_visits'].sum(),
                    'avg_daily_visits': df['unique_visits'].mean(),
                    'median_daily_visits': df['unique_visits'].median(),
                    'p95_daily_visits': df['unique_visits'].quantile(0.95),
                    'p98_daily_visits': df['unique_visits'].quantile(0.98),
                    'max_daily_visits': df['unique_visits'].max(),
                    'min_daily_visits': df['unique_visits'].min(),
                    'std_daily_visits': df['unique_visits'].std(),
                })

            if 'feature_id' in df.columns:
                stats['num_polygons'] = df['feature_id'].nunique()

            summary_data.append(stats)

        except Exception as e:
            print(f"Error processing {file.name}: {e}")
            continue

    summary_df = pd.DataFrame(summary_data)
    if 'city' in summary_df.columns and 'layer' in summary_df.columns:
        summary_df = summary_df.sort_values(['city', 'layer'])

    output_file = Path(output_file).with_suffix('.parquet')
    summary_df.to_parquet(output_file, index=False, engine='pyarrow')

    print(f"\n✅ Summary statistics saved to {output_file}")
    print(f"   Total files summarized: {len(summary_df)}")

    if all_visits:
        print_percentiles(pd.Series(all_visits), label="unique_visits (across all files)")

    if 'total_unique_visits' in summary_df.columns:
        print(f"\n📊 Top 10 by Total Visits:")
        top_10 = summary_df.nlargest(10, 'total_unique_visits')[['city', 'layer', 'total_unique_visits']]
        print(top_10.to_string(index=False))


def compare_layers_by_city(input_dir, output_file):
    """
    Create a pivot comparing each layer's total visits by city.

    Args:
        input_dir: Directory containing individual parquet files
        output_file: Path for comparison parquet
    """
    files = find_output_files(input_dir)

    if not files:
        print(f"No parquet files found in {input_dir}")
        return

    print(f"Creating layer comparison from {len(files)} files...")

    dfs = []
    for file in files:
        try:
            dfs.append(pd.read_parquet(file))
        except Exception as e:
            print(f"Error reading {file.name}: {e}")
            continue

    if not dfs:
        print("No data to analyze")
        return

    combined = pd.concat(dfs, ignore_index=True)

    city_col = 'city_name' if 'city_name' in combined.columns else None
    if city_col is None or 'layer' not in combined.columns or 'unique_visits' not in combined.columns:
        print("❌ Combined data is missing required columns (city_name, layer, unique_visits)")
        return

    comparison = combined.groupby([city_col, 'layer']).agg(
        total_visits=('unique_visits', 'sum'),
        avg_daily_visits=('unique_visits', 'mean'),
        num_days=('date', 'nunique')
    ).reset_index()

    pivot = comparison.pivot(index=city_col, columns='layer', values='total_visits')
    pivot = pivot.fillna(0)
    pivot['total_all_layers'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('total_all_layers', ascending=False)

    output_file = Path(output_file).with_suffix('.parquet')
    pivot.to_parquet(output_file, engine='pyarrow')

    print(f"\n✅ Layer comparison saved to {output_file}")
    print(f"\nTop 10 cities by total visits across all layers:")
    print(pivot.head(10).to_string())

    print_percentiles(pivot['total_all_layers'], label="total_all_layers (per city)")


def main():
    parser = argparse.ArgumentParser(
        description='Combine and analyze multilayer processing results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--input-dir', type=Path, required=True,
                       help='Directory containing individual city parquet files')
    parser.add_argument('--output-file', type=Path,
                       help='Output file path for all operations')
    parser.add_argument('--combine-all', action='store_true',
                       help='Combine all files into one master file (default)')
    parser.add_argument('--summarize', action='store_true',
                       help='Generate summary statistics with percentiles')
    parser.add_argument('--compare-layers', action='store_true',
                       help='Compare layer performance by city')

    args = parser.parse_args()

    if not any([args.combine_all, args.summarize, args.compare_layers]):
        args.combine_all = True

    if args.combine_all:
        output = args.output_file or Path(args.input_dir) / 'combined_all_cities.parquet'
        combine_all_files(args.input_dir, output)

    if args.summarize:
        output = args.output_file or Path(args.input_dir) / 'summary_statistics.parquet'
        create_summary_statistics(args.input_dir, output)

    if args.compare_layers:
        output = args.output_file or Path(args.input_dir) / 'layer_comparison.parquet'
        compare_layers_by_city(args.input_dir, output)

    print("\n✅ All operations complete!")


if __name__ == '__main__':
    main()
