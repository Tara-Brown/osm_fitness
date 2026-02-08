#!/usr/bin/env python3
"""
Analyze visits with smoke data.

This script:
1. Reads visit CSV files (one per city-layer-day)
2. Merges with smoke lookup table (city-date-smoke)
3. Categorizes smoke levels
4. Generates summary statistics per city:
   - Total visits
   - Days per smoke category
   - Average visits per smoke category
   - City name

Output: city_smoke_visit_summary.csv
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from glob import glob
import warnings

warnings.filterwarnings('ignore')


# Smoke category thresholds (PM2.5 µg/m³)
# Adjust these based on your analysis needs
SMOKE_CATEGORIES = {
    'Good': (0, 12),
    'Moderate': (12, 35.5),
    'Unhealthy_Sensitive': (35.5, 55.5),
    'Unhealthy': (55.5, 150.5),
    'Very_Unhealthy': (150.5, 250.5),
    'Hazardous': (250.5, float('inf'))
}

# Alternative simpler categorization
SIMPLE_CATEGORIES = {
    'Low': (0, 10),
    'Medium': (10, 25),
    'High': (25, 50),
    'Very_High': (50, float('inf'))
}


def categorize_smoke(smoke_level, categories=SMOKE_CATEGORIES):
    """Assign smoke level to category."""
    for category, (low, high) in categories.items():
        if low <= smoke_level < high:
            return category
    return 'Unknown'


def load_visit_files(visit_dir, pattern='*_daily.csv'):
    """
    Load all visit CSV files from a directory.
    Returns a single DataFrame with all visits.
    """
    print(f"Loading visit files from {visit_dir}...")
    
    visit_files = list(Path(visit_dir).glob(pattern))
    
    if not visit_files:
        print(f"❌ No visit files found matching pattern: {pattern}")
        sys.exit(1)
    
    print(f"  Found {len(visit_files)} visit files")
    
    # Load and concatenate
    dfs = []
    for file in visit_files:
        try:
            df = pd.read_csv(file)
            dfs.append(df)
        except Exception as e:
            print(f"  ⚠️  Error loading {file.name}: {e}")
            continue
    
    if not dfs:
        print(f"❌ No files could be loaded")
        sys.exit(1)
    
    combined = pd.concat(dfs, ignore_index=True)
    
    print(f"  ✅ Loaded {len(combined):,} visit records")
    print(f"  Cities: {combined['city_code'].nunique() if 'city_code' in combined.columns else 'N/A'}")
    print(f"  Layers: {combined['layer'].nunique() if 'layer' in combined.columns else 'N/A'}")
    
    return combined


def load_smoke_lookup(smoke_csv):
    """Load smoke lookup table."""
    print(f"\nLoading smoke lookup from {smoke_csv}...")
    
    smoke_df = pd.read_csv(smoke_csv)
    
    print(f"  ✅ Loaded {len(smoke_df):,} smoke records")
    print(f"  Cities: {smoke_df['city_code'].nunique()}")
    print(f"  Date range: {smoke_df['date'].min()} to {smoke_df['date'].max()}")
    
    return smoke_df


def merge_visits_with_smoke(visits_df, smoke_df, smoke_column='weighted_avg_smoke_pm'):
    """
    Merge visit data with smoke data on city_code and date.
    """
    print(f"\nMerging visits with smoke data...")
    
    # Ensure date columns are in same format
    visits_df['date'] = pd.to_datetime(visits_df['date']).dt.strftime('%Y-%m-%d')
    smoke_df['date'] = pd.to_datetime(smoke_df['date']).dt.strftime('%Y-%m-%d')
    
    # Merge
    merged = visits_df.merge(
        smoke_df[['city_code', 'date', smoke_column, 'n_grids']],
        on=['city_code', 'date'],
        how='left'
    )
    
    # Report matching
    total_records = len(merged)
    matched_records = merged[smoke_column].notna().sum()
    match_rate = matched_records / total_records * 100
    
    print(f"  ✅ Merged {len(merged):,} records")
    print(f"  Matched with smoke: {matched_records:,} ({match_rate:.1f}%)")
    
    if match_rate < 50:
        print(f"  ⚠️  Warning: Low match rate! Check that city_code and date formats match")
    
    return merged


def generate_city_summaries(merged_df, smoke_column='weighted_avg_smoke_pm', 
                           categories=SMOKE_CATEGORIES):
    """
    Generate summary statistics per city.
    
    Returns DataFrame with:
    - city_code, city_name
    - total_visits
    - days_Good, days_Moderate, etc. (one per category)
    - avg_visits_Good, avg_visits_Moderate, etc.
    """
    print(f"\nGenerating city summaries...")
    
    # Add smoke category
    merged_df['smoke_category'] = merged_df[smoke_column].apply(
        lambda x: categorize_smoke(x, categories) if pd.notna(x) else 'Unknown'
    )
    
    # Remove records without smoke data
    analysis_df = merged_df[merged_df['smoke_category'] != 'Unknown'].copy()
    
    print(f"  Analyzing {len(analysis_df):,} records with smoke data")
    
    summaries = []
    
    for city_code in analysis_df['city_code'].unique():
        city_data = analysis_df[analysis_df['city_code'] == city_code]
        
        # Get city name
        city_name = city_data['city'].iloc[0] if 'city' in city_data.columns else str(city_code)
        
        # Total visits
        total_visits = city_data['unique_visits'].sum()
        
        # Days and average visits per category
        summary = {
            'city_code': city_code,
            'city_name': city_name,
            'total_visits': total_visits,
            'total_days': city_data['date'].nunique()
        }
        
        for category in categories.keys():
            cat_data = city_data[city_data['smoke_category'] == category]
            
            # Days in this category
            summary[f'days_{category}'] = cat_data['date'].nunique()
            
            # Average visits in this category
            if len(cat_data) > 0:
                summary[f'avg_visits_{category}'] = cat_data['unique_visits'].mean()
            else:
                summary[f'avg_visits_{category}'] = 0
        
        # Add overall statistics
        summary['avg_smoke_pm'] = city_data[smoke_column].mean()
        summary['max_smoke_pm'] = city_data[smoke_column].max()
        summary['min_smoke_pm'] = city_data[smoke_column].min()
        
        summaries.append(summary)
    
    summary_df = pd.DataFrame(summaries)
    summary_df = summary_df.sort_values('total_visits', ascending=False)
    
    print(f"  ✅ Generated summaries for {len(summary_df)} cities")
    
    return summary_df


def save_detailed_data(merged_df, output_file):
    """Save the full merged dataset for further analysis."""
    print(f"\nSaving detailed data to {output_file}...")
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(output_file, index=False)
    
    print(f"  ✅ Saved {len(merged_df):,} records")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze visit data with smoke levels',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python analyze_visits_smoke.py \\
        --visit-dir multilayer_outputs \\
        --smoke-lookup city_smoke_lookup.csv \\
        --output-summary city_smoke_summary.csv \\
        --output-detailed visits_with_smoke.csv

Smoke categories (EPA AQI-based):
    Good: 0-12 µg/m³
    Moderate: 12-35.5 µg/m³
    Unhealthy for Sensitive: 35.5-55.5 µg/m³
    Unhealthy: 55.5-150.5 µg/m³
    Very Unhealthy: 150.5-250.5 µg/m³
    Hazardous: 250.5+ µg/m³

Use --simple-categories for simpler grouping:
    Low: 0-10 µg/m³
    Medium: 10-25 µg/m³
    High: 25-50 µg/m³
    Very High: 50+ µg/m³
        """
    )
    
    parser.add_argument('--visit-dir', type=Path, required=True,
                       help='Directory with visit CSV files')
    parser.add_argument('--visit-pattern', type=str, default='*_daily.csv',
                       help='Pattern to match visit files (default: *_daily.csv)')
    parser.add_argument('--smoke-lookup', type=Path, required=True,
                       help='City-smoke lookup CSV file')
    parser.add_argument('--smoke-column', type=str, default='weighted_avg_smoke_pm',
                       help='Column to use for smoke level (default: weighted_avg_smoke_pm)')
    parser.add_argument('--output-summary', type=Path, required=True,
                       help='Output file for city summaries')
    parser.add_argument('--output-detailed', type=Path,
                       help='Optional: Save detailed merged data')
    parser.add_argument('--simple-categories', action='store_true',
                       help='Use simpler smoke categories (Low/Medium/High/Very_High)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.visit_dir.exists():
        print(f"❌ Visit directory not found: {args.visit_dir}")
        sys.exit(1)
    
    if not args.smoke_lookup.exists():
        print(f"❌ Smoke lookup file not found: {args.smoke_lookup}")
        sys.exit(1)
    
    # Choose category scheme
    categories = SIMPLE_CATEGORIES if args.simple_categories else SMOKE_CATEGORIES
    
    print(f"\n{'='*80}")
    print(f"Visit-Smoke Analysis")
    print(f"{'='*80}")
    print(f"Smoke categories: {list(categories.keys())}")
    print(f"{'='*80}\n")
    
    # Load data
    visits_df = load_visit_files(args.visit_dir, args.visit_pattern)
    smoke_df = load_smoke_lookup(args.smoke_lookup)
    
    # Merge
    merged_df = merge_visits_with_smoke(visits_df, smoke_df, args.smoke_column)
    
    # Generate summaries
    summary_df = generate_city_summaries(merged_df, args.smoke_column, categories)
    
    # Save summary
    print(f"\nSaving summary to {args.output_summary}...")
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.output_summary, index=False)
    print(f"  ✅ Saved summary")
    
    # Save detailed if requested
    if args.output_detailed:
        save_detailed_data(merged_df, args.output_detailed)
    
    # Display sample results
    print(f"\n{'='*80}")
    print(f"✅ ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nTop 10 cities by total visits:")
    print(summary_df.head(10)[['city_name', 'total_visits', 'total_days', 'avg_smoke_pm']].to_string())
    print(f"\nFull results saved to: {args.output_summary}")
    
    if args.output_detailed:
        print(f"Detailed data saved to: {args.output_detailed}")


if __name__ == "__main__":
    main()
