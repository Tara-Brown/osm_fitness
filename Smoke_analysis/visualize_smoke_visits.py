#!/usr/bin/env python3
"""
Generate visualizations from smoke-visit analysis results.

Creates publication-ready plots showing:
- Visits by smoke category
- City comparisons
- Time series with smoke overlays
- Correlation plots
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")


def plot_visits_by_smoke_category(summary_df, output_dir):
    """
    Bar plot showing average visits per smoke category.
    """
    print("Creating visits by smoke category plot...")
    
    # Get smoke categories from column names
    cat_cols = [col for col in summary_df.columns if col.startswith('avg_visits_')]
    categories = [col.replace('avg_visits_', '') for col in cat_cols]
    
    # Calculate overall averages
    avg_visits = summary_df[cat_cols].mean()
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    bars = ax.bar(categories, avg_visits, color=sns.color_palette("RdYlGn_r", len(categories)))
    
    ax.set_xlabel('Smoke Category', fontsize=12, fontweight='bold')
    ax.set_ylabel('Average Daily Visits', fontsize=12, fontweight='bold')
    ax.set_title('Average Park Visits by Smoke Level Category', fontsize=14, fontweight='bold')
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{height:.0f}',
               ha='center', va='bottom', fontweight='bold')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    output_file = output_dir / 'visits_by_smoke_category.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_city_comparison(summary_df, output_dir, top_n=10):
    """
    Compare top cities across smoke categories.
    """
    print(f"Creating top {top_n} cities comparison plot...")
    
    # Get top N cities by total visits
    top_cities = summary_df.nlargest(top_n, 'total_visits')
    
    # Get category columns
    cat_cols = [col for col in summary_df.columns if col.startswith('avg_visits_')]
    categories = [col.replace('avg_visits_', '') for col in cat_cols]
    
    # Prepare data
    city_names = top_cities['city_name'].values
    data = top_cities[cat_cols].values
    
    # Create grouped bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(city_names))
    width = 0.8 / len(categories)
    
    for i, (cat, col) in enumerate(zip(categories, cat_cols)):
        offset = width * i - (width * len(categories) / 2)
        ax.bar(x + offset, top_cities[col], width, label=cat)
    
    ax.set_xlabel('City', fontsize=12, fontweight='bold')
    ax.set_ylabel('Average Daily Visits', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Cities: Visits by Smoke Category', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(city_names, rotation=45, ha='right')
    ax.legend(title='Smoke Category', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    
    output_file = output_dir / f'top_{top_n}_cities_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_smoke_vs_visits_scatter(detailed_df, output_dir, sample_frac=0.1):
    """
    Scatter plot of smoke levels vs visits with trendline.
    """
    print("Creating smoke vs visits scatter plot...")
    
    # Sample data if too large
    if len(detailed_df) > 10000:
        plot_df = detailed_df.sample(frac=sample_frac, random_state=42)
    else:
        plot_df = detailed_df
    
    # Remove missing values
    plot_df = plot_df.dropna(subset=['weighted_avg_smoke_pm', 'unique_visits'])
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot
    ax.scatter(plot_df['weighted_avg_smoke_pm'], plot_df['unique_visits'], 
              alpha=0.3, s=20, color='steelblue')
    
    # Add trendline
    z = np.polyfit(plot_df['weighted_avg_smoke_pm'], plot_df['unique_visits'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(plot_df['weighted_avg_smoke_pm'].min(), 
                        plot_df['weighted_avg_smoke_pm'].max(), 100)
    ax.plot(x_line, p(x_line), "r--", linewidth=2, label=f'Trend: y={z[0]:.2f}x+{z[1]:.0f}')
    
    # Calculate correlation
    corr = plot_df[['weighted_avg_smoke_pm', 'unique_visits']].corr().iloc[0, 1]
    
    ax.set_xlabel('Smoke PM2.5 (µg/m³)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Unique Daily Visits', fontsize=12, fontweight='bold')
    ax.set_title(f'Smoke Level vs Park Visits\n(Correlation: {corr:.3f})', 
                fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = output_dir / 'smoke_vs_visits_scatter.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_time_series_with_smoke(detailed_df, output_dir, city_name=None, layer=None):
    """
    Time series plot showing visits and smoke levels over time.
    """
    print("Creating time series plot...")
    
    # Filter to specific city and layer if provided
    plot_df = detailed_df.copy()
    
    if city_name:
        plot_df = plot_df[plot_df['city'] == city_name]
    if layer:
        plot_df = plot_df[plot_df['layer'] == layer]
    
    if len(plot_df) == 0:
        print(f"  ⚠️  No data found for city={city_name}, layer={layer}")
        return
    
    # Convert date to datetime
    plot_df['date'] = pd.to_datetime(plot_df['date'])
    plot_df = plot_df.sort_values('date')
    
    # Aggregate by date (in case multiple layers)
    daily = plot_df.groupby('date').agg({
        'unique_visits': 'sum',
        'weighted_avg_smoke_pm': 'mean'
    }).reset_index()
    
    # Create plot with dual y-axis
    fig, ax1 = plt.subplots(figsize=(14, 6))
    
    # Visits on left axis
    color = 'tab:blue'
    ax1.set_xlabel('Date', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Unique Daily Visits', color=color, fontsize=12, fontweight='bold')
    ax1.plot(daily['date'], daily['unique_visits'], color=color, linewidth=1.5, label='Visits')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)
    
    # Smoke on right axis
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Smoke PM2.5 (µg/m³)', color=color, fontsize=12, fontweight='bold')
    ax2.plot(daily['date'], daily['weighted_avg_smoke_pm'], color=color, 
            linewidth=1.5, alpha=0.7, label='Smoke')
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Title
    title = 'Daily Visits and Smoke Levels Over Time'
    if city_name:
        title += f'\n{city_name}'
    if layer:
        title += f' - {layer}'
    ax1.set_title(title, fontsize=14, fontweight='bold')
    
    # Legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    suffix = ''
    if city_name:
        suffix += f'_{city_name.replace(" ", "_")}'
    if layer:
        suffix += f'_{layer}'
    
    output_file = output_dir / f'time_series{suffix}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_smoke_distribution(summary_df, output_dir):
    """
    Distribution of smoke levels across cities.
    """
    print("Creating smoke distribution plot...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    ax1.hist(summary_df['avg_smoke_pm'], bins=30, color='coral', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Average Smoke PM2.5 (µg/m³)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Cities', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of Average Smoke Levels', fontsize=13, fontweight='bold')
    ax1.axvline(summary_df['avg_smoke_pm'].mean(), color='red', linestyle='--', 
               linewidth=2, label=f"Mean: {summary_df['avg_smoke_pm'].mean():.2f}")
    ax1.axvline(summary_df['avg_smoke_pm'].median(), color='blue', linestyle='--', 
               linewidth=2, label=f"Median: {summary_df['avg_smoke_pm'].median():.2f}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    ax2.boxplot(summary_df['avg_smoke_pm'], vert=True)
    ax2.set_ylabel('Average Smoke PM2.5 (µg/m³)', fontsize=12, fontweight='bold')
    ax2.set_title('Smoke Level Box Plot', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = output_dir / 'smoke_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def create_summary_table(summary_df, output_dir):
    """
    Create a summary statistics table as an image.
    """
    print("Creating summary statistics table...")
    
    # Calculate statistics
    stats = {
        'Metric': [
            'Total Cities',
            'Total Visits',
            'Mean Visits per City',
            'Median Visits per City',
            'Mean Smoke PM2.5',
            'Median Smoke PM2.5',
            'Max Smoke PM2.5',
            'Min Smoke PM2.5'
        ],
        'Value': [
            f"{len(summary_df):,}",
            f"{summary_df['total_visits'].sum():,}",
            f"{summary_df['total_visits'].mean():,.0f}",
            f"{summary_df['total_visits'].median():,.0f}",
            f"{summary_df['avg_smoke_pm'].mean():.2f} µg/m³",
            f"{summary_df['avg_smoke_pm'].median():.2f} µg/m³",
            f"{summary_df['max_smoke_pm'].max():.2f} µg/m³",
            f"{summary_df['min_smoke_pm'].min():.2f} µg/m³"
        ]
    }
    
    stats_df = pd.DataFrame(stats)
    
    # Create table plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('tight')
    ax.axis('off')
    
    table = ax.table(cellText=stats_df.values, colLabels=stats_df.columns,
                    cellLoc='left', loc='center',
                    colWidths=[0.6, 0.4])
    
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2)
    
    # Style header
    for i in range(len(stats_df.columns)):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style rows
    for i in range(1, len(stats_df) + 1):
        if i % 2 == 0:
            table[(i, 0)].set_facecolor('#f0f0f0')
            table[(i, 1)].set_facecolor('#f0f0f0')
    
    plt.title('Summary Statistics', fontsize=14, fontweight='bold', pad=20)
    
    output_file = output_dir / 'summary_statistics_table.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Generate visualizations from smoke-visit analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--summary-file', type=Path, required=True,
                       help='City summary CSV file')
    parser.add_argument('--detailed-file', type=Path,
                       help='Detailed visits CSV file (optional, for time series)')
    parser.add_argument('--output-dir', type=Path, required=True,
                       help='Output directory for plots')
    parser.add_argument('--city', type=str,
                       help='Specific city for time series plot')
    parser.add_argument('--layer', type=str,
                       help='Specific layer for time series plot')
    parser.add_argument('--top-n', type=int, default=10,
                       help='Number of top cities to compare (default: 10)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.summary_file.exists():
        print(f"❌ Summary file not found: {args.summary_file}")
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*80}")
    print(f"Generating Visualizations")
    print(f"{'='*80}\n")
    
    # Load summary data
    print(f"Loading summary data from {args.summary_file}...")
    summary_df = pd.read_csv(args.summary_file)
    print(f"  ✅ Loaded {len(summary_df)} cities\n")
    
    # Generate plots
    plot_visits_by_smoke_category(summary_df, args.output_dir)
    plot_city_comparison(summary_df, args.output_dir, args.top_n)
    plot_smoke_distribution(summary_df, args.output_dir)
    create_summary_table(summary_df, args.output_dir)
    
    # Detailed plots if file provided
    if args.detailed_file and args.detailed_file.exists():
        print(f"\nLoading detailed data from {args.detailed_file}...")
        detailed_df = pd.read_csv(args.detailed_file)
        print(f"  ✅ Loaded {len(detailed_df):,} records\n")
        
        plot_smoke_vs_visits_scatter(detailed_df, args.output_dir)
        plot_time_series_with_smoke(detailed_df, args.output_dir, 
                                    args.city, args.layer)
    
    print(f"\n{'='*80}")
    print(f"✅ ALL VISUALIZATIONS COMPLETE")
    print(f"{'='*80}")
    print(f"\nOutput files saved to: {args.output_dir}")
    print(f"\nGenerated plots:")
    for plot_file in sorted(args.output_dir.glob('*.png')):
        print(f"  - {plot_file.name}")
    print()


if __name__ == "__main__":
    main()
