#!/bin/bash
#
# Complete Smoke-Visit Analysis Workflow
# 
# This script runs all three steps of the analysis:
# 1. Create city-smoke lookup table
# 2. Check visit data (assumed already processed)
# 3. Analyze visits with smoke
#
# Usage: bash run_smoke_visit_analysis.sh

set -e  # Exit on error

# =============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# =============================================================================

# Input data
SMOKE_CSV="/mnt/beegfs/projects/cellphone_mobility/testing_files/smokePM2pt5_predictions_daily_10km_20060101-20231231.csv"
GRID_SHP="/mnt/beegfs/projects/cellphone_mobility/testing_files/smoke_grid_10km.shp"
CITY_GPKG="/mnt/beegfs/projects/cellphone_mobility/testing_files/geopackage_layers/state_53_components.gpkg"
VISIT_DIR="/mnt/beegfs/projects/cellphone_mobility/multilayer_outputs"

# Output directory
OUTPUT_DIR="/mnt/beegfs/projects/cellphone_mobility/smoke_analysis"

# Processing options
CITY_LAYER="parks"  # Layer in geopackage to use for city boundaries
SIMPLE_CATEGORIES=""  # Set to "--simple-categories" for simplified smoke categories

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Change to output directory
cd "$OUTPUT_DIR"

# =============================================================================
# STEP 1: CREATE CITY-SMOKE LOOKUP
# =============================================================================

echo ""
echo "=============================================================================="
echo "STEP 1: Creating City-Smoke Lookup Table"
echo "=============================================================================="
echo ""

LOOKUP_FILE="${OUTPUT_DIR}/city_smoke_lookup.csv"
OVERLAP_FILE="${OUTPUT_DIR}/city_grid_overlaps.csv"

if [ -f "$LOOKUP_FILE" ]; then
    echo "⚠️  Lookup file already exists: $LOOKUP_FILE"
    read -p "Recreate it? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping Step 1, using existing lookup file"
    else
        python create_city_smoke_lookup.py \
            --smoke-csv "$SMOKE_CSV" \
            --grid-shp "$GRID_SHP" \
            --city-gpkg "$CITY_GPKG" \
            --city-layer "$CITY_LAYER" \
            --output "$LOOKUP_FILE" \
            --save-overlaps "$OVERLAP_FILE"
    fi
else
    python create_city_smoke_lookup.py \
        --smoke-csv "$SMOKE_CSV" \
        --grid-shp "$GRID_SHP" \
        --city-gpkg "$CITY_GPKG" \
        --city-layer "$CITY_LAYER" \
        --output "$LOOKUP_FILE" \
        --save-overlaps "$OVERLAP_FILE"
fi

if [ ! -f "$LOOKUP_FILE" ]; then
    echo "❌ Error: Lookup file was not created"
    exit 1
fi

# =============================================================================
# STEP 2: VERIFY VISIT DATA
# =============================================================================

echo "=============================================================================="
echo "STEP 2: Verifying Visit Data"
echo "=============================================================================="
echo ""

if [ ! -d "$VISIT_DIR" ]; then
    echo "❌ Error: Visit directory not found: $VISIT_DIR"
    echo ""
    echo "Please run your visit processing first:"
    echo "  python process_multilayer_visits.py \\"
    echo "    --tsv-dir <your_tsv_dir> \\"
    echo "    --gpkg-dir <your_gpkg_dir> \\"
    echo "    --output-dir $VISIT_DIR"
    exit 1
fi

# Count visit files
VISIT_COUNT=$(ls -1 ${VISIT_DIR}/*_daily.csv 2>/dev/null | wc -l)

if [ $VISIT_COUNT -eq 0 ]; then
    echo "❌ Error: No visit files found in $VISIT_DIR"
    echo ""
    echo "Expected files like: 10056774_Arlington_parks_daily.csv"
    exit 1
fi

echo "✅ Found $VISIT_COUNT visit files"
echo ""
echo "Sample files:"
ls -1 ${VISIT_DIR}/*_daily.csv | head -5
echo ""

# =============================================================================
# STEP 3: ANALYZE VISITS WITH SMOKE
# =============================================================================

echo "=============================================================================="
echo "STEP 3: Analyzing Visits with Smoke Data"
echo "=============================================================================="
echo ""

SUMMARY_FILE="${OUTPUT_DIR}/city_smoke_visit_summary.csv"
DETAILED_FILE="${OUTPUT_DIR}/visits_with_smoke_detailed.csv"

python analyze_visits_smoke.py \
    --visit-dir "$VISIT_DIR" \
    --smoke-lookup "$LOOKUP_FILE" \
    --output-summary "$SUMMARY_FILE" \
    --output-detailed "$DETAILED_FILE" \
    $SIMPLE_CATEGORIES

if [ ! -f "$SUMMARY_FILE" ]; then
    echo "❌ Error: Summary file was not created"
    exit 1
fi

echo ""
echo "✅ Step 3 Complete"
echo "   Summary: $SUMMARY_FILE"
echo "   Detailed: $DETAILED_FILE"
echo ""

# =============================================================================
# FINAL SUMMARY
# =============================================================================

echo "=============================================================================="
echo "✅ WORKFLOW COMPLETE!"
echo "=============================================================================="
echo ""
echo "Output files:"
echo "  1. City-Smoke Lookup:  $LOOKUP_FILE"
echo "  2. Grid Overlaps:      $OVERLAP_FILE"
echo "  3. City Summary:       $SUMMARY_FILE"
echo "  4. Detailed Data:      $DETAILED_FILE"
echo ""
echo "=============================================================================="

# Generate quick preview
echo ""
echo "Preview of results (top 10 cities by total visits):"
echo "----------------------------------------------------------------------"
head -11 "$SUMMARY_FILE" | column -t -s,
echo "----------------------------------------------------------------------"
echo ""

# Generate basic statistics
echo "Basic statistics:"
echo "----------------------------------------------------------------------"
python3 << EOF
import pandas as pd
import numpy as np

df = pd.read_csv('$SUMMARY_FILE')

print(f"Total cities analyzed: {len(df)}")
print(f"Total visits across all cities: {df['total_visits'].sum():,}")
print(f"Average visits per city: {df['total_visits'].mean():,.0f}")
print(f"\nSmoke levels:")
print(f"  Mean smoke PM: {df['avg_smoke_pm'].mean():.2f} µg/m³")
print(f"  Min smoke PM: {df['min_smoke_pm'].min():.2f} µg/m³")
print(f"  Max smoke PM: {df['max_smoke_pm'].max():.2f} µg/m³")
EOF

echo "----------------------------------------------------------------------"
echo ""
echo "Analysis complete! Check $OUTPUT_DIR for all output files."
