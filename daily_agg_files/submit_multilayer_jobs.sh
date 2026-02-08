#!/bin/bash
#SBATCH --job-name=cellphone_multilayer
#SBATCH --output=logs/cellphone_%A_%a.out
#SBATCH --error=logs/cellphone_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --array=0-N  # Replace N with (number of files - 1)

# =============================================================================
# SLURM Array Job Script for Processing Cellphone Mobility Data
# =============================================================================
# This script processes multiple TSV files in parallel (for city groups) using SLURM job arrays.
# Each array task processes one unique TSV file group.
#
# Usage:
#   1. Update the paths below
#   2. Update --array=0-N where N = number_of_unique_tsv_files - 1
#   3. Run: sbatch submit_multilayer_jobs.sh
# =============================================================================

# =============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# =============================================================================

# Base directory
BASE_DIR="/mnt/beegfs/projects/cellphone_mobility"
TSV_DIR="${BASE_DIR}/testing_files/unzipped_Washington"
GPKG_DIR="${BASE_DIR}/testing_files/geopackage_layers/state_53_components.gpkg"
OUTPUT_DIR="${BASE_DIR}/multilayer_outputs"
SCRIPT_DIR="${BASE_DIR}/testing_files/daily_agg_and_smoke"
PROCESS_SCRIPT="${SCRIPT_DIR}/process_single_tsv_multilayer.py"

# Processing options
BY_POLYGON=""  # Set to "--by-polygon" to enable polygon-level aggregation

# =============================================================================
# SETUP
# =============================================================================

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p logs

# Get list of unique TSV files (handling multi-part files)
# This removes the _XXX suffix to get unique base names
TSV_FILES=($(ls ${TSV_DIR}/*_pin_report*.tsv | \
             sed 's/_pin_report_[0-9]\{3\}\.tsv/_pin_report.tsv/' | \
             sort -u))

# Get the TSV file for this array task
TSV_FILE="${TSV_FILES[$SLURM_ARRAY_TASK_ID]}"
# If base file doesn't exist, try to find the first part
if [ ! -f "$TSV_FILE" ]; then
    BASE_NAME=$(basename "$TSV_FILE" .tsv)
    FIRST_PART="${TSV_DIR}/${BASE_NAME}_000.tsv"
    
    if [ -f "$FIRST_PART" ]; then
        echo "Using first part: $FIRST_PART"
        TSV_FILE="$FIRST_PART"
    else
        echo "ERROR: Neither $TSV_FILE nor $FIRST_PART found"
        exit 1
    fi
fi

# =============================================================================
# RUN PROCESSING
# =============================================================================

echo "Starting processing..."
echo ""

python3 "$PROCESS_SCRIPT" "$TSV_FILE" "$GPKG_DIR" "$OUTPUT_DIR" $BY_POLYGON

EXIT_CODE=$?

echo ""
echo "=============================================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Processing completed successfully"
else
    echo "❌ Processing failed with exit code: $EXIT_CODE"
fi
echo "End Time: $(date)"
echo "=============================================================================="

exit $EXIT_CODE
