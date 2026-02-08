#!/bin/bash
# =============================================================================
# Setup and Submit Script for Multilayer Processing
# =============================================================================
# This script counts the TSV files, updates the SLURM script, and submits it.
#
# Usage: bash setup_and_submit_multilayer.sh
# =============================================================================

# Configuration
# Change these depending on where your files are
BASE_DIR="/mnt/beegfs/projects/cellphone_mobility/testing_files"
# Only Washington, for now
TSV_DIR="${BASE_DIR}/unzipped_Washington"
SCRIPT_NAME="daily_agg_files/submit_multilayer_jobs.sh"

echo "=============================================================================="
echo "Cellphone Mobility Multilayer Processing Setup"
echo "=============================================================================="
echo ""

# Navigate to base directory
cd "$BASE_DIR" || exit 1

# Create logs directory
mkdir -p logs

# Count unique TSV file groups
echo "Counting TSV files in: $TSV_DIR"
NUM_FILES=$(ls ${TSV_DIR}/*_pin_report*.tsv 2>/dev/null | \
            sed 's/_pin_report_[0-9]\{3\}\.tsv/_pin_report.tsv/' | \
            sort -u | wc -l)

if [ $NUM_FILES -eq 0 ]; then
    echo "❌ Error: No TSV files found in $TSV_DIR"
    exit 1
fi

echo "✅ Found $NUM_FILES unique TSV file groups to process"
echo ""

# Calculate array index (0-based)
MAX_INDEX=$((NUM_FILES - 1))

# Update the SLURM script with correct array size
echo "Updating $SCRIPT_NAME with array size: 0-$MAX_INDEX"
sed "s/#SBATCH --array=0-N/#SBATCH --array=0-$MAX_INDEX/" "$SCRIPT_NAME" > "${SCRIPT_NAME}.tmp"

# Make it executable
chmod +x "${SCRIPT_NAME}.tmp"

# Show sample of files that will be processed
echo ""
echo "Sample of files to be processed:"
echo "----------------------------------------------------------------------"
ls ${TSV_DIR}/*_pin_report*.tsv 2>/dev/null | \
   sed 's/_pin_report_[0-9]\{3\}\.tsv/_pin_report.tsv/' | \
   sort -u | head -10
echo ""
if [ $NUM_FILES -gt 10 ]; then
    echo "... and $((NUM_FILES - 10)) more files"
    echo ""
fi

# Ask for confirmation
echo "=============================================================================="
read -p "Submit $NUM_FILES processing jobs? (y/n): " -n 1 -r
echo ""
echo "=============================================================================="

if [[ $REPLY =~ ^[Yy]$ ]]; then
    # Submit the job array
    echo "Submitting job array..."
    JOB_ID=$(sbatch "${SCRIPT_NAME}.tmp" | awk '{print $NF}')
    
    if [ $? -eq 0 ]; then
        echo "Check logs in: logs/"
        echo "   logs/cellphone_${JOB_ID}_*.out"
        echo "   logs/cellphone_${JOB_ID}_*.err"
    else
        echo "❌ Error: Job submission failed"
        exit 1
    fi
else
    echo "Submission cancelled"
    rm -f "${SCRIPT_NAME}.tmp"
    exit 0
fi

# Clean up temp file
rm -f "${SCRIPT_NAME}.tmp"

echo ""
echo "=============================================================================="
echo "Setup complete!"
echo "=============================================================================="
