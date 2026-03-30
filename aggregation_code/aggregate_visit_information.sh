#!/bin/bash
# =============================================================================
# Setup and Submit Script for Visit Statistics
# =============================================================================
# This script counts the TSV files, updates the SLURM script, and submits it.
#
# Usage: bash aggregate_visit_information
# =============================================================================

# Configuration
# Change these depending on where your files are
BASE_DIR="/mnt/beegfs/projects/cellphone_mobility"
# Only Oregon, for now
TSV_DIR="${BASE_DIR}/testing_files/unzipped_files"
SCRIPT_NAME="aggregation_code/submit_visit_jobs.sh"

# Job submission limit
MAX_CONCURRENT_JOBS=80

echo "=============================================================================="
echo "Understanding Cellphone Ping Data Better"
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

# Update the SLURM script with correct array size AND concurrency limit
echo "Updating $SCRIPT_NAME with array size: 0-$MAX_INDEX%$MAX_CONCURRENT_JOBS"
sed "s/#SBATCH --array=0-N/#SBATCH --array=0-$MAX_INDEX%$MAX_CONCURRENT_JOBS/" "$SCRIPT_NAME" > "${SCRIPT_NAME}.tmp"

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
echo "Job Configuration:"
echo "  Total jobs: $NUM_FILES"
echo "  Max concurrent: $MAX_CONCURRENT_JOBS"
echo "  Remaining jobs will queue automatically"
echo "=============================================================================="
read -p "Submit $NUM_FILES processing jobs (max $MAX_CONCURRENT_JOBS at a time)? (y/n): " -n 1 -r
echo ""
echo "=============================================================================="

if [[ $REPLY =~ ^[Yy]$ ]]; then
    # Submit the job array
    echo "Submitting job array..."
    JOB_ID=$(sbatch "${SCRIPT_NAME}.tmp" | awk '{print $NF}')
    
    if [ $? -eq 0 ]; then
        echo "✅ Job array submitted successfully: $JOB_ID"
        echo ""
        echo "Check status with: squeue -j $JOB_ID"
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
