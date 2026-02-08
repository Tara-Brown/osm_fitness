#!/bin/bash
# Setup script to submit parallel processing jobs

cd /mnt/beegfs/projects/cellphone_mobility

# Create logs directory
mkdir -p logs

# Get list of unique TSV files (handling multi-part files)
TSV_DIR="test_cities"
NUM_FILES=$(ls ${TSV_DIR}/*.tsv | sed 's/_pin_report_[0-9]\{3\}\.tsv/_pin_report.tsv/' | sed 's/_pin_report\.tsv/_pin_report.tsv/' | sort -u | wc -l)

echo "Found $NUM_FILES unique TSV file groups to process"

# Update the sbatch script with correct array size
sed "s/#SBATCH --array=0-N/#SBATCH --array=0-$((NUM_FILES-1))/" submit_parallel_jobs.sh > submit_parallel_jobs_temp.sh

# Make it executable
chmod +x submit_parallel_jobs_temp.sh

# Submit the job array
echo "Submitting job array..."
sbatch submit_parallel_jobs_temp.sh

echo "Jobs submitted! Check status with: squeue -u $USER"
echo "Monitor logs in: logs/"
