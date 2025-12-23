#!/bin/bash
#SBATCH --job-name=azira_dl
#SBATCH --output=/mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira/logs/%x-%j.out
#SBATCH --error=/mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira/logs/%x-%j.err
#SBATCH --time=01:00:00

set -euo pipefail

BASE="/mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira"
mkdir -p "$BASE/logs"
cd "$BASE"

# (Optional) activate conda / modules
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate your_env

# Run the downloader
bash run_download_urls.sh >> "$BASE/logs/download_$(date +%F).log" 2>&1

# Schedule the next run at 00:15 tomorrow (scheduler time)
#NEXT="$(date -d 'tomorrow 00:15' '+%Y-%m-%dT%H:%M:%S')"
#sbatch --begin="$NEXT" "$BASE/azira_daily.sbatch" >/dev/null