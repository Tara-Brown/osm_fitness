#!/bin/bash

# Load conda and activate environment
module load conda

# init conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activate env
conda activate ox

# Run script
python /mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira/download_urls.py