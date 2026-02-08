# Cellphone Mobility Multilayer Processing

This toolkit processes cellphone mobility data and matches visits to multiple GeoPackage layers (parks, gyms, sports centres) with daily aggregation of unique device visits.

## Overview

For each city, the output shows the number of unique device visits to each layer per day.

## Files

- `process_multilayer_visits.py` - Main processing script with full CLI
- `process_single_tsv_multilayer.py` - Wrapper for processing single files (for parallel execution)
- `submit_multilayer_jobs.sh` - SLURM batch script for HPC parallel processing
- `setup_and_submit_multilayer.sh` - Helper script to count files and submit jobs
- `combine_results.py` - Utility to combine individual city outputs into summary files


usage: process_multilayer_visits.py [-h] (--tsv-dir TSV_DIR | --tsv-file TSV_FILE)
                                     --gpkg-dir GPKG_DIR [--gpkg-file GPKG_FILE]
                                     --output-dir OUTPUT_DIR [--by-polygon]
                                     [--layers LAYERS [LAYERS ...]]
                                     [--chunk-size CHUNK_SIZE]

Required arguments:
  --tsv-dir TSV_DIR         Directory containing TSV files
  --tsv-file TSV_FILE       Single TSV file to process
  --gpkg-dir GPKG_DIR       Directory containing GeoPackage files
  --output-dir OUTPUT_DIR   Output directory for CSV files

Optional arguments:
  --gpkg-file GPKG_FILE     Specific GeoPackage file (auto-detected if not specified)
  --by-polygon              Aggregate by individual polygons within layers
  --layers LAYER [LAYER...] Specify which layers to process
  --chunk-size SIZE         Rows to process at once (default: 100000)
```

## Output Files

### Layer-Level Aggregation (default)

For each city and layer, creates a CSV file:
```
{city_code}_{city_name}_{layer_name}_daily.csv
```

Columns:
- `date`: Date of visit
- `layer`: Layer name (parkserve, osm_parks, etc.)
- `unique_visits`: Number of unique devices that visited this layer on this date
- `city`: City name
- `city_code`: City numeric code

Example: `10056774_Arlington_parkserve_daily.csv`

### Polygon-Level Aggregation (--by-polygon flag)

For each city and layer, creates a CSV file:
```
{city_code}_{city_name}_{layer_name}_by_polygon_daily.csv
```

Columns:
- `date`: Date of visit
- `layer`: Layer name
- `polygon_id`: Individual polygon/feature ID within the layer
- `unique_visits`: Number of unique devices that visited this polygon on this date
- `city`: City name
- `city_code`: City numeric code

Example: `10056774_Arlington_parkserve_by_polygon_daily.csv`

### Summary Report

A summary CSV is created showing processing statistics:
```
processing_summary.csv
```

Columns:
- `city_code`: City numeric code
- `city_name`: City name
- `total_points`: Total cellphone pings processed
- `total_chunks`: Number of chunks processed
- `{layer}_points`: Number of points matched to each layer


