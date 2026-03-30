#!/bin/bash
# test_single_file.sh - Test processing on one TSV file

# Configuration
BASE_DIR="/mnt/beegfs/projects/cellphone_mobility"
TSV_DIR="${BASE_DIR}/testing_files/unzipped_Oregon"
GPKG_FILE="/mnt/beegfs/hellgate/home/tb208541/osm_fitness/unified_geopackages_polygon_level/state_41_components.gpkg"
OUTPUT_DIR="${BASE_DIR}/multilayer_outputs/test_run"
SCRIPT_DIR="${BASE_DIR}/testing_files/daily_agg_and_smoke"
PROCESS_SCRIPT="${SCRIPT_DIR}/process_single_tsv_multilayer.py"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Find the first TSV file
TSV_FILE=$(ls ${TSV_DIR}/*_pin_report*.tsv | head -1)

if [ -z "$TSV_FILE" ]; then
    echo "❌ No TSV files found in $TSV_DIR"
    exit 1
fi

echo "=============================================================================="
echo "Testing with single file"
echo "=============================================================================="
echo "TSV File: $TSV_FILE"
echo "GPKG File: $GPKG_FILE"
echo "Output Dir: $OUTPUT_DIR"
echo ""
echo "Starting processing..."
echo "=============================================================================="

# Run the processing
python3 "$PROCESS_SCRIPT" "$TSV_FILE" "$GPKG_FILE" "$OUTPUT_DIR"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=============================================================================="
    echo "✅ Test completed successfully!"
    echo "=============================================================================="
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}"/*.csv 2>/dev/null || echo "No CSV files generated"
else
    echo ""
    echo "=============================================================================="
    echo "❌ Test failed with exit code: $EXIT_CODE"
    echo "=============================================================================="
fi

exit $EXIT_CODE
