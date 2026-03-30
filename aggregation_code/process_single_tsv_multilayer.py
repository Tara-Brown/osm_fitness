#/usr/bin/env python3
"""Wrapper for processing single TSV file."""
import sys
from pathlib import Path

if len(sys.argv) < 4:
    print("Usage: python process_single_tsv_multilayer.py <tsv_file> <gpkg_file> <output_dir> [--by-polygon]")
    sys.exit(1)

tsv_file = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
by_polygon = '--by-polygon' in sys.argv

if not tsv_file.exists():
    print(f"Error: TSV file not found: {tsv_file}")
    sys.exit(1)

print(f"Processing: {tsv_file.name}")
print(f"Output directory: {output_dir}")
print(f"Polygon-level aggregation: {by_polygon}")

# Import and run
from process_multilayer_visits import process_tsv_file

try:
    stats = process_tsv_file(tsv_file, output_dir, by_polygon=by_polygon)
    
    if stats.get('saved_files'):
        print(f"\n✅ Success! Processed {stats['total_points']:,} points")
        sys.exit(0)
    else:
        print(f"\n⚠️  No output files created")
        sys.exit(1)
        
except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
