#!/usr/bin/env python3
"""Wrapper for processing single TSV file — raw daily device visit counts, no spatial join."""
import sys
import pandas as pd
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage: python process_tsvs_visit.py <tsv_file> <output_dir>")
    sys.exit(1)

tsv_file = Path(sys.argv[1])
output_dir = Path(sys.argv[2])

if not tsv_file.exists():
    print(f"Error: TSV file not found: {tsv_file}")
    sys.exit(1)

print(f"Processing: {tsv_file.name}")
print(f"Output directory: {output_dir}")

output_dir.mkdir(parents=True, exist_ok=True)

from process_multilayer_visits import find_tsv_parts, chunk_tsv_by_day, TIME_COL, CHUNK_SIZE, DEVICE_COL


def aggregate_daily_visits(gdf_chunk):
    """Count how many times each device was spotted per day."""
    if len(gdf_chunk) == 0:
        return pd.DataFrame(), pd.DataFrame()

    if 'date' not in gdf_chunk.columns:
        if 'datetime' in gdf_chunk.columns:
            gdf_chunk['date'] = gdf_chunk['datetime'].dt.date
        elif TIME_COL in gdf_chunk.columns:
            gdf_chunk['datetime'] = pd.to_datetime(gdf_chunk[TIME_COL], unit='s', errors='coerce')
            gdf_chunk['date'] = gdf_chunk['datetime'].dt.date
        else:
            return pd.DataFrame(), pd.DataFrame()

    daily = gdf_chunk.groupby(['date', DEVICE_COL]).size().reset_index(name='visit_count')

    ping_spans = (
        gdf_chunk.groupby(['date', DEVICE_COL])['datetime']
        .agg(first_ping='min', last_ping='max')
        .reset_index()
    )

    return daily, ping_spans

try:
    tsv_parts = find_tsv_parts(tsv_file)
    print(f"  Found {len(tsv_parts)} TSV part(s)")

    stats = {
        'total_chunks': 0,
        'total_points': 0,
        'saved_files': []
    }

    all_daily_device_counts = []
    all_daily_device_pings = []
    for part_num, tsv_part in enumerate(tsv_parts, 1):
        if not tsv_part.exists():
            print(f"  ⚠️  Part {part_num} does not exist: {tsv_part.name}")
            continue

        print(f"  📄 Processing part {part_num}/{len(tsv_parts)}: {tsv_part.name}")
        chunk_count = 0

        for gdf_chunk in chunk_tsv_by_day(tsv_part, target_rows=CHUNK_SIZE):
            chunk_count += 1
            stats['total_chunks'] += 1
            stats['total_points'] += len(gdf_chunk)

            daily_counts, ping_intervals = aggregate_daily_visits(gdf_chunk)
            if len(daily_counts) > 0:
                all_daily_device_counts.append(daily_counts)
                all_daily_device_pings.append(ping_intervals)

        print(f"    Processed {chunk_count} chunks from this part")

    if all_daily_device_counts:
        print(f"\n  Combining {len(all_daily_device_counts)} chunk results...")
        combined = pd.concat(all_daily_device_counts, ignore_index=True)
        
        print(f"  Aggregating {len(combined):,} records...")
        final_counts = (
            combined
            .groupby(['date', DEVICE_COL])['visit_count']
            .sum()
            .reset_index()
        )

        combined_pings = pd.concat(all_daily_device_pings, ignore_index=True)
        final_pings = (
            combined_pings
            .groupby(['date', DEVICE_COL])
            .agg(first_ping=('first_ping', 'min'), last_ping=('last_ping', 'max'))
            .reset_index()
        )      
        
        # Merge counts and ping spans
        final_counts = final_counts.merge(final_pings, on=['date', DEVICE_COL], how='left')

        # Ensure date is datetime for parquet compatibility
        final_counts['date'] = pd.to_datetime(final_counts['date'])

        output_file = output_dir / f"daily_device_instances_{tsv_file.stem}.parquet"
        final_counts.to_parquet(output_file, index=False, engine='pyarrow')

        print(f"\n📊 Daily device instance counts saved to: {output_file}")
        print(f"   Total unique (date, device) pairs: {len(final_counts):,}")
        stats['saved_files'].append(str(output_file))

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
