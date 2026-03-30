'''
Per-state CEL processing job.
Reads CEL files for a given state, performs spatial join to block groups,
joins ACS statistics, and saves a per-device lookup parquet.

Appends fully resolved records from unrecognized/truncated TSV filenames
(precomputed once by collect_unrecognized_ids.py) so no device is dropped.

Requires the following to have been run first:
  - precompute_shared.py       → shared/acs_all_states.parquet, shared/bg_all_states.gpkg
  - collect_unrecognized_ids.py → shared/unrecognized_lookup.parquet

Usage:
    python process_cel.py --state_abbr LA --state_name Louisiana
    python process_cel.py --state_abbr DC --state_name District_of_Columbia
'''
import argparse
import glob
import geopandas as gpd
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--state_abbr", required=True, help="Two-letter state abbreviation, e.g. LA")
parser.add_argument("--state_name", required=True, help="Full state name used in filenames, e.g. Louisiana or New_Hampshire")
args = parser.parse_args()

STATE_ABBR = args.state_abbr.upper()
STATE_NAME = args.state_name  # used verbatim in glob and parquet filename

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE_DIR = "/mnt/beegfs/projects/cellphone_mobility"
TSV_DIR  = f"{BASE_DIR}/testing_files/unzipped_files"
SHARED   = Path(f"{BASE_DIR}/shared")

STATE_FIPS = {
    'AL': '01', 'AZ': '04', 'AR': '05', 'CA': '06', 'CO': '08',
    'CT': '09', 'DE': '10', 'DC': '11', 'FL': '12', 'GA': '13',
    'ID': '16', 'IL': '17', 'IN': '18', 'IA': '19', 'KS': '20',
    'KY': '21', 'LA': '22', 'ME': '23', 'MD': '24', 'MA': '25',
    'MI': '26', 'MN': '27', 'MS': '28', 'MO': '29', 'MT': '30',
    'NE': '31', 'NV': '32', 'NH': '33', 'NJ': '34', 'NM': '35',
    'NY': '36', 'NC': '37', 'ND': '38', 'OH': '39', 'OK': '40',
    'OR': '41', 'PA': '42', 'RI': '44', 'SC': '45', 'SD': '46',
    'TN': '47', 'TX': '48', 'UT': '49', 'VT': '50', 'VA': '51',
    'WA': '53', 'WV': '54', 'WI': '55', 'WY': '56',
}

if STATE_ABBR not in STATE_FIPS:
    raise ValueError(f"Unknown or unsupported state abbreviation: {STATE_ABBR}")

# ---------------------------------------------------------------------------
# Load precomputed shared files
# ---------------------------------------------------------------------------
acs_path          = SHARED / "acs_all_states.parquet"
bg_path           = SHARED / "bg_all_states.gpkg"
unrecognized_path = SHARED / "unrecognized_lookup.parquet"

for p in [acs_path, bg_path, unrecognized_path]:
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. Run precompute_shared.py and collect_unrecognized_ids.py before submitting jobs."
        )

print("Loading precomputed ACS data...")
acs_df = pd.read_parquet(acs_path)
acs_df["GEOID"] = acs_df["GEOID"].astype(str)
print(f"  {len(acs_df):,} block groups loaded")

print("Loading precomputed block group geometries...")
bg = gpd.read_file(bg_path)
print(f"  {len(bg):,} block groups loaded")

print("Loading unrecognized file lookup...")
unrecognized_lookup = pd.read_parquet(unrecognized_path)
print(f"  {len(unrecognized_lookup):,} unrecognized devices loaded")

# ---------------------------------------------------------------------------
# Main chunk processing for this state's TSV files
# ---------------------------------------------------------------------------
dtype_dict = {
    "Hashed Ubermedia Id":    str,
    "Common Evening Lat":     float,
    "Common Evening Long":    float,
    "Polygon Id":             str,
    "First Visit Timestamp":  str,
    "Common Evening Country": str,
    "Common Evening Census":  str,
    "Common Evening Postal1": str,
    "Common Evening Postal2": str,
}

ubermedia_lookup = {}

suffix      = "_expanded_cel_cdl_report_cel.tsv"
tsv_pattern = f"{TSV_DIR}/*{STATE_NAME}{suffix}"
tsv_files   = glob.glob(tsv_pattern)
print(f"Found {len(tsv_files)} TSV file(s) matching: {tsv_pattern}")

for file in tqdm(tsv_files, desc=f"Processing {STATE_ABBR}"):
    for chunk in pd.read_csv(file, sep="\t", chunksize=500_000, dtype=dtype_dict):

        chunk = chunk.dropna(subset=["Common Evening Lat", "Common Evening Long"])
        if chunk.empty:
            continue

        # Skip devices already seen — avoids spatial join cost for known devices
        new_only = chunk[~chunk["Hashed Ubermedia Id"].isin(ubermedia_lookup)]
        if new_only.empty:
            continue

        gdf = gpd.GeoDataFrame(
            new_only,
            geometry=gpd.points_from_xy(new_only["Common Evening Long"], new_only["Common Evening Lat"]),
            crs="EPSG:4326"
        )

        joined = gpd.sjoin(gdf, bg, how="left", predicate="within")
        joined = joined.merge(acs_df, on="GEOID", how="left")

        # Deduplicate within chunk then update lookup — no iterrows
        joined_deduped = joined.drop_duplicates(subset="Hashed Ubermedia Id", keep="first")
        new_records = (
            joined_deduped
            .set_index("Hashed Ubermedia Id")[["GEOID", "median_income", "avg_income", "pct_bachelors"]]
            .rename(columns={"avg_income": "average_income", "pct_bachelors": "pct_education"})
        )
        ubermedia_lookup.update(new_records.to_dict(orient="index"))

# ---------------------------------------------------------------------------
# Append records from unrecognized files, skipping any already resolved above
# ---------------------------------------------------------------------------
existing_ids = set(ubermedia_lookup.keys())
new_unrecognized = unrecognized_lookup[
    ~unrecognized_lookup["Hashed Ubermedia Id"].isin(existing_ids)
]
new_records = new_unrecognized.set_index("Hashed Ubermedia Id").to_dict(orient="index")
ubermedia_lookup.update(new_records)
print(f"Appended {len(new_records):,} devices from unrecognized files")

# ---------------------------------------------------------------------------
# Save output
# ---------------------------------------------------------------------------
Path("lookups").mkdir(exist_ok=True)

lookup_df = pd.DataFrame.from_dict(ubermedia_lookup, orient="index").reset_index()
lookup_df.rename(columns={"index": "Hashed Ubermedia Id"}, inplace=True)

# e.g. lookups/ubermedia_lookup_Louisiana.parquet
#      lookups/ubermedia_lookup_New_Hampshire.parquet
out_path = f"lookups/ubermedia_lookup_{STATE_NAME}.parquet"
lookup_df.to_parquet(out_path, index=False)

print(f"Saved {len(lookup_df):,} devices → {out_path}")
