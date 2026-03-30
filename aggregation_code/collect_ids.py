'''
Collects Ubermedia IDs from TSV files whose names don't match any known state
pattern (e.g. truncated names like *New_expanded... instead of *New_Mexico_expanded...)
and saves them to shared/unrecognized_ids.parquet for use by per-state jobs.

Run this once before submitting per-state jobs if precompute_shared.py has
already been run without this step.

Usage:
    python collect_unrecognized_ids.py
'''
import glob
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import geopandas as gpd

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE_DIR = "/mnt/beegfs/projects/cellphone_mobility"
TSV_DIR  = f"{BASE_DIR}/testing_files/unzipped_files"
SHARED   = Path(f"{BASE_DIR}/shared")

STATE_NAMES = {
    'AL': 'Alabama', 'AZ': 'Arizona', 'AR': 'Arkansas', 'CA': 'California',
    'CO': 'Colorado', 'CT': 'Connecticut', 'DE': 'Delaware', 'DC': 'Washington_DC',
    'FL': 'Florida', 'GA': 'Georgia', 'ID': '_ID', 'IL': 'Illinois',
    'IN': 'Indiana', 'IA': 'Iowa', 'KS': 'Kansas', 'KY': 'Kentucky',
    'LA': 'Louisiana', 'ME': 'Maine', 'MD': 'Maryland', 'MA': 'Massachusetts',
    'MI': 'Michigan', 'MN': 'Minnesota', 'MS': 'Mississippi', 'MO': 'Missouri',
    'MT': 'Montana', 'NE': 'Nebraska', 'NV': 'Nevada', 'NH': 'New_Hampshire',
    'NJ': 'New_Jersey', 'NM': 'New_Mexico', 'NY': 'New_York', 'NC': 'North_Carolina',
    'ND': '_ND', 'OH': 'Ohio', 'OK': 'Oklahoma', 'OR': '_OR',
    'PA': 'Pennsylvania', 'RI': 'Rhode_Island', 'SC': 'South_Carolina',
    'SD': 'South_Dakota', 'TN': 'Tennesee', 'TX': 'Texas', 'UT': 'Utah', #Tennessee was spelled wrong earlier-propogating the cycle
    'VT': 'Vermont', 'VA': 'Virginia', 'WA': 'Washington', 'WV': 'West_Virginia',
    'WI': 'Wisconsin', 'WY': 'Wyoming', 
}

suffix = "_expanded_cel_cdl_report_cel.tsv"

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

# ---------------------------------------------------------------------------
# Identify unrecognized files
# ---------------------------------------------------------------------------
all_tsv_files = set(glob.glob(f"{TSV_DIR}/*{suffix}"))

known_tsv_files = set()
for name in STATE_NAMES.values():
    known_tsv_files.update(glob.glob(f"{TSV_DIR}/*{name}{suffix}"))

unrecognized_files = sorted(all_tsv_files - known_tsv_files)
print(f"Found {len(unrecognized_files)} unrecognized TSV file(s):")
for f in unrecognized_files:
    print(f"  {Path(f).name}")

if not unrecognized_files:
    print("Nothing to do.")
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# Load precomputed shared files
# ---------------------------------------------------------------------------
acs_path = SHARED / "acs_all_states.parquet"
bg_path  = SHARED / "bg_all_states.gpkg"

for p in [acs_path, bg_path]:
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. Run precompute_shared.py before this script."
        )

print("Loading precomputed ACS data...")
acs_df = pd.read_parquet(acs_path)
acs_df["GEOID"] = acs_df["GEOID"].astype(str)

print("Loading precomputed block group geometries...")
bg = gpd.read_file(bg_path)
print(f"  {len(bg):,} block groups loaded")

# ---------------------------------------------------------------------------
# Full spatial join pipeline — same as process_cel.py
# ---------------------------------------------------------------------------
ubermedia_lookup = {}

for file in tqdm(unrecognized_files, desc="Processing unrecognized files"):
    for chunk in pd.read_csv(file, sep="\t", chunksize=500_000, dtype=dtype_dict):

        chunk = chunk.dropna(subset=["Common Evening Lat", "Common Evening Long"])
        if chunk.empty:
            continue

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

        joined_deduped = joined.drop_duplicates(subset="Hashed Ubermedia Id", keep="first")
        new_records = (
            joined_deduped
            .set_index("Hashed Ubermedia Id")[["GEOID", "median_income", "avg_income", "pct_bachelors"]]
            .rename(columns={"avg_income": "average_income", "pct_bachelors": "pct_education"})
        )
        ubermedia_lookup.update(new_records.to_dict(orient="index"))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
lookup_df = pd.DataFrame.from_dict(ubermedia_lookup, orient="index").reset_index()
lookup_df.rename(columns={"index": "Hashed Ubermedia Id"}, inplace=True)

out_path = SHARED / "unrecognized_lookup.parquet"
lookup_df.to_parquet(out_path, index=False)
print(f"Saved {len(lookup_df):,} devices → {out_path}")