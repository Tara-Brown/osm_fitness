'''
One-time setup script. Run this ONCE before submitting per-state jobs.

Produces two shared files that every per-state job reads at startup:
  - shared/acs_all_states.parquet   : ACS block-group statistics for all states
  - shared/bg_all_states.gpkg       : merged block-group geometries with spatial index

Usage:
    python precompute_shared.py
'''
import glob
import time
import numpy as np
import geopandas as gpd
import pandas as pd
import requests
from tqdm import tqdm
from pathlib import Path

# ---------------------------------------------------------------------------
# Config — keep in sync with process_cel.py
# ---------------------------------------------------------------------------
YEAR     = 2022
BASE_URL = f"https://api.census.gov/data/{YEAR}/acs/acs5"

VARS = [
    "B19013_001E",  # median income
    "B19001_001E",  # households
    "B19025_001E",  # aggregate income
    "B15003_001E",  # pop 25+
    "B15003_022E", "B15003_023E",
    "B15003_024E", "B15003_025E"
]

BASE_DIR = "/mnt/beegfs/projects/cellphone_mobility"
BG_DIR   = f"{BASE_DIR}/testing_files/shapefiles"
OUT_DIR  = Path(f"{BASE_DIR}/shared")

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

OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. ACS — fetch all states, save as parquet
# ---------------------------------------------------------------------------
def fetch_state_blockgroups(state_fips):
    params = {
        "get": ",".join(VARS),
        "for": "block group:*",
        "in": f"state:{state_fips} county:* tract:*",
        "key": "0ee12a5952120a2347aa982ff6209cec33ee5053",
    }
    r = requests.get(BASE_URL, params=params, timeout=60)
    r.raise_for_status()

    data = r.json()
    df   = pd.DataFrame(data[1:], columns=data[0])
    df["GEOID"] = df["state"] + df["county"] + df["tract"] + df["block group"]

    for v in VARS:
        df[v] = pd.to_numeric(df[v], errors="coerce").replace(-666666666, np.nan)

    df["avg_income"]    = df["B19025_001E"] / df["B19001_001E"].replace(0, np.nan)
    bachelors           = df[["B15003_022E","B15003_023E","B15003_024E","B15003_025E"]].sum(axis=1, skipna=False)
    df["pct_bachelors"] = bachelors / df["B15003_001E"].replace(0, np.nan)

    return df[["GEOID","B19013_001E","avg_income","pct_bachelors"]].rename(
        columns={"B19013_001E": "median_income"}
    )

acs_path = OUT_DIR / "acs_all_states.parquet"
print("--- Step 1: Fetching ACS data for all states ---")
dfs = []
for abbr, fips in tqdm(STATE_FIPS.items(), desc="ACS fetch"):
    try:
        df = fetch_state_blockgroups(fips)
        dfs.append(df)
        time.sleep(0.5)  # polite pause
    except Exception as e:
        print(f"  Warning: ACS fetch failed for {abbr}: {e}")

acs_df = pd.concat(dfs, ignore_index=True)
acs_df["GEOID"] = acs_df["GEOID"].astype(str)
acs_df.to_parquet(acs_path, index=False)
print(f"Saved {len(acs_df):,} block groups → {acs_path}\n")

# ---------------------------------------------------------------------------
# 2. Shapefiles — merge all states, save as GeoPackage (includes spatial index)
# ---------------------------------------------------------------------------
bg_path = OUT_DIR / "bg_all_states.gpkg"
print("--- Step 2: Merging block group shapefiles ---")
gdfs = []
for abbr, fips in tqdm(STATE_FIPS.items(), desc="Shapefiles"):
    shp_files = glob.glob(f"{BG_DIR}/tl_202*_{fips}_bg.shp")
    if not shp_files:
        print(f"  Warning: no shapefile for {abbr} (FIPS {fips}), skipping")
        continue
    gdf = gpd.read_file(shp_files[0])[["GEOID","geometry"]].to_crs(epsg=4326)
    gdfs.append(gdf)

bg = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index=True), geometry="geometry", crs="EPSG:4326")
bg.to_file(bg_path, driver="GPKG")
print(f"Saved {len(bg):,} block groups → {bg_path}\n")

print("Precompute complete. You can now submit per-state jobs.")
