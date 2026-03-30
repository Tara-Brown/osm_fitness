import subprocess
import os

# Full state names mapped to abbreviations, contiguous US + DC (no AK, HI)
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
    'VT': 'Vermont', 'VA': 'Virginia',  'WV': 'West_Virginia',
    'WI': 'Wisconsin', 'WY': 'Wyoming',
}

for abbr, full_name in STATE_NAMES.items():
    job_name = f"cel_{abbr}"
    output_log = f"logs/cel_{abbr}_%j.out"
    error_log  = f"logs/cel_{abbr}_%j.err"

    output_file = f"lookups/ubermedia_lookup_{full_name}.parquet"
    
    if os.path.exists(output_file):
        print(f"Skipping {abbr} ({full_name}) — output file already exists.")
        continue

    cmd = [
        "sbatch",
        f"--job-name={job_name}",
        f"--output={output_log}",
        f"--error={error_log}",
        "--time=48:00:00",
        "--mem=128G",
        "--cpus-per-task=4",
        "--wrap",
        f"python income_and_education.py --state_abbr {abbr} --state_name {full_name}"
    ]

    print(f"Submitting {abbr} ({full_name})...")
    subprocess.run(cmd)
