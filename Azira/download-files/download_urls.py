import os
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from turtle import end_fill
from urllib.parse import urlparse, parse_qs, unquote
import re
import requests
from datetime import datetime

# ---------------- config ----------------
URLS_FILE = "/mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira/get-urls/urls.txt"
OUT_DIR = "/mnt/beegfs/hellgate/home/vc149353/azira_downloads"
MAX_WORKERS = 10
TIMEOUT = 60          # seconds per request
RETRIES = 5
CHUNK_SIZE = 1024 * 1024  # 1 MB
DOWNLOAD_COUNT = 0
# ----------------------------------------

os.makedirs(OUT_DIR, exist_ok=True)
lock = threading.Lock()


def filename_from_url(url):
    qs = parse_qs(urlparse(url).query)
    cd = qs.get("response-content-disposition", [None])[0]
    if not cd:
        return None
    cd = unquote(cd)  # attachment; filename="..."
    m = re.search(r'filename="?([^"]+)"?', cd)
    return m.group(1) if m else None


def download_one(url: str):
    global DOWNLOAD_COUNT       
    name = filename_from_url(url)
    out_path = os.path.join(OUT_DIR, name)

    if os.path.exists(out_path):
        with lock:
            print(f"[skip] {name}")
        return

    for attempt in range(1, RETRIES + 1):
        try:
            with requests.get(
                url,
                stream=True,
                timeout=TIMEOUT,
                allow_redirects=True,
            ) as r:
                r.raise_for_status()

                # Honor Content-Disposition if present
                cd = r.headers.get("Content-Disposition", "")
                if "filename=" in cd:
                    fname = cd.split("filename=")[-1].strip("\"'")
                    out_path = os.path.join(OUT_DIR, fname)
                    if os.path.exists(out_path):
                        with lock:
                            print(f"[skip] {fname}")
                        return

                tmp = out_path + ".part"
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)

                os.replace(tmp, out_path)

                with lock:
                    DOWNLOAD_COUNT += 1
                    print(f"[ok]   {os.path.basename(out_path)}")
                return

        except Exception as e:
            if attempt == RETRIES:
                with lock:
                    print(f"[FAIL] {name}: {e}")
            else:
                time.sleep(2 ** attempt)


def main():
    start = datetime.now()
    print(f"Starting downloads at {start.isoformat()}...")

    if not os.path.exists(URLS_FILE):
        print(f"missing {URLS_FILE}")
        sys.exit(1)

    with open(URLS_FILE) as f:
        urls = [l.strip() for l in f if l.strip()]

    if not urls:
        print("no URLs found")
        return

    print(f"Attempting {len(urls)} files with {MAX_WORKERS} threads")

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(download_one, url) for url in urls]
        for _ in as_completed(futures):
            pass

    end = datetime.now()
    print(f"Downloaded {DOWNLOAD_COUNT} files.")
    print(f"Finished at {end.isoformat()}, duration: {end - start}")
    

if __name__ == "__main__":
    main()