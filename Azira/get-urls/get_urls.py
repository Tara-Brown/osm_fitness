#from playwright.sync_api import sync_playwright
from playwright.async_api import async_playwright, expect, TimeoutError
from pathlib import Path
from datetime import datetime
import asyncio
from urllib.parse import urlparse, parse_qs, unquote
import re
import subprocess


# ---------------CONFIG----------------------
# Azira login info 
USER, PW = 'katrina.mullan@mso.umt.edu', 'wvp6rau6rqb_bwy1EYT!'
URL = 'https://pinnacle.azira.com/'

# states to download
MY_STATES = ['Colorado', 'Oklahoma', 'Utah', 'California', 'Nevada', "Arizona", "North", "New_Mexico", "New"] # north is north carolina
TWO_WORD_STATES = ['New', 'South', 'Rhode', 'North']
# my_states = None   # all states

# url file locations
LOCAL_URL_FILE = '/home/vince/Documents/SmartFires/osm_fitness/Azira/get-urls/urls.txt'
REMOTE_URL_FILE =  "vc149353@um:/mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira/get-urls/urls.txt"

# ssh key for hellgate
HELLGATE_SSH_KEY = "/home/vince/.ssh/um_vpn"
#----------------------------------------------

async def start_browser(headless=False):
    playwright = await async_playwright().start()

    browser = await playwright.chromium.launch(headless=headless)

    context = await browser.new_context(
        geolocation={"latitude": 46.8721, "longitude": -113.9940},
        permissions=["geolocation"],
    )

    page = await context.new_page()

    return {
        "playwright": playwright,
        "browser": browser,
        "context": context,
        "page": page,
    }

async def login(page, user, pw):
    await page.goto(URL) 
    await page.get_by_role("textbox", name="Email Address").fill(user)
    await page.get_by_role("textbox", name="Password").fill(pw)
    await page.get_by_role("button", name="Log in").click()

    # wait for page load
    await page.wait_for_url("**/home/feeds/datasets**")
    load_more = page.get_by_role("button", name="Load More")
    await load_more.wait_for(state="visible")

# get links on a page
async def get_links(page, restrict_state=None):
    links = []      

    rows = page.locator("table tr")
    row_count = await rows.count()

    for i in range(row_count):
        links_cell = rows.nth(i).locator("td:nth-child(4)")

        # skip empty cells
        if await links_cell.count() == 0:
            continue    

        # get row metadata
        job_id = await rows.nth(i).locator("td:nth-child(1)").inner_text()   
        job_name = await rows.nth(i).locator("td:nth-child(2)").inner_text()

        if restrict_state:
            job_name_items = job_name.split("_")
            state = job_name_items[0] if len(job_name_items) < 2 else ("_".join(job_name_items[-2:]) if job_name_items[-2] in TWO_WORD_STATES else job_name_items[-1])
            if not (state in restrict_state):
               continue 

        anchors = links_cell.locator("a")
        link_count = await anchors.count()
        for j in range(link_count):
            link = anchors.nth(j)
            text = (await link.inner_text()).strip()   # "Pin Dataset" or "Expanded Standard Dataset"

            if "Pin Dataset" in text or "Expanded Standard" in text:
                text_link = await link.get_attribute("href")
                if text_link:
                    links.append((job_id, job_name, text, text_link))

    return links


# filter for already downloaded files
def filter_downloaded(all_links):

    new_links = []

    for a in all_links:
        filename = filename_from_url(a[3])
        if not Path(f"azira_downloads/{filename}").exists():
            print(f"adding {filename} to download list")
            new_links.append(a)

    return new_links    

def _safe(s: str) -> str:
    return "".join(c if c.isalnum() or c in "._- " else "_" for c in s).strip()

def filename_from_url(url):
    qs = parse_qs(urlparse(url).query)
    cd = qs.get("response-content-disposition", [None])[0]
    if not cd:
        return None
    cd = unquote(cd)  # attachment; filename="..."
    m = re.search(r'filename="?([^"]+)"?', cd)
    return m.group(1) if m else None


async def click_load_more_button(page):   
    while True:
        load_more = page.get_by_role("button", name="Load More")

        if await load_more.is_visible():
            await load_more.click()
        else:
            break
        # Wait for new rows to load before next iteration
        await page.wait_for_timeout(1000) 

async def main():
    # login 
    browser_handle = await start_browser(headless=True)
    page = browser_handle["page"]

    print(datetime.now())

    try:
        # get all links from the page
        print("Logging in...")
        await login(page, USER, PW)
        print("Loading all rows...")
        await click_load_more_button(page)
        print("Extracting download links...")
        all_links = await get_links(page, MY_STATES)

        # write to file
        with open("urls.txt", "w") as f:
            for item in all_links:
                curr_url = item[3]
                f.write(f"{curr_url}\n")
        print(f"{len(all_links)} urls written to urls.txt")

    finally:
        await browser_handle["context"].close()
        await browser_handle["browser"].close()
        await browser_handle["playwright"].stop() 
    

    # transfer urls.txt to HPC

    subprocess.run(
        [
            "scp",
            "-i", HELLGATE_SSH_KEY,
            LOCAL_URL_FILE,
            REMOTE_URL_FILE
        ],
        check=True
    )
    print("urls.txt transferred to hellgate")

if __name__ == "__main__":  
    asyncio.run(main())