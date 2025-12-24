#from playwright.sync_api import sync_playwright
from playwright.async_api import async_playwright, expect, TimeoutError
import json
from pathlib import Path
import time
from datetime import date, timedelta
import asyncio
from urllib.parse import urljoin
import os
from urllib.parse import urlparse, parse_qs, unquote
import re
from pydoc import cli
import subprocess

# login info 
user, pw = 'katrina.mullan@mso.umt.edu', 'wvp6rau6rqb_bwy1EYT!'
url = 'https://pinnacle.azira.com/'


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
    await page.goto(url) 
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
            for name_item in job_name.split("_")[::-1]:
                if not name_item in ['a','b']:
                    state = name_item
                    break
            if state not in restrict_state:
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

    try:
        await login(page, user, pw)
        print("logged in")

        await click_load_more_button(page)
        print("loaded all entries")

        my_states = ['Colorado', 'Oklahoma', 'Utah', 'California', 'Nevada', "New", "Arizona", "North"]
        all_links = await get_links(page, my_states)

        print(f"found {len(all_links)} links")

        with open("urls.txt", "w") as f:
            for item in all_links:
                url = item[3]
                f.write(f"{url}\n")
        print("urls written to urls.txt")
    finally:
        await browser_handle["context"].close()
        await browser_handle["browser"].close()
        await browser_handle["playwright"].stop() 
    

    # transfer to HPC

    local_file = "urls.txt"
    remote =  "vc149353@um:/mnt/beegfs/hellgate/home/vc149353/osm_fitness/Azira/urls.txt"

    subprocess.run(
        [
            "scp",
            "-i", "/home/vince/.ssh/um_vpn",
            local_file,
            remote
        ],
        check=True
    )
    print("urls.txt transferred to hellgate")

if __name__ == "__main__":  
    asyncio.run(main())