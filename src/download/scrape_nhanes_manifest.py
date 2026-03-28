"""Scrape the NHANES data page to build a manifest of all available XPT files.

This script crawls the CDC NHANES data pages for every available cycle and
component (Demographics, Dietary, Examination, Laboratory, Questionnaire),
extracts all .xpt download links, and writes a CSV manifest.

Usage:
    uv run python -m src.download.scrape_nhanes_manifest

Output:
    data/nhanes_file_manifest.csv — one row per downloadable XPT file with columns:
        cycle, component, file_name, url, local_path
"""

from __future__ import annotations

import csv
import re
import sys
import time
from pathlib import Path
from urllib.request import Request, urlopen

# ---------------------------------------------------------------------------
# NHANES cycles and components to scrape
# ---------------------------------------------------------------------------

# All public continuous NHANES cycles, plus the pre-pandemic combined cycle.
# 2019-2020 data collection was interrupted by COVID and is not released
# separately; the usable pre-pandemic portion is folded into 2017-2020.
# 2025-2026 is listed on CDC but has no data files yet.
CYCLES: list[str] = [
    "1999-2000",
    "2001-2002",
    "2003-2004",
    "2005-2006",
    "2007-2008",
    "2009-2010",
    "2011-2012",
    "2013-2014",
    "2015-2016",
    "2017-2018",
    "2017-2020",  # pre-pandemic combined (P_ prefix files)
    "2021-2023",
]

COMPONENTS: list[str] = [
    "Demographics",
    "Dietary",
    "Examination",
    "Laboratory",
    "Questionnaire",
]

BASE_URL = "https://wwwn.cdc.gov"
DATA_PAGE_URL = BASE_URL + "/nchs/nhanes/search/datapage.aspx?Component={component}&Cycle={cycle}"

# Regex to pull .xpt links from the HTML
XPT_LINK_RE = re.compile(r'href="(/Nchs/Data/Nhanes/[^"]*\.xpt)"', re.IGNORECASE)


def cycle_to_local_dir(cycle: str) -> str:
    """Map cycle label to local directory name under data/raw/."""
    return cycle


def url_to_local_path(url_path: str, cycle: str) -> str:
    """Convert a CDC URL path to a local file path under data/raw/.

    CDC URLs look like:
        /Nchs/Data/Nhanes/Public/2003/DataFiles/DEMO_C.xpt
        /Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.xpt

    We store them as:
        data/raw/2003-2004/DEMO_C.XPT
        data/raw/2017-2020/P_DEMO.XPT
    """
    file_name = url_path.rsplit("/", 1)[-1]
    # Normalize extension to uppercase to match existing convention
    stem, _ext = file_name.rsplit(".", 1)
    file_name = f"{stem}.XPT"
    local_dir = cycle_to_local_dir(cycle)
    return f"data/raw/{local_dir}/{file_name}"


def scrape_cycle_component(cycle: str, component: str) -> list[dict]:
    """Scrape a single cycle+component page and return list of file records."""
    url = DATA_PAGE_URL.format(component=component, cycle=cycle)
    req = Request(url, headers={"User-Agent": "NHANES-download-script/1.0"})
    try:
        with urlopen(req, timeout=30) as resp:
            html = resp.read().decode("utf-8", errors="replace")
    except Exception as e:
        print(f"  ⚠ Error fetching {cycle}/{component}: {e}", file=sys.stderr)
        return []

    matches = XPT_LINK_RE.findall(html)
    records = []
    for url_path in matches:
        file_name = url_path.rsplit("/", 1)[-1]
        stem = file_name.rsplit(".", 1)[0]
        full_url = BASE_URL + url_path
        local_path = url_to_local_path(url_path, cycle)
        records.append(
            {
                "cycle": cycle,
                "component": component,
                "file_name": stem,
                "url": full_url,
                "local_path": local_path,
            }
        )
    return records


def main() -> None:
    output_path = Path("data/nhanes_file_manifest.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_records: list[dict] = []
    for cycle in CYCLES:
        for component in COMPONENTS:
            print(f"Scraping {cycle} / {component}...", end=" ", flush=True)
            records = scrape_cycle_component(cycle, component)
            print(f"{len(records)} files")
            all_records.extend(records)
            # Be polite to CDC servers
            time.sleep(0.3)

    # Deduplicate by URL (same file can theoretically appear in multiple pages)
    seen_urls: set[str] = set()
    unique_records: list[dict] = []
    for r in all_records:
        if r["url"] not in seen_urls:
            seen_urls.add(r["url"])
            unique_records.append(r)

    # Sort for stable output
    unique_records.sort(key=lambda r: (r["cycle"], r["component"], r["file_name"]))

    # Write CSV
    fieldnames = ["cycle", "component", "file_name", "url", "local_path"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(unique_records)

    print(f"\nWrote {len(unique_records)} files to {output_path}")
    print(f"Cycles: {sorted(set(r['cycle'] for r in unique_records))}")
    print(f"Components: {sorted(set(r['component'] for r in unique_records))}")


if __name__ == "__main__":
    main()
