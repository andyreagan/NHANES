"""Download NHANES data files from a manifest CSV.

Reads the manifest produced by scrape_nhanes_manifest.py and downloads all
XPT files that are not already present locally. Downloads are parallelized
with a configurable number of workers and include retry logic.

Usage:
    # Download everything:
    uv run python -m src.download.download_nhanes_files

    # Download only specific cycles:
    uv run python -m src.download.download_nhanes_files --cycles 2003-2004 2005-2006

    # Download only specific components:
    uv run python -m src.download.download_nhanes_files --components Laboratory Demographics

    # Dry run (show what would be downloaded):
    uv run python -m src.download.download_nhanes_files --dry-run

    # Control parallelism:
    uv run python -m src.download.download_nhanes_files --workers 4
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

MANIFEST_PATH = Path("data/nhanes_file_manifest.csv")
MAX_RETRIES = 3
RETRY_DELAY = 2.0  # seconds between retries
USER_AGENT = "NHANES-download-script/1.0"


def _file_exists(dest: Path) -> Path | None:
    """Check if a file exists, trying both .XPT and .xpt extensions.

    Some files were previously downloaded with lowercase .xpt.
    Returns the actual path if found, None otherwise.
    """
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    # Try alternate case for extension
    alt = dest.with_suffix(dest.suffix.swapcase())
    if alt.exists() and alt.stat().st_size > 100:
        return alt
    return None


def download_file(url: str, dest: Path, retries: int = MAX_RETRIES) -> tuple[str, bool, str]:
    """Download a single file with retry logic.

    Returns:
        (url, success, message)
    """
    existing = _file_exists(dest)
    if existing is not None:
        return (url, True, f"already exists ({existing.stat().st_size:,d} bytes)")

    dest.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=60) as resp:
                data = resp.read()
            dest.write_bytes(data)
            return (url, True, f"downloaded ({len(data):,d} bytes)")
        except (HTTPError, URLError, TimeoutError, OSError) as e:
            if attempt < retries:
                time.sleep(RETRY_DELAY * attempt)
            else:
                return (url, False, f"FAILED after {retries} attempts: {e}")

    return (url, False, "FAILED: unknown error")


def load_manifest(
    manifest_path: Path,
    cycles: list[str] | None = None,
    components: list[str] | None = None,
    files: list[str] | None = None,
) -> list[dict]:
    """Load and optionally filter the manifest CSV.

    Args:
        manifest_path: Path to the manifest CSV.
        cycles: Only include these cycles (e.g., ["2003-2004", "2005-2006"]).
        components: Only include these components (e.g., ["Laboratory"]).
        files: Only include files whose base name (without suffix) starts with
               one of these prefixes (e.g., ["DEMO", "BIOPRO", "L40"]).
               This matches "DEMO" against "DEMO", "DEMO_B", "DEMO_C", etc.
    """
    if not manifest_path.exists():
        print(
            f"Manifest not found at {manifest_path}.\n"
            "Run first:  uv run python -m src.download.scrape_nhanes_manifest",
            file=sys.stderr,
        )
        sys.exit(1)

    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        records = list(reader)

    if cycles:
        cycle_set = set(cycles)
        records = [r for r in records if r["cycle"] in cycle_set]

    if components:
        comp_set = set(components)
        records = [r for r in records if r["component"] in comp_set]

    if files:
        file_set = set(files)

        def _matches_file(record: dict) -> bool:
            fname = record["file_name"]
            # Match exact name or name without cycle suffix (e.g., DEMO_C → DEMO)
            if fname in file_set:
                return True
            # Strip the suffix letter: "BIOPRO_D" → "BIOPRO", "L40_C" → "L40"
            # Also handle P_ prefix: "P_DEMO" → "DEMO"
            base = fname
            if base.startswith("P_"):
                base = base[2:]
            # Remove trailing _X suffix (single letter)
            parts = base.rsplit("_", 1)
            if len(parts) == 2 and len(parts[1]) == 1 and parts[1].isalpha():
                base = parts[0]
            return base in file_set

        records = [r for r in records if _matches_file(r)]

    return records


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download NHANES XPT files from CDC.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=MANIFEST_PATH,
        help=f"Path to manifest CSV (default: {MANIFEST_PATH})",
    )
    parser.add_argument(
        "--cycles",
        nargs="+",
        default=None,
        help="Only download these cycles (e.g., 2003-2004 2005-2006)",
    )
    parser.add_argument(
        "--components",
        nargs="+",
        default=None,
        help="Only download these components (e.g., Laboratory Demographics)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel download workers (default: 4)",
    )
    parser.add_argument(
        "--files",
        nargs="+",
        default=None,
        help="Only download files matching these base names (e.g., DEMO BIOPRO L40)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without actually downloading",
    )
    args = parser.parse_args()

    records = load_manifest(args.manifest, args.cycles, args.components, args.files)

    if not records:
        print("No files match the specified filters.")
        return

    # Partition into already-downloaded and to-download
    to_download = []
    already_have = []
    for r in records:
        dest = Path(r["local_path"])
        if _file_exists(dest) is not None:
            already_have.append(r)
        else:
            to_download.append(r)

    total = len(records)
    print(f"Manifest: {total} files total")
    print(f"  Already downloaded: {len(already_have)}")
    print(f"  To download:        {len(to_download)}")

    if args.dry_run:
        print("\nDry run — files that would be downloaded:")
        for r in to_download:
            print(f"  {r['local_path']}  ←  {r['url']}")
        return

    if not to_download:
        print("\nAll files already downloaded!")
        return

    print(f"\nDownloading {len(to_download)} files with {args.workers} workers...\n")

    succeeded = 0
    failed = 0

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_file, r["url"], Path(r["local_path"])): r for r in to_download
        }
        for i, future in enumerate(as_completed(futures), 1):
            record = futures[future]
            url, success, message = future.result()
            status = "✓" if success else "✗"
            if success:
                succeeded += 1
            else:
                failed += 1
            print(
                f"  [{i}/{len(to_download)}] {status} {record['cycle']}/{record['file_name']}: "
                f"{message}"
            )

    print(f"\nDone: {succeeded} succeeded, {failed} failed out of {len(to_download)} attempted")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
