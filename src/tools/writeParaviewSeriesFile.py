import os
import json
import shutil
import re
from pathlib import Path
from typing import List, Dict, Any, Tuple
from collections import defaultdict

TIME_RE = re.compile(r"_([0-9]*\.[0-9]+)_")        # e.g. _0.250003_
SUBDOMAIN_RE = re.compile(r"(subdomain\d+)")       # subdomain0, subdomain1, ...

def extract_time(name: str) -> float:
    m = TIME_RE.search(name)
    return float(m.group(1)) if m else 0.0

def extract_subdomain(name: str) -> str:
    m = SUBDOMAIN_RE.search(name)
    return m.group(1) if m else ""

def series_key(prefix: str, subdomain: str) -> str:
    return f"{prefix}_{subdomain}" if subdomain else prefix

# build series funciton
def generate_series(root_directory: str = "./") -> None:
    """
    1. Scan root directory for *.vtk files
    2. Group by (prefix, subdomain) → independent series
    3. For each series:
         • Create results_data/<series_id>/
         • MOVE original .vtk files there (no rename)
         • Write <series_id>.series with:
           {
               "file-series-version": "1.0",
               "files": [
                   { "name": "results_data/<series_id>/<original_name>.vtk", "time": X.XXXXXX },
                   ...
               ]
           }
    """
    root_path = Path(root_directory).resolve()

    # 1. Find all .vtk files (non-recursive)
    vtk_files = [p for p in root_path.iterdir() if p.suffix.lower() == ".vtk"]
    if not vtk_files:
        print("No .vtk files found.")
        return

    # 2. Group: (prefix, subdomain) → list of Path
    groups: Dict[Tuple[str, str], List[Path]] = defaultdict(list)

    for f in vtk_files:
        name = f.name
        parts = name.split("_", 2)
        #if len(parts) < 3:
        #    continue
        prefix = parts[0]
        subdomain = extract_subdomain(name)
        groups[(prefix, subdomain)].append(f)

    print(f"Detected {len(groups)} independent series.")

    # 3. Process each series
    for (prefix, subdomain), files in groups.items():
        series_id = series_key(prefix, subdomain)
        series_dir = root_path / "results_data" / series_id
        series_dir.mkdir(parents=True, exist_ok=True)

        # sort by embedded time
        files_sorted = sorted(files, key=lambda p: extract_time(p.name))

        series_entries: List[Dict[str, Any]] = []
        moved = 0

        for src_path in files_sorted:
            dest_path = series_dir / src_path.name

            # move only if not already there (same content)
            if dest_path.exists() and src_path.read_bytes() == dest_path.read_bytes():
                pass
            else:
                shutil.move(str(src_path), str(dest_path))
                moved += 1

            rel_path = os.path.relpath(dest_path, root_path).replace(os.sep, "/")
            time_val = extract_time(src_path.name)

            series_entries.append({
                "name": rel_path,
                "time": time_val
            })

        # 4. Write .series file
        series_file = root_path / f"{series_id}.vtk.series"
        data = {
            "file-series-version": "1.0",
            "files": series_entries
        }
        with series_file.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)

        print(f"[{series_id}] {moved} files → {series_dir} | {series_file.name} written")

if __name__ == "__main__":

    generate_series()
