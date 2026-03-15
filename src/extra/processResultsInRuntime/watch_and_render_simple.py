#!/usr/bin/env python3
"""
Watches a results directory for completed VTK files and renders them
offscreen with ParaView, then deletes the source file.

Usage:
    python watch_and_render.py [results_dir] [screenshots_dir]

Defaults:
    results_dir    = ./results
    screenshots_dir = ./screenshots
"""

import sys
import os
import re
import time
import shutil
import logging
import tempfile
import subprocess
from pathlib import Path
from threading import Thread
from queue import Queue, Empty

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler, FileCreatedEvent, FileMovedEvent
except ImportError:
    sys.exit(
        "watchdog not installed. Run: pip install watchdog --user\n"
        "Or if using a venv/conda: pip install watchdog"
    )

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

#SCRIPT_DIR       = Path(__file__).parent.resolve()
SCRIPT_DIR       = Path("/home/tomas/Work/sph-projects/fresh/tnl-sph/examples/WCSPH-BI/damBreak2D_WCSPH-BI/template")
VTK_PATTERN      = re.compile(r"fluid_[\d.]+_particles\.vtk$")
SETTLE_TIME      = 2.0          # seconds to wait after file appears before processing
WORKER_THREADS   = 2            # parallel pvpython renders (careful: RAM heavy)
PVPYTHON         = "pvpython"   # adjust if not on PATH, e.g. "/opt/paraview/bin/pvpython"
OFFSCREEN_FLAG   = "--force-offscreen-rendering"

STATES = {
    "pressure": SCRIPT_DIR / "results_pressure_vertical_legend_snapshot.pvsm",
    "velocity": SCRIPT_DIR / "results_velocity_vertical_legend_snapshot.pvsm",
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(SCRIPT_DIR / "render.log"),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# ParaView rendering (called inside a worker thread, runs pvpython subprocess)
# ---------------------------------------------------------------------------

PVPYTHON_SCRIPT = SCRIPT_DIR / "process_vtk.py"   # same script as before


def render_vtk(vtk_path: Path, screenshots_dir: Path) -> bool:
    """Render all states for a single VTK file. Returns True if all succeeded."""
    base = vtk_path.stem  # e.g. fluid_0.500006_particles
    all_ok = True

    for label, state_path in STATES.items():
        if not state_path.exists():
            log.warning("State file not found, skipping %s: %s", label, state_path)
            continue

        out_png = screenshots_dir / f"{base}_{label}.png"
        cmd = [
            PVPYTHON, OFFSCREEN_FLAG,
            str(PVPYTHON_SCRIPT),
            str(vtk_path),
            str(state_path),
            str(out_png),
        ]

        log.info("Rendering %-10s → %s", label, out_png.name)
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0 and out_png.exists():
            log.info("  ✓ %s", out_png.name)
        else:
            log.error("  ✗ %s render failed (exit %d)", label, result.returncode)
            if result.stdout:
                log.error("    stdout: %s", result.stdout.strip())
            if result.stderr:
                log.error("    stderr: %s", result.stderr.strip())
            all_ok = False

    return all_ok


# ---------------------------------------------------------------------------
# Worker thread — consumes the queue
# ---------------------------------------------------------------------------

def worker(queue: Queue, screenshots_dir: Path):
    while True:
        item = queue.get()
        if item is None:            # poison pill → shut down
            queue.task_done()
            break

        vtk_path: Path = item
        log.info("Processing: %s", vtk_path.name)

        # Small settle delay (file may still be flushed by the OS)
        time.sleep(SETTLE_TIME)

        if not vtk_path.exists():
            log.warning("File vanished before processing: %s", vtk_path)
            queue.task_done()
            continue

        success = render_vtk(vtk_path, screenshots_dir)

        if success:
            vtk_path.unlink()
            log.info("Deleted: %s", vtk_path.name)
        else:
            log.error("Keeping %s for inspection (render failed)", vtk_path.name)

        queue.task_done()


# ---------------------------------------------------------------------------
# Filesystem event handler
# ---------------------------------------------------------------------------

class VtkHandler(FileSystemEventHandler):
    """Enqueues VTK files when they appear (created or renamed into place)."""

    def __init__(self, queue: Queue):
        self.queue = queue
        self._seen: set[str] = set()

    def _enqueue(self, path: str):
        if VTK_PATTERN.search(path) and path not in self._seen:
            self._seen.add(path)
            self.queue.put(Path(path))
            log.info("Queued: %s (queue depth: %d)", Path(path).name, self.queue.qsize())

    def on_created(self, event):
        if not event.is_directory:
            self._enqueue(event.src_path)

    def on_moved(self, event):
        # Catches the atomic rename trick from C++ side
        if not event.is_directory:
            self._enqueue(event.dest_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    results_dir     = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./results")
    screenshots_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("./screenshots")

    results_dir.mkdir(parents=True, exist_ok=True)
    screenshots_dir.mkdir(parents=True, exist_ok=True)

    log.info("Watching   : %s", results_dir.resolve())
    log.info("Screenshots: %s", screenshots_dir.resolve())
    log.info("States     : %s", ", ".join(STATES))
    log.info("Workers    : %d", WORKER_THREADS)

    queue: Queue = Queue()

    # Start worker threads
    threads = []
    for _ in range(WORKER_THREADS):
        t = Thread(target=worker, args=(queue, screenshots_dir), daemon=True)
        t.start()
        threads.append(t)

    # Pick up any VTK files that arrived before we started watching
    for existing in sorted(results_dir.glob("fluid_*_particles.vtk")):
        queue.put(existing)
        log.info("Queued pre-existing file: %s", existing.name)

    # Start filesystem watcher
    handler  = VtkHandler(queue)
    observer = Observer()
    observer.schedule(handler, str(results_dir), recursive=False)
    observer.start()
    log.info("Watcher running — press Ctrl+C to stop.")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        log.info("Shutting down...")

    observer.stop()
    observer.join()

    # Drain the queue gracefully
    queue.join()
    for _ in threads:
        queue.put(None)     # poison pills
    for t in threads:
        t.join()

    log.info("All done.")


if __name__ == "__main__":
    main()
