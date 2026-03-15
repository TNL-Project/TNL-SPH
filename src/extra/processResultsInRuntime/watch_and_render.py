#!/usr/bin/env python3
"""
Generic offscreen ParaView watcher.

Usage:
    python watch_and_render.py --config path/to/render_config.toml [overrides]
    python watch_and_render.py --results ./results --screenshots ./shots \
                                --states pressure:state_p.pvsm velocity:state_v.pvsm

Run with --help for all options.
"""

import sys
import re
import time
import logging
import argparse
import subprocess
from pathlib import Path
from threading import Thread
from queue import Queue

# Optional TOML support (stdlib in Python 3.11+, else install tomli)
def _load_toml(path: Path) -> dict:
    try:
        import tomllib                          # Python 3.11+
        with open(path, "rb") as f:
            return tomllib.load(f)
    except ImportError:
        try:
            import tomli as tomllib             # pip install tomli
            with open(path, "rb") as f:
                return tomllib.load(f)
        except ImportError:
            sys.exit(
                "TOML support requires Python 3.11+, or: pip install tomli --user"
            )

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler
except ImportError:
    sys.exit("watchdog not installed. Run: pip install watchdog --user")


# ---------------------------------------------------------------------------
# Defaults (overridden by config file and/or CLI args)
# ---------------------------------------------------------------------------

DEFAULTS = {
    "results_dir":      "./results",
    "screenshots_dir":  "./screenshots",
    "script_dir":       None,           # directory containing .pvsm files; defaults to config file dir
    "pvpython":         "pvpython",
    "offscreen_flag":   "--force-offscreen-rendering",
    "settle_time":      2.0,
    "worker_threads":   2,
    "patterns":         [
        r"fluid_[\d.]+_particles\.vtk$",
        r"boundary_[\d.]+_particles\.vtk$",
    ],
    "states":           {},             # label -> pvsm path (relative to script_dir or absolute)
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generic offscreen ParaView watcher/renderer.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use a config file (recommended)
  python watch_and_render.py --config ./sim_A/render_config.toml

  # Override results dir from config
  python watch_and_render.py --config ./render_config.toml --results /scratch/run42/results

  # Fully CLI-driven (no config file)
  python watch_and_render.py \\
      --results ./results \\
      --screenshots ./screenshots \\
      --script-dir ./states \\
      --states pressure:results_pressure.pvsm velocity:results_velocity.pvsm \\
      --patterns "fluid_[\\d.]+_particles\\.vtk$" "boundary_[\\d.]+_particles\\.vtk$"
""",
    )

    p.add_argument("--config",        metavar="FILE",    help="Path to render_config.toml")
    p.add_argument("--results",       metavar="DIR",     help="Directory to watch for VTK files")
    p.add_argument("--screenshots",   metavar="DIR",     help="Directory to save PNG screenshots")
    p.add_argument("--script-dir",    metavar="DIR",     help="Base directory for .pvsm state files")
    p.add_argument("--pvpython",      metavar="CMD",     help="pvpython executable (default: pvpython)")
    p.add_argument("--workers",       metavar="N", type=int, help="Parallel render threads")
    p.add_argument("--settle-time",   metavar="S", type=float, help="Seconds to wait after file appears")
    p.add_argument(
        "--states", nargs="+", metavar="LABEL:FILE",
        help="State files as label:filename pairs, e.g. pressure:state_p.pvsm"
    )
    p.add_argument(
        "--patterns", nargs="+", metavar="REGEX",
        help="Filename regex patterns to match (default: fluid + boundary particles)"
    )
    p.add_argument(
        "--process-script", metavar="FILE",
        help="Path to process_vtk.py (default: process_vtk.py next to this script)"
    )

    return p.parse_args()


# ---------------------------------------------------------------------------
# Config assembly: defaults → toml → CLI  (later sources win)
# ---------------------------------------------------------------------------

def build_config(args: argparse.Namespace) -> dict:
    cfg = dict(DEFAULTS)

    # Layer 1: TOML config file
    config_dir = None
    if args.config:
        config_path = Path(args.config).resolve()
        if not config_path.exists():
            sys.exit(f"Config file not found: {config_path}")
        toml = _load_toml(config_path)
        config_dir = config_path.parent
        cfg.update(toml)

    # Layer 2: CLI overrides
    if args.results:        cfg["results_dir"]     = args.results
    if args.screenshots:    cfg["screenshots_dir"] = args.screenshots
    if args.script_dir:     cfg["script_dir"]      = args.script_dir
    if args.pvpython:       cfg["pvpython"]        = args.pvpython
    if args.workers:        cfg["worker_threads"]  = args.workers
    if args.settle_time:    cfg["settle_time"]     = args.settle_time
    if args.patterns:       cfg["patterns"]        = args.patterns

    if args.states:
        parsed = {}
        for item in args.states:
            if ":" not in item:
                sys.exit(f"--states entries must be label:file, got: {item!r}")
            label, path = item.split(":", 1)
            parsed[label] = path
        cfg["states"] = parsed

    # Resolve script_dir: CLI > TOML > directory of config file > cwd
    if cfg["script_dir"] is None:
        cfg["script_dir"] = str(config_dir) if config_dir else "."
    cfg["script_dir"] = Path(cfg["script_dir"]).resolve()

    # Resolve state paths relative to script_dir
    script_dir = cfg["script_dir"]
    resolved_states = {}
    for label, pvsm in cfg["states"].items():
        p = Path(pvsm)
        resolved = p if p.is_absolute() else (script_dir / p)
        if not resolved.exists():
            sys.exit(f"State file for '{label}' not found: {resolved}")
        resolved_states[label] = resolved
    cfg["states"] = resolved_states

    if not cfg["states"]:
        sys.exit("No states defined. Use --states or add [states] to your config file.")

    # Compile patterns
    cfg["patterns"] = [re.compile(pat) for pat in cfg["patterns"]]

    # Resolve process_vtk.py
    if args.process_script:
        pv_script = Path(args.process_script).resolve()
    else:
        pv_script = Path(__file__).parent / "process_vtk.py"
    if not pv_script.exists():
        sys.exit(f"process_vtk.py not found: {pv_script}")
    cfg["process_script"] = pv_script

    return cfg


# ---------------------------------------------------------------------------
# Logging setup (called after config is known)
# ---------------------------------------------------------------------------

def setup_logging(screenshots_dir: Path):
    log_path = screenshots_dir / "render.log"
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path),
        ],
    )


log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_vtk(vtk_path: Path, cfg: dict) -> bool:
    screenshots_dir = Path(cfg["screenshots_dir"])
    base = vtk_path.stem
    all_ok = True

    for label, state_path in cfg["states"].items():
        out_png = screenshots_dir / f"{base}_{label}.png"
        cmd = [
            cfg["pvpython"],
            cfg["offscreen_flag"],
            str(cfg["process_script"]),
            str(vtk_path),
            str(state_path),
            str(out_png),
        ]

        log.info("  Rendering %-12s → %s", label, out_png.name)
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0 and out_png.exists():
            log.info("    ✓ saved %s", out_png.name)
        else:
            log.error("    ✗ %s failed (exit %d)", label, result.returncode)
            for line in (result.stdout + result.stderr).strip().splitlines():
                log.error("      %s", line)
            all_ok = False

    return all_ok


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def worker(queue: Queue, cfg: dict):
    while True:
        item = queue.get()
        if item is None:
            queue.task_done()
            break

        vtk_path: Path = item
        log.info("Processing: %s", vtk_path.name)
        time.sleep(cfg["settle_time"])

        if not vtk_path.exists():
            log.warning("File vanished before processing: %s", vtk_path)
            queue.task_done()
            continue

        success = render_vtk(vtk_path, cfg)

        if success:
            vtk_path.unlink()
            log.info("Deleted: %s", vtk_path.name)
        else:
            log.error("Keeping %s for inspection", vtk_path.name)

        queue.task_done()


# ---------------------------------------------------------------------------
# Filesystem watcher
# ---------------------------------------------------------------------------

class VtkHandler(FileSystemEventHandler):
    def __init__(self, queue: Queue, patterns: list):
        self.queue    = queue
        self.patterns = patterns
        self._seen: set = set()

    def _matches(self, path: str) -> bool:
        name = Path(path).name
        return any(pat.search(name) for pat in self.patterns)

    def _enqueue(self, path: str):
        if self._matches(path) and path not in self._seen:
            self._seen.add(path)
            self.queue.put(Path(path))
            log.info("Queued: %s  (queue depth: %d)", Path(path).name, self.queue.qsize())

    def on_created(self, event):
        if not event.is_directory:
            self._enqueue(event.src_path)

    def on_moved(self, event):           # catches atomic rename from C++
        if not event.is_directory:
            self._enqueue(event.dest_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    cfg  = build_config(args)

    results_dir     = Path(cfg["results_dir"]).resolve()
    screenshots_dir = Path(cfg["screenshots_dir"]).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    screenshots_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(screenshots_dir)

    log.info("=== ParaView watcher starting ===")
    log.info("Results dir  : %s", results_dir)
    log.info("Screenshots  : %s", screenshots_dir)
    log.info("Script dir   : %s", cfg['script_dir'])
    log.info("States       : %s", {k: str(v) for k, v in cfg['states'].items()})
    log.info("Patterns     : %s", [p.pattern for p in cfg['patterns']])
    log.info("Workers      : %d", cfg['worker_threads'])
    log.info("pvpython     : %s", cfg['pvpython'])

    queue: Queue = Queue()

    threads = [
        Thread(target=worker, args=(queue, cfg), daemon=True)
        for _ in range(cfg["worker_threads"])
    ]
    for t in threads:
        t.start()

    # Catch files that already exist before the observer starts
    for existing in sorted(results_dir.iterdir()):
        if not existing.is_file():
            continue
        name = existing.name
        if any(pat.search(name) for pat in cfg["patterns"]):
            queue.put(existing)
            log.info("Queued pre-existing: %s", name)

    handler  = VtkHandler(queue, cfg["patterns"])
    observer = Observer()
    observer.schedule(handler, str(results_dir), recursive=False)
    observer.start()
    log.info("Watching for new files — press Ctrl+C to stop.")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        log.info("Shutting down...")

    observer.stop()
    observer.join()
    queue.join()

    for _ in threads:
        queue.put(None)
    for t in threads:
        t.join()

    log.info("=== Watcher stopped ===")


if __name__ == "__main__":
    main()
