#!/usr/bin/env python3
"""
Generic offscreen ParaView watcher.

Usage:
    python watch_and_render.py --config path/to/render_config.toml [overrides]
    python watch_and_render.py --results ./results --screenshots ./shots \
                                --states pressure:state_p.pvsm velocity:state_v.pvsm

Run with --help for all options.
"""

from __future__ import annotations

import sys
import re
import time
import logging
import argparse
import subprocess
import threading
from pathlib import Path
from threading import Thread
from queue import Queue

# Import template_to_regex from process_vtk (same directory)
sys.path.insert(0, str(Path(__file__).parent))
from process_vtk import template_to_regex  # noqa: E402

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
    "states":               {},             # label -> pvsm path (relative to script_dir or absolute)
    "particle_filename_template": "fluid_{time}_particles.vtk",
    "expected_files_per_timestep": 1,       # 1 = no grouping (legacy single-resolution behaviour)
    "group_timeout":        60.0,           # seconds before flushing an incomplete group
    "project_dir":          None,           # project root for resolving {project_dir} in .pvsm
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
    p.add_argument(
        "--template", metavar="TEMPLATE",
        help="Filename template for patching, e.g. 'fluid_subdomain{id}_{time}_particles.vtk'"
    )
    p.add_argument(
        "--project-dir", metavar="DIR",
        help="Project root for resolving {project_dir} placeholders in .pvsm files"
    )
    p.add_argument(
        "--expected-files", metavar="N", type=int,
        help="Number of VTK files per timestep (1 = no grouping, default: 1)"
    )
    p.add_argument(
        "--group-timeout", metavar="S", type=float,
        help="Seconds before flushing an incomplete group (default: 60)"
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
    if args.template:       cfg["particle_filename_template"] = args.template
    if args.project_dir:    cfg["project_dir"]     = args.project_dir
    if args.expected_files: cfg["expected_files_per_timestep"] = args.expected_files
    if args.group_timeout:  cfg["group_timeout"]   = args.group_timeout

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

    # Compile template regex for timestep extraction and grouping
    cfg["_tmpl_regex"] = template_to_regex(cfg["particle_filename_template"])

    # Resolve project_dir relative to config file directory (or cwd)
    if cfg.get("project_dir"):
        p = Path(cfg["project_dir"])
        if not p.is_absolute():
            base = config_dir if config_dir else Path(".")
            p = (base / p).resolve()
        else:
            p = p.resolve()
        cfg["project_dir"] = str(p)

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

def render_vtk(vtk_paths: list[Path], cfg: dict) -> bool:
    screenshots_dir = Path(cfg["screenshots_dir"])
    tmpl_regex = cfg["_tmpl_regex"]

    # Derive a common base name from the timestep shared by all files in the group
    m = tmpl_regex.search(vtk_paths[0].name)
    if m and "time" in m.groupdict():
        base = f"fluid_{m.group('time')}_particles"
    else:
        base = vtk_paths[0].stem

    all_ok = True

    for label, state_path in cfg["states"].items():
        out_png = screenshots_dir / f"{base}_{label}.png"
        cmd = [
            cfg["pvpython"],
            cfg["offscreen_flag"],
            str(cfg["process_script"]),
            str(state_path),
            str(out_png),
            "--template", cfg["particle_filename_template"],
        ]
        if cfg.get("project_dir"):
            cmd += ["--project-dir", str(cfg["project_dir"])]
        cmd += [str(p) for p in vtk_paths]

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

        vtk_paths: list[Path] = item
        names = [p.name for p in vtk_paths]
        log.info("Processing: %s", names)
        time.sleep(cfg["settle_time"])

        missing = [p for p in vtk_paths if not p.exists()]
        if missing:
            log.warning("Files vanished before processing: %s", [p.name for p in missing])
            queue.task_done()
            continue

        success = render_vtk(vtk_paths, cfg)

        if success:
            for p in vtk_paths:
                p.unlink()
            log.info("Deleted: %s", names)
        else:
            log.error("Keeping %s for inspection", names)

        queue.task_done()


# ---------------------------------------------------------------------------
# Filesystem watcher
# ---------------------------------------------------------------------------

class VtkHandler(FileSystemEventHandler):
    def __init__(
        self,
        queue: Queue,
        patterns: list,
        template_regex: re.Pattern,
        expected_files: int,
        group_timeout: float,
    ):
        self.queue          = queue
        self.patterns       = patterns
        self.template_regex = template_regex
        self.expected_files = expected_files
        self.group_timeout  = group_timeout
        self._seen: set[str] = set()
        self._buffer: dict[str, dict[str, Path]] = {}
        self._buffer_times: dict[str, float] = {}
        self._lock = threading.Lock()

    def _matches(self, path: str) -> bool:
        name = Path(path).name
        return any(pat.search(name) for pat in self.patterns)

    def _extract_timestep(self, path: str) -> str | None:
        name = Path(path).name
        m = self.template_regex.search(name)
        if m and "time" in m.groupdict():
            return m.group("time")
        return None

    def _enqueue_group(self, timestep: str):
        group = self._buffer.pop(timestep, {})
        self._buffer_times.pop(timestep, None)
        if group:
            paths = list(group.values())
            self.queue.put(paths)
            log.info(
                "Queued group (t=%s): %s  (queue depth: %d)",
                timestep, [p.name for p in paths], self.queue.qsize(),
            )

    def _buffer_file(self, path: str):
        if not self._matches(path) or path in self._seen:
            return
        self._seen.add(path)

        if self.expected_files <= 1:
            self.queue.put([Path(path)])
            log.info("Queued: %s  (queue depth: %d)", Path(path).name, self.queue.qsize())
            return

        timestep = self._extract_timestep(path)
        if timestep is None:
            log.warning("Could not extract timestep from %s, enqueueing individually", path)
            self.queue.put([Path(path)])
            return

        with self._lock:
            if timestep not in self._buffer:
                self._buffer[timestep] = {}
                self._buffer_times[timestep] = time.time()
            self._buffer[timestep][path] = Path(path)
            count = len(self._buffer[timestep])
            log.info("Buffered (t=%s, %d/%d): %s", timestep, count, self.expected_files, Path(path).name)

            if count >= self.expected_files:
                self._enqueue_group(timestep)

    def flush_stale_groups(self):
        with self._lock:
            now = time.time()
            stale = [
                ts for ts, t0 in self._buffer_times.items()
                if now - t0 > self.group_timeout
            ]
            for ts in stale:
                log.warning(
                    "Flushing incomplete group (t=%s, %d/%d files) after %.0fs timeout",
                    ts, len(self._buffer.get(ts, {})), self.expected_files, self.group_timeout,
                )
                self._enqueue_group(ts)

    def on_created(self, event):
        if not event.is_directory:
            self._buffer_file(str(event.src_path))

    def on_moved(self, event):           # catches atomic rename from C++
        if not event.is_directory:
            self._buffer_file(str(event.dest_path))


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
    log.info("Template     : %s", cfg['particle_filename_template'])
    log.info("Expected/timestep: %d", cfg['expected_files_per_timestep'])
    log.info("Group timeout: %.0fs", cfg['group_timeout'])
    log.info("Project dir  : %s", cfg.get('project_dir') or '(none)')
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
    tmpl_regex = cfg["_tmpl_regex"]
    expected   = cfg["expected_files_per_timestep"]

    if expected <= 1:
        for existing in sorted(results_dir.iterdir()):
            if not existing.is_file():
                continue
            name = existing.name
            if any(pat.search(name) for pat in cfg["patterns"]):
                queue.put([existing])
                log.info("Queued pre-existing: %s", name)
    else:
        groups: dict[str, list[Path]] = {}
        for existing in sorted(results_dir.iterdir()):
            if not existing.is_file():
                continue
            name = existing.name
            if not any(pat.search(name) for pat in cfg["patterns"]):
                continue
            m = tmpl_regex.search(name)
            timestep = m.group("time") if (m and "time" in m.groupdict()) else None
            if timestep is None:
                queue.put([existing])
                log.info("Queued pre-existing (no timestep): %s", name)
                continue
            groups.setdefault(timestep, []).append(existing)

        for timestep, paths in sorted(groups.items()):
            if len(paths) >= expected:
                queue.put(paths)
                log.info("Queued pre-existing group (t=%s): %s", timestep, [p.name for p in paths])
            else:
                log.warning(
                    "Incomplete pre-existing group (t=%s, %d/%d): %s — enqueuing anyway",
                    timestep, len(paths), expected, [p.name for p in paths],
                )
                queue.put(paths)

    handler = VtkHandler(
        queue,
        cfg["patterns"],
        tmpl_regex,
        expected,
        cfg["group_timeout"],
    )
    observer = Observer()
    observer.schedule(handler, str(results_dir), recursive=False)
    observer.start()
    log.info("Watching for new files — press Ctrl+C to stop.")

    # Flush thread: periodically flush incomplete groups past their timeout
    def _flush_loop():
        while True:
            time.sleep(10)
            handler.flush_stale_groups()

    flush_thread = Thread(target=_flush_loop, daemon=True)
    flush_thread.start()

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
