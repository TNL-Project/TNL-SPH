#!/usr/bin/env python3
"""Benchmark TNL-SPH examples across particle-system classes and resolutions.

For each particle-system class (variant), the script:
  1. Patches ``template/config.h`` in the example to select the class.
  2. Compiles the example's CUDA target (single core by default).
  3. For every resolution (``--dp`` value):
       a. Runs ``init.py --dp <res>`` to (re)generate ``sources/``.
       b. Optionally overrides ``final-time`` in ``sources/config.jsonc``.
       c. Runs the compiled binary and captures its stdout.
  4. Parses the ``Computation time`` block printed by the solver and records
     the integrate / interact / search / total timings plus particle counts.
  5. Stores per-run logs, a ``summary.json`` and ``results.csv`` under
     ``benchmark_results/<timestamp>/``, and prints a comparison table with
     speedup relative to the first (baseline) variant.

The script is example-agnostic: pass ``--example`` with a directory (or a
predefined key like ``damBreak2D`` / ``damBreak3D``) and the target/resolutions
are auto-detected from the example's ``CMakeLists.txt`` and ``init.py``.

Examples
--------
  # Default: 2D dam break, base + Warp variant, default resolutions
  ./benchmark_particles.py

  # 3D dam break, custom resolutions, quick runs (0.2s final time)
  ./benchmark_particles.py --example damBreak3D \\
      --resolutions 0.03,0.02,0.015 --end-time 0.2

  # Only the Warp variant, one resolution
  ./benchmark_particles.py --variants Warp --resolutions 0.002
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

# --------------------------------------------------------------------------- #
# Variant registry
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Variant:
    """A particle-system class selectable via ``template/config.h``."""
    name: str
    include: str        # angle-bracketed header path
    using_expr: str     # RHS of ``using ParticlesSys = ...;``


VARIANTS: dict[str, Variant] = {
    "CLL": Variant(
        "CLL",
        "<TNL/Particles/ParticlesLinkedList.h>",
        "TNL::Particles::ParticlesLinkedList< ParticlesConfig, Device >",
    ),
    "CLLWithList": Variant(
        "CLLWithList",
        "<TNL/Particles/ParticlesLinkedListWithList.h>",
        "TNL::Particles::ParticlesLinkedListWithList< ParticlesConfig, Device >",
    ),
    # --- Experimental variants (build/_deps/tnl-src/.../experimental/) ---
    "Ellpack": Variant(
        "Ellpack",
        "<TNL/Particles/experimental/ParticlesLinkedListWithListEllpack.h>",
        "TNL::Particles::ParticlesLinkedListWithListEllpack< ParticlesConfig, Device >",
    ),
    "Warp": Variant(
        "Warp",
        "<TNL/Particles/experimental/ParticlesLinkedListWithListWarp.h>",
        "TNL::Particles::ParticlesLinkedListWithListWarp< ParticlesConfig, Device >",
    ),
    "SharedMem": Variant(
        "SharedMem",
        "<TNL/Particles/experimental/ParticlesLinkedListWithListSharedMem.h>",
        "TNL::Particles::ParticlesLinkedListWithListSharedMem< ParticlesConfig, Device >",
    ),
    "Optimized": Variant(
        "Optimized",
        "<TNL/Particles/experimental/ParticlesLinkedListWithListOptimized.h>",
        "TNL::Particles::ParticlesLinkedListWithListOptimized< ParticlesConfig, Device >",
    ),
    "Atomic": Variant(
        "Atomic",
        "<TNL/Particles/experimental/ParticlesLinkedListWithListAtomic.h>",
        "TNL::Particles::ParticlesLinkedListWithListAtomic< ParticlesConfig, Device >",
    ),
}

# --------------------------------------------------------------------------- #
# Predefined examples (auto-detection fills in target / default resolutions)
# --------------------------------------------------------------------------- #

EXAMPLES: dict[str, dict] = {
    "damBreak2D": {
        "dir": "examples/WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "resolutions": [0.0025, 0.002, 0.0015],
    },
    "damBreak3D": {
        "dir": "examples/WCSPH-DBC/damBreak3D_WCSPH-DBC",
        "resolutions": [0.03, 0.02, 0.015],
    },
}

# --------------------------------------------------------------------------- #
# Result dataclasses
# --------------------------------------------------------------------------- #

@dataclass
class RunResult:
    variant: str
    resolution: float
    particle_system_type: str = ""
    fluid_particles: int = 0
    boundary_particles: int = 0
    integrate: Optional[float] = None
    interact: Optional[float] = None
    search: Optional[float] = None
    total_time: Optional[float] = None
    simulation_steps: int = 0
    final_sim_time: Optional[float] = None
    build_seconds: Optional[float] = None
    run_seconds: Optional[float] = None
    ok: bool = True
    error: str = ""
    log_path: str = ""


# --------------------------------------------------------------------------- #
# config.h patching
# --------------------------------------------------------------------------- #

# Matches the active (uncommented) particle-system include + using pair.
# Commented alternatives (``//#include`` / ``//using``) are skipped because
# the line must start with ``#include``.
_ACTIVE_PAIR = re.compile(
    r'^(#include\s+<TNL/Particles[^>]+>\s*\nusing\s+ParticlesSys\s*=\s*[^;\n]+;)',
    re.MULTILINE,
)


def patch_config_h(config_h: Path, variant: Variant) -> bool:
    """Swap the active particle-system class in ``config.h``.

    Returns True if a pair was found and replaced.
    """
    text = config_h.read_text()
    replacement = f"#include {variant.include}\nusing ParticlesSys = {variant.using_expr};"
    new_text, n = _ACTIVE_PAIR.subn(replacement, text, count=1)
    if n == 0:
        return False
    config_h.write_text(new_text)
    return True


def detect_active_variant(config_h: Path) -> Optional[str]:
    """Return the name of the currently active variant in config.h, if known."""
    text = config_h.read_text()
    m = _ACTIVE_PAIR.search(text)
    if not m:
        return None
    for name, v in VARIANTS.items():
        if v.include in m.group(1) and v.using_expr in m.group(1):
            return name
    return None


# --------------------------------------------------------------------------- #
# Example auto-detection
# --------------------------------------------------------------------------- #

def resolve_example_dir(arg: str) -> Path:
    """Resolve a --example argument to an example directory."""
    if arg in EXAMPLES:
        return Path(EXAMPLES[arg]["dir"])
    p = Path(arg)
    if not p.is_absolute():
        p = Path.cwd() / p
    if not p.is_dir():
        sys.exit(f"error: example directory not found: {arg}")
    return p


def detect_target(example_dir: Path) -> str:
    """Read the first add_executable target from the example's CMakeLists.txt."""
    cmake = example_dir / "CMakeLists.txt"
    if not cmake.is_file():
        sys.exit(f"error: no CMakeLists.txt in {example_dir}")
    m = re.search(r'add_executable\(\s*(\S+)\s+', cmake.read_text())
    if not m:
        sys.exit(f"error: could not detect build target in {cmake}")
    return m.group(1)


def detect_init_dp_default(example_dir: Path) -> float:
    """Return the default --dp value declared in init.py (for informative use)."""
    init = example_dir / "init.py"
    if not init.is_file():
        return 0.0
    m = re.search(r'add_argument\(\s*["\']--dp["\'].*?default\s*=\s*([0-9.eE+-]+)',
                 init.read_text(), re.DOTALL)
    return float(m.group(1)) if m else 0.0


# --------------------------------------------------------------------------- #
# Build / init / run
# --------------------------------------------------------------------------- #

def run_cmd(cmd: list, cwd: Path, log_file: Optional[Path] = None,
            capture: bool = True, check: bool = True) -> str:
    """Run a command, streaming + capturing stdout. Returns combined output."""
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.Popen(
        cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        lines.append(line)
        if not capture:
            sys.stdout.write(line)
    proc.wait()
    output = "".join(lines)
    if log_file:
        log_file.write_text(output)
    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output)
    return output


def build_target(project_dir: Path, build_dir: Path, target: str, jobs: int) -> float:
    """Compile the given target. Returns elapsed seconds."""
    t0 = time.perf_counter()
    run_cmd(["cmake", "--build", str(build_dir), "--target", target, "-j", str(jobs)],
            cwd=project_dir, capture=False, check=True)
    return time.perf_counter() - t0


def run_init(example_dir: Path, dp: float, extra_args: list[str]) -> None:
    """Run init.py to regenerate sources/ for a given resolution."""
    run_cmd([sys.executable, "init.py", "--dp", str(dp), *extra_args],
            cwd=example_dir, capture=False, check=True)


def patch_final_time(config_jsonc: Path, end_time: float) -> None:
    """Override final-time in the generated config.jsonc (quick-benchmark mode)."""
    if not config_jsonc.is_file():
        return
    text = config_jsonc.read_text()
    text = re.sub(
        r'("final-time"\s*:\s*)[0-9.eE+-]+',
        rf'\g<1>{end_time}',
        text,
    )
    config_jsonc.write_text(text)


def binary_path(project_dir: Path, build_dir: Path, example_dir: Path, target: str) -> Path:
    """Locate the compiled binary mirroring the source tree layout."""
    rel = example_dir.resolve().relative_to(project_dir.resolve())
    return build_dir / rel / target


# --------------------------------------------------------------------------- #
# Log parsing
# --------------------------------------------------------------------------- #

def _field(text: str, key: str) -> Optional[str]:
    """Extract the value of ``| key: value |`` or ``| key  value |`` rows."""
    m = re.search(rf'\|\s*{re.escape(key)}\s*:?\s*([^|\n]+?)\s*\|', text)
    return m.group(1).strip() if m else None


def _timing(text: str, key: str) -> Optional[float]:
    """Extract a timing value from the ``Computation time`` block.

    Keys like ``integrate``/``search`` must not match ``integrate-average``
    or ``search_cellIndices`` — enforced by requiring whitespace immediately
    after the key.
    """
    m = re.search(rf'^\|\s*{re.escape(key)}\s+([0-9.eE+-]+)\s*\|', text, re.MULTILINE)
    return float(m.group(1)) if m else None


def parse_run_log(text: str) -> dict:
    """Parse solver stdout for configuration + timing fields."""
    out: dict = {}
    out["particle_system_type"] = _field(text, "Particle system type") or ""
    dp = _field(text, "Initial particle distance (dp)")
    out["dp_parsed"] = float(dp) if dp else None
    fp = _field(text, "Number of allocated fluid particles")
    out["fluid_particles"] = int(fp.replace(",", "")) if fp else 0
    bp = _field(text, "Number of allocated boundary particles")
    out["boundary_particles"] = int(bp.replace(",", "")) if bp else 0
    out["integrate"] = _timing(text, "integrate")
    out["interact"] = _timing(text, "interact")
    out["search"] = _timing(text, "search")
    out["total_time"] = _timing(text, "Total time:")
    # Last "Simulation time: X s, simulation step: N" line → final values
    steps = re.findall(r'simulation step:\s*(\d+)', text)
    out["simulation_steps"] = int(steps[-1]) if steps else 0
    simtimes = re.findall(r'Simulation time:\s*([0-9.eE+-]+)\s*s', text)
    out["final_sim_time"] = float(simtimes[-1]) if simtimes else None
    return out


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #

def _fmt(v: Optional[float], spec: str = ".3f") -> str:
    return format(v, spec) if v is not None else "—"


def print_table(results: list[RunResult], baseline: str) -> None:
    """Print a comparison table; rich if available, else plain text."""
    if not results:
        print("(no results)")
        return
    base_total = next((r.total_time for r in results
                       if r.variant == baseline and r.total_time), None)

    headers = ["Variant", "dp", "Particles", "Total [s]",
               "Search [s]", "Interact [s]", "Integrate [s]", "Steps",
               "Speedup", "Status"]
    rows = []
    for r in results:
        speedup = (base_total / r.total_time) if (base_total and r.total_time) else None
        rows.append([
            r.variant,
            f"{r.resolution:g}",
            f"{r.fluid_particles:,}" if r.fluid_particles else "—",
            _fmt(r.total_time),
            _fmt(r.search),
            _fmt(r.interact),
            _fmt(r.integrate),
            str(r.simulation_steps) if r.simulation_steps else "—",
            f"{speedup:.2f}x" if speedup else "—",
            "ok" if r.ok else f"FAIL: {r.error[:40]}",
        ])

    try:
        from rich.console import Console
        from rich.table import Table
        table = Table(title="TNL-SPH particle-class benchmark", show_lines=True)
        for h in headers:
            table.add_column(h)
        for row in rows:
            table.add_row(*row)
        Console().print(table)
        return
    except ImportError:
        pass

    # Plain-text fallback
    widths = [max(len(str(h)), *(len(str(row[i])) for row in rows))
              for i, h in enumerate(headers)]
    sep = "+".join("-" * (w + 2) for w in widths)
    print("TNL-SPH particle-class benchmark")
    print(sep)
    print(" | ".join(h.ljust(w) for h, w in zip(headers, widths)))
    print(sep.replace("-", "="))
    for row in rows:
        print(" | ".join(str(c).ljust(w) for c, w in zip(row, widths)))
    print(sep)


def write_summary(out_dir: Path, results: list[RunResult]) -> None:
    """Persist results as JSON + CSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(
        json.dumps([asdict(r) for r in results], indent=2, default=str))
    with (out_dir / "results.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(asdict(results[0]).keys())
        for r in results:
            w.writerow(asdict(r).values())


# --------------------------------------------------------------------------- #
# Main orchestration
# --------------------------------------------------------------------------- #

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Benchmark TNL-SPH examples across particle classes and resolutions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--example", default="damBreak2D",
                    help="example directory or predefined key "
                         f"({', '.join(EXAMPLES)}). Default: damBreak2D")
    ap.add_argument("--variants", default="CLLWithList,Warp",
                    help="comma-separated variant names from: "
                         f"{', '.join(VARIANTS)}. Default: CLLWithList,Warp")
    ap.add_argument("--resolutions", default=None,
                    help="comma-separated dp values (overrides example defaults)")
    ap.add_argument("--target", default=None,
                    help="build target (auto-detected from CMakeLists.txt if omitted)")
    ap.add_argument("--build-dir", default=None,
                    help="CMake build directory (default: <repo>/build)")
    ap.add_argument("--jobs", type=int, default=1,
                    help="compile parallelism (default: 1 = single core)")
    ap.add_argument("--end-time", type=float, default=None,
                    help="override final-time in config.jsonc for quick runs "
                         "(default: run the full simulation)")
    ap.add_argument("--init-args", default="",
                    help="extra args passed to init.py, e.g. '--h-coef 2.0'")
    ap.add_argument("--results-dir", default="benchmark_results",
                    help="directory for stored results (default: benchmark_results)")
    ap.add_argument("--rebuild", action="store_true",
                    help="rebuild before every variant (default: build once per variant)")
    ap.add_argument("--no-restore", action="store_true",
                    help="do not restore the original config.h on exit")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the plan without building/running")
    args = ap.parse_args()

    # Resolve variant list
    variant_names = [v.strip() for v in args.variants.split(",") if v.strip()]
    for name in variant_names:
        if name not in VARIANTS:
            sys.exit(f"error: unknown variant '{name}'. Known: {', '.join(VARIANTS)}")
    variants = [VARIANTS[name] for name in variant_names]
    baseline = variant_names[0]

    # Resolve example
    example_dir = resolve_example_dir(args.example)
    project_dir = Path(__file__).resolve().parent
    build_dir = Path(args.build_dir) if args.build_dir else (project_dir / "build")
    target = args.target or detect_target(example_dir)

    # Resolve resolutions
    if args.resolutions:
        resolutions = [float(x) for x in args.resolutions.split(",") if x.strip()]
    elif args.example in EXAMPLES:
        resolutions = EXAMPLES[args.example]["resolutions"]
    else:
        default_dp = detect_init_dp_default(example_dir)
        resolutions = [default_dp] if default_dp else [0.002]

    config_h = example_dir / "template" / "config.h"
    if not config_h.is_file():
        sys.exit(f"error: {config_h} not found")

    init_extra = args.init_args.split() if args.init_args else []

    # Plan
    print(f"Example      : {example_dir}")
    print(f"Target       : {target}")
    print(f"Build dir    : {build_dir}")
    print(f"Variants     : {variant_names}  (baseline = {baseline})")
    print(f"Resolutions  : {resolutions}")
    print(f"Jobs         : {args.jobs}")
    print(f"End time     : {args.end_time if args.end_time else '(full simulation)'}")
    print()

    if args.dry_run:
        print("Dry run — no build/run performed.")
        return 0

    # Verify build dir
    if not (build_dir / "CMakeCache.txt").is_file():
        sys.exit(f"error: not a CMake build directory: {build_dir} "
                 "(run `cmake -B build -S .` first)")

    # Backup config.h
    backup = config_h.read_text()
    print(f"Backed up {config_h}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.results_dir) / timestamp

    results: list[RunResult] = []

    try:
        for variant in variants:
            print(f"\n{'='*70}\nVariant: {variant.name}\n{'='*70}")

            # 1. Patch config.h
            if not patch_config_h(config_h, variant):
                print(f"warning: could not patch active particle-system pair in {config_h}")
            active = detect_active_variant(config_h)
            if active != variant.name:
                print(f"warning: config.h active variant is {active}, expected {variant.name}")
            else:
                print(f"patched config.h -> {variant.name}")

            # 2. Build
            try:
                build_s = build_target(project_dir, build_dir, target, args.jobs)
                print(f"built {target} in {build_s:.1f}s")
            except subprocess.CalledProcessError as e:
                print(f"BUILD FAILED for {variant.name}: {e}")
                for res in resolutions:
                    results.append(RunResult(
                        variant=variant.name, resolution=res, ok=False,
                        error="build failed"))
                continue

            binp = binary_path(project_dir, build_dir, example_dir, target)
            if not binp.is_file():
                print(f"error: binary not found at {binp}")
                for res in resolutions:
                    results.append(RunResult(
                        variant=variant.name, resolution=res, ok=False,
                        error="binary missing"))
                continue

            # 3. For each resolution: init, (patch end-time), run, parse
            for res in resolutions:
                print(f"\n--- {variant.name} @ dp={res:g} ---")
                log_file = out_dir / f"{variant.name}__dp{res:g}.log"
                r = RunResult(variant=variant.name, resolution=res, log_path=str(log_file))

                try:
                    run_init(example_dir, res, init_extra)
                    if args.end_time is not None:
                        patch_final_time(example_dir / "sources" / "config.jsonc",
                                         args.end_time)

                    # Re-patch config.h in case init.py regenerated it from a
                    # template (3D case). The binary is already compiled, so
                    # this only keeps config.h consistent for inspection.
                    patch_config_h(config_h, variant)

                    t0 = time.perf_counter()
                    output = run_cmd(
                        [str(binp), "--config", "sources/config.jsonc"],
                        cwd=example_dir, log_file=log_file, capture=False, check=True)
                    r.run_seconds = time.perf_counter() - t0
                    r.build_seconds = build_s

                    parsed = parse_run_log(output)
                    r.particle_system_type = parsed["particle_system_type"]
                    r.fluid_particles = parsed["fluid_particles"]
                    r.boundary_particles = parsed["boundary_particles"]
                    r.integrate = parsed["integrate"]
                    r.interact = parsed["interact"]
                    r.search = parsed["search"]
                    r.total_time = parsed["total_time"]
                    r.simulation_steps = parsed["simulation_steps"]
                    r.final_sim_time = parsed["final_sim_time"]
                    print(f"total={_fmt(r.total_time)}  search={_fmt(r.search)}  "
                          f"interact={_fmt(r.interact)}  particles={r.fluid_particles}")
                except subprocess.CalledProcessError as e:
                    r.ok = False
                    r.error = f"command failed (rc={e.returncode})"
                    print(f"RUN FAILED: {r.error}  (see {log_file})")
                except Exception as e:  # noqa: BLE001
                    r.ok = False
                    r.error = f"{type(e).__name__}: {e}"
                    print(f"ERROR: {r.error}")

                results.append(r)
                write_summary(out_dir, results)  # incremental save

    finally:
        if not args.no_restore:
            config_h.write_text(backup)
            print(f"\nRestored original {config_h}")

    # Final report
    print(f"\nResults stored in: {out_dir}")
    print_table(results, baseline)
    return 0 if all(r.ok for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
