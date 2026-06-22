#! /usr/bin/env python3

"""
CFL parametric study for damBreak2D_WCSPH-BI.

For each configuration x CFL combination:
  1. init.py     - regenerate particles + config files with the given CFL
  2. cmake build - rebuild the solver (time-integration is a compiled template param)
  3. solver      - run the simulation
  4. collect     - success/crash + computational time from time_measurements.json

Results are written incrementally to a CSV file after each run, so partial
results survive a mid-run crash.

Usage:
  ./cfl_study.py
  ./cfl_study.py --cfls 0.005 0.01 0.05 0.1 0.2 0.5 1.0 2.0
  ./cfl_study.py --cfl-range 0.005 2.0 --cfl-count 20
  ./cfl_study.py --cfl-range 0.005 2.0 --cfl-step 0.1
  ./cfl_study.py --cfl-range 0.005 2.0 --cfl-count 20 --log-spacing
  ./cfl_study.py --skip-build   # skip cmake build, reuse existing binary
  ./cfl_study.py --backup-results  # rename results/ after each run
"""

import csv
import json
import shutil
import argparse
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

script_dir = Path(__file__).parent.resolve()
project_dir = (script_dir / ".." / ".." / "..").resolve()
example_dir = project_dir / "examples" / "WCSPH-BI" / "damBreak2D_WCSPH-BI"
bin_dir = project_dir / "build" / "examples" / "WCSPH-BI" / "damBreak2D_WCSPH-BI"
results_dir = example_dir / "results"
build_dir = project_dir / "build"
target_name = "damBreak2D_WCSPH-BI_cuda"


def emptyTest(case_dir):
    return 0, 0


CFLTest = [
    {
        "case-tag": "damBreak2D_WCSPH-BI_Verlet",
        "case": "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type": "BIConservative_numeric",
        "viscous-term": "None",
        "diffusive-term": "None",
        "dp": 0.00075,
        "h-coef": 4,
        "time-integration": "Verlet",
        "evaluation-function": emptyTest,
    },
    {
        "case-tag": "damBreak2D_WCSPH-BI_Midpoint",
        "case": "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type": "BIConservative_numeric",
        "viscous-term": "None",
        "diffusive-term": "None",
        "dp": 0.00075,
        "h-coef": 4,
        "time-integration": "MidpointScheme",
        "evaluation-function": emptyTest,
    },
    {
        "case-tag": "damBreak2D_WCSPH-BI_MidpointAnderson",
        "case": "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type": "BIConservative_numeric",
        "viscous-term": "None",
        "diffusive-term": "None",
        "dp": 0.00075,
        "h-coef": 4,
        "time-integration": "MidpointSchemeWithAnderson",
        "evaluation-function": emptyTest,
    },
]

DEFAULT_CFLS = [0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]


def init_case(cfl, conf):
    args = [str(example_dir / "init.py")]
    for key, value in conf.items():
        if key in ("case", "case-tag", "evaluation-function"):
            continue
        args += [f"--{key}", str(value)]
    args += ["--cfl", str(cfl)]
    subprocess.run(
        args,
        check=True,
        cwd=example_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )


def build_case():
    p = subprocess.run(
        ["cmake", "--build", str(build_dir), "--target", target_name],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)


def run_case():
    solver_path = bin_dir / target_name
    config_path = example_dir / "sources" / "config.ini"
    p = subprocess.run(
        [str(solver_path), "--config", str(config_path)],
        cwd=example_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return p.returncode == 0


def parse_computational_time():
    filename = results_dir / "time_measurements.json"
    try:
        with open(filename) as f:
            data = json.load(f)
            return float(data["total"])
    except Exception:
        return None


def backup_results(case_tag, cfl):
    cfl_str = f"{cfl}".replace(".", "p")
    dest = example_dir / f"results_{case_tag}_CFL{cfl_str}"
    if dest.exists():
        shutil.rmtree(dest)
    if results_dir.exists():
        shutil.move(str(results_dir), str(dest))


def run_cfl_study(cfls, skip_build=False, csv_path=None, backup=False):
    """
    Run the full CFL study.

    Returns a dict:  { case_tag: { cfl: {"success": bool, "time": float|None} } }

    If csv_path is given, each result is appended to the CSV immediately
    with a flush, so partial results survive a mid-run crash.
    If backup is True, results/ is renamed after each run.
    """
    results = {}

    csv_file = None
    csv_writer = None
    if csv_path:
        csv_file = open(csv_path, "w", newline="")
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["case_tag", "cfl", "success", "computational_time"])

    try:
        for conf in CFLTest:
            case_tag = conf["case-tag"]
            results[case_tag] = {}
            print(f"\n{'=' * 60}")
            print(f"Configuration: {case_tag}")
            print(f"  time-integration: {conf.get('time-integration', 'default')}")
            print(f"{'=' * 60}")

            if not skip_build:
                print("  [build] init + cmake --build ...")
                init_case(cfls[0], conf)
                build_case()
                print("  [build] done.")
            else:
                print("  [build] skipped (--skip-build).")

            for cfl in cfls:
                print(f"  CFL={cfl} ... ", end="", flush=True)

                try:
                    init_case(cfl, conf)
                except subprocess.CalledProcessError:
                    results[case_tag][cfl] = {"success": False, "time": None}
                    if csv_writer:
                        assert csv_file is not None
                        csv_writer.writerow([case_tag, cfl, False, ""])
                        csv_file.flush()
                    if backup:
                        backup_results(case_tag, cfl)
                    print("INIT FAILED")
                    continue

                success = run_case()
                comp_time = parse_computational_time() if success else None

                results[case_tag][cfl] = {"success": success, "time": comp_time}

                if csv_writer:
                    assert csv_file is not None
                    csv_writer.writerow(
                        [case_tag, cfl, success, comp_time if comp_time is not None else ""]
                    )
                    csv_file.flush()

                if backup:
                    backup_results(case_tag, cfl)

                if success and comp_time is not None:
                    print(f"OK ({comp_time:.2f}s)")
                elif success:
                    print("OK (no time data)")
                else:
                    print("CRASHED")
    finally:
        if csv_file:
            csv_file.close()

    return results


def plot_results(results, output_path):
    fig, (ax_time, ax_status) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    for case_tag, cfl_data in results.items():
        xs = sorted(cfl_data.keys())
        times = [cfl_data[c]["time"] for c in xs]
        successes = [cfl_data[c]["success"] for c in xs]

        valid_x = [c for c, t in zip(xs, times) if t is not None]
        valid_t = [t for t in times if t is not None]
        if valid_x:
            ax_time.plot(
                valid_x, valid_t, "o-", label=case_tag, markersize=7, linewidth=1.5
            )

        for cfl, ok in zip(xs, successes):
            y = 1 if ok else 0
            marker = "o" if ok else "X"
            color = "tab:green" if ok else "tab:red"
            ax_status.scatter(cfl, y, marker=marker, color=color, s=80, zorder=3)

    labels_done = set()
    for case_tag, cfl_data in results.items():
        for cfl, info in cfl_data.items():
            ok = info["success"]
            y = 1 if ok else 0
            if (y, case_tag) not in labels_done:
                ax_status.text(
                    cfl, y + 0.03, case_tag.split("_")[-1],
                    fontsize=7, alpha=0.7,
                )
                labels_done.add((y, case_tag))

    ax_time.set_ylabel("Computational Time [s]")
    ax_time.set_title("Computational Time vs CFL Number")
    ax_time.legend(fontsize=8)
    ax_time.grid(True, linestyle="--", alpha=0.5)
    ax_time.set_xscale("log")
    ax_time.set_yscale("log")

    ax_status.set_xlabel("CFL Number")
    ax_status.set_ylabel("Status")
    ax_status.set_yticks([0, 1])
    ax_status.set_yticklabels(["Crashed", "Success"])
    ax_status.set_title("Simulation Success / Crash")
    ax_status.grid(True, linestyle="--", alpha=0.5)
    ax_status.set_xscale("log")
    ax_status.set_ylim(-0.3, 1.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"\nPlot saved to: {output_path}")


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="CFL parametric study for damBreak2D_WCSPH-BI"
    )
    argparser.add_argument(
        "--cfls",
        type=float,
        nargs="+",
        default=None,
        help="explicit list of CFL values to test",
    )
    argparser.add_argument(
        "--cfl-range",
        type=float,
        nargs=2,
        metavar=("START", "END"),
        default=None,
        help="generate CFL values in [START, END] interval "
        "(use with --cfl-count or --cfl-step)",
    )
    argparser.add_argument(
        "--cfl-count",
        type=int,
        default=20,
        help="number of CFL values when using --cfl-range (default: %(default)s)",
    )
    argparser.add_argument(
        "--cfl-step",
        type=float,
        default=None,
        help="spacing between CFL values when using --cfl-range "
        "(overrides --cfl-count)",
    )
    argparser.add_argument(
        "--log-spacing",
        action="store_true",
        help="use logarithmic spacing for --cfl-range "
        "(useful for wide intervals like 0.005-2)",
    )
    argparser.add_argument(
        "--skip-build",
        action="store_true",
        help="skip cmake build, reuse the existing binary",
    )
    argparser.add_argument(
        "--backup-results",
        action="store_true",
        help="rename results/ to results_{tag}_CFL{value} after each run",
    )
    args = argparser.parse_args()

    if args.cfl_range:
        start, end = args.cfl_range
        if args.cfl_step is not None:
            count = int(round((end - start) / args.cfl_step)) + 1
            cfls = list(np.linspace(start, end, count))
        else:
            if args.log_spacing:
                cfls = list(
                    np.logspace(
                        np.log10(start),
                        np.log10(end),
                        args.cfl_count,
                    )
                )
            else:
                cfls = list(np.linspace(start, end, args.cfl_count))
        print(f"Generated {len(cfls)} CFL values from {start} to {end}"
              f" ({'log' if args.log_spacing else 'linear'} spacing):")
        print(f"  {[round(c, 6) for c in cfls]}")
    elif args.cfls:
        cfls = args.cfls
    else:
        cfls = DEFAULT_CFLS

    cfls = sorted(cfls)

    if not args.backup_results:
        results_dir.mkdir(exist_ok=True)

    csv_path = example_dir / "cfl_study_results.csv"

    results = run_cfl_study(
        cfls,
        skip_build=args.skip_build,
        csv_path=csv_path,
        backup=args.backup_results,
    )
    print(f"\nRaw results saved to: {csv_path}")

    plot_path = example_dir / "cfl_study.png"
    plot_results(results, plot_path)
