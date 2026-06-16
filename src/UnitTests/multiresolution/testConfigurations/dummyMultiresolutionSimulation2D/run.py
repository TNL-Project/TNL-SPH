#! /usr/bin/env python3

import os
import sys
import argparse
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from configurations import CONFIGURATIONS as _ALL_CONFIGURATIONS
CONFIGURATIONS = {k: v for k, v in _ALL_CONFIGURATIONS.items() if v.get("dimension") == 2}

example_dir = Path(__file__).parent
project_dir = (example_dir / ".." / ".." / ".." / ".." / ".." ).resolve()
bin_dir = project_dir / "build" / example_dir.relative_to(project_dir)


def solve(config_path: Path):
    solver_path = bin_dir / "dummyMultiresolutionSimulation2D"

    args = [
        solver_path,
        "--config", config_path,
    ]
    print(solver_path)
    print(args)

    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          bufsize=1, cwd=example_dir, text=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)

def move_results(name: str):
    results_dir = example_dir / "results"
    for f in os.listdir(example_dir / "results"):
        f_path = results_dir / f
        if os.path.isfile(f_path):
            os.rename(f_path, results_dir / f"{name}" / f)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Multiresolution test example")
    argparser.add_argument("--config-name", type=str, default="dummy-center",
                           choices=list(CONFIGURATIONS.keys()),
                           help="named configuration to run")
    argparser.add_argument("--all", action="store_true",
                           help="run all configurations sequentially")

    args = argparser.parse_args()

    if args.all:
        for name in CONFIGURATIONS:
            config_path = example_dir / f"sources/{name}/dummyConfig2D.ini"
            print(f"\n{'='*60}")
            print(f"Running configuration: {name}")
            print(f"{'='*60}\n")
            solve(config_path)
            move_results(name)
    else:
        config_path = example_dir / f"sources/{args.config_name}/dummyConfig2D.ini"
        solve(config_path)
        move_results(args.config_name)
