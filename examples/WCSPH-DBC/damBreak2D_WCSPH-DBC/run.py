#! /usr/bin/env python3

import argparse
import configparser
import subprocess
from pathlib import Path

# initialize directories
example_dir = Path(__file__).parent
project_dir = (example_dir / ".." / ".." / ".." ).resolve()
bin_dir = project_dir / "build" / example_dir.relative_to(project_dir)

def solve():
    solver_path = bin_dir / f"damBreak2D_WCSPH-DBC_cuda"

    args = []
    args += [ solver_path ]

    # run the process and print its output as it is being executed
    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          bufsize=1, cwd=example_dir, text=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)

if __name__ == "__main__":

    solve()
