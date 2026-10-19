"""Canonical filesystem locations used across the TNL-SPH unified test suite."""

from pathlib import Path

# Repository root: src/Tests/suite/paths.py -> parents[3] is the repo root.
ROOT = Path(__file__).resolve().parents[3]

EXAMPLES_DIR = ROOT / "examples"
TOOLS_DIR = ROOT / "src" / "tools"
BUILD_DIR = ROOT / "build"
TESTS_DIR = ROOT / "src" / "Tests"

assert ROOT.joinpath("CMakeLists.txt").exists(), (
    f"suite.paths: failed to locate repository root (resolved {ROOT})"
)
