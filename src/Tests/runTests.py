#!/usr/bin/env python3
"""Unified TNL-SPH test-suite orchestrator (skeleton: tier selection only, no runners yet)."""

import argparse
import sys
from pathlib import Path

# Locate the suite package relative to this file (same file-relative pattern as suite.paths).
_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))

import suite.bootstrap  # noqa: F401  (side effect: guarantees the suite sys.path setup)


def main(argv=None) -> int:
    """Parse the requested tiers and report that no runner is implemented yet."""
    parser = argparse.ArgumentParser(description="TNL-SPH unified test suite orchestrator.")
    parser.add_argument(
        "--tier",
        action="append",
        choices=("all", "unit", "infrastructure", "regression"),
        default=None,
        help="Test tier(s) to run; may be given multiple times (default: all).",
    )
    args = parser.parse_args(argv)

    tiers = args.tier or ["all"]
    print(f"runTests.py: tiers {tiers} not implemented yet.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
