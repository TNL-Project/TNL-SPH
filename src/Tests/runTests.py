#!/usr/bin/env python3
"""Unified TNL-SPH test-suite orchestrator (thin entry; all logic in suite.cli)."""

import sys
from pathlib import Path

# Locate the suite package relative to this file (same file-relative pattern as suite.paths).
_TESTS_DIR = Path(__file__).resolve().parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))

import suite.bootstrap  # noqa: F401  (side effect: guarantees the suite sys.path setup)
from suite.cli import main

if __name__ == "__main__":
    sys.exit(main())
