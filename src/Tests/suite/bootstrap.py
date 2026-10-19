"""Single sys.path hook of the TNL-SPH unified test suite (runs ensure_paths() on import)."""

import sys

from .paths import TESTS_DIR, TOOLS_DIR


def ensure_paths() -> None:
    """Insert the test-suite and tools directories into sys.path if absent."""
    for path in (str(TESTS_DIR), str(TOOLS_DIR)):
        if path not in sys.path:
            sys.path.insert(0, path)


# Import-time side effect: importing suite.bootstrap guarantees the path setup.
ensure_paths()
