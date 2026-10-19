"""Selftest entry point: ``PYTHONPATH=src/Tests python3 -m regression.metrics``.

Runs :func:`regression.metrics.selftest` (the optional ``--selftest`` token
is accepted for plan-acceptance compatibility and ignored).
"""

from __future__ import annotations

from . import selftest

if __name__ == "__main__":
    raise SystemExit( selftest() )
