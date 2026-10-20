"""
compare.py

Deep structural comparison ported verbatim (semantics) from
``src/Tests/infrastructure/run_sph_tests.py:110-162`` plus a loader for
frozen init-report references.

Behavioral parity with the original engine is asserted by
``infrastructure.validators.tests`` against the legacy function itself.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional


# ── Deep comparison (port of run_sph_tests.deep_compare) ─────────────────────

class Diff:
    def __init__(self, path: str, expected: Any, actual: Any):
        self.path = path
        self.expected = expected
        self.actual = actual

    def __str__(self):
        return (
            f"  {self.path}\n"
            f"    expected : {self.expected!r}\n"
            f"    actual   : {self.actual!r}"
        )


def deep_compare(expected: Any, actual: Any, path: str = '', tol: float = 1e-9) -> list[Diff]:
    """Recursively compare two structures. Returns list of Diff objects."""
    diffs = []

    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            return [Diff(path, expected, actual)]
        all_keys = set(expected) | set(actual)
        for k in sorted(all_keys):
            child = f'{path}.{k}' if path else k
            if k not in expected:
                continue
            elif k not in actual:
                diffs.append(Diff(child, expected[k], '<missing>'))
            else:
                diffs.extend(deep_compare(expected[k], actual[k], child, tol))

    elif isinstance(expected, list):
        if not isinstance(actual, list):
            return [Diff(path, expected, actual)]
        if len(expected) != len(actual):
            diffs.append(Diff(path + '.length', len(expected), len(actual)))
        for i, (e, a) in enumerate(zip(expected, actual)):
            diffs.extend(deep_compare(e, a, f'{path}[{i}]', tol))

    elif isinstance(expected, float) or isinstance(actual, float):
        try:
            e, a = float(expected), float(actual)
            if abs(e - a) > tol * max(1.0, abs(e)):
                diffs.append(Diff(path, expected, actual))
        except (TypeError, ValueError):
            if expected != actual:
                diffs.append(Diff(path, expected, actual))

    else:
        if expected != actual:
            diffs.append(Diff(path, expected, actual))

    return diffs


# ── Frozen reference store ───────────────────────────────────────────────────

_REFERENCES_DIR = Path(__file__).resolve().parents[1] / "references"


def sanitize_tag(case_tag: str) -> str:
    """Same convention as regression.references: ':' -> '__', '/' kept."""
    return case_tag.replace(":", "__")


def load_frozen_reference(case_tag: str) -> Optional[dict]:
    """
    Load a frozen init-report reference for *case_tag* from
    ``infrastructure/references/<sanitized>.json``.

    Missing reference is skip-friendly: returns None, never raises.
    """
    ref_path = _REFERENCES_DIR / f"{sanitize_tag(case_tag)}.json"
    if not ref_path.exists():
        return None
    return json.loads(ref_path.read_text(encoding="utf-8"))
