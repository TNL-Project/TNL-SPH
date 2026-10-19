"""Reference store of the TNL-SPH regression metrics.

Reference JSON files live under ``src/Tests/regression/references/`` and hold
the verified metric values a case is compared against.  The on-disk layout
mirrors the sanitized case tag (todo 10's ``--generate-references`` stages
unverified refs under ``references/_staging/``, which comparators MUST
ignore -- see :func:`staging_refs_only`)::

    references/WCSPH-DBC/damBreak2D_WCSPH-DBC.json           (committed, verified)
    references/_staging/WCSPH-DBC/damBreak2D_WCSPH-DBC__MGVT.json   (unverified)

Tag sanitization (D6): ``:`` becomes ``__`` and ``/`` is kept as path
separation, e.g. ``WCSPH-DBC/damBreak2D_WCSPH-DBC:MGVT`` maps to
``WCSPH-DBC/damBreak2D_WCSPH-DBC__MGVT.json``.

Schema of a reference file::

    {"case-tag": <str>,
     "source": <str>,            # provenance, e.g. "ported from evaluateExamplesMetrics@<sha>"
     "date": <str>,              # ISO date the reference was produced/verified
     "verified": <bool>,         # only verified refs may live outside _staging/
     "metrics": {<name>: {"value": <float | list[float]>,
                          "atol": <float>,
                          "rtol": <float>}, ...}}

Metric names follow ``<sensor-stem>_avg`` / ``<sensor-stem>_std`` (see
:mod:`regression.metrics.sensors`).

This module MUST NOT print anything.  Writes are guarded:
:func:`save_reference` refuses any target outside the ``references/`` tree.
"""

from __future__ import annotations

import json
from pathlib import Path

from suite import paths

# Root of the committed reference store (staging lives in _staging/ below it).
REFERENCES_DIR = paths.TESTS_DIR / "regression" / "references"

_STAGING_DIRNAME = "_staging"


def sanitize_tag( case_tag: str ) -> str:
    """Map a case tag to its on-disk relative path (``:`` -> ``__``, ``/`` kept)."""
    return case_tag.replace( ":", "__" )


def reference_path( case_tag: str, *, staging: bool = False ) -> Path:
    """Absolute path of the reference JSON for ``case_tag``.

    With ``staging=True`` the path inside ``references/_staging/`` is
    returned; otherwise the committed-reference path.
    """
    base = REFERENCES_DIR / _STAGING_DIRNAME if staging else REFERENCES_DIR
    return base / ( sanitize_tag( case_tag ) + ".json" )


def load_reference( path: Path ) -> dict:
    """Load and return a reference JSON dict; errors always name ``path``."""
    if not path.is_file():
        raise FileNotFoundError( f"load_reference: reference file not found: {path}" )
    try:
        with path.open( encoding="utf-8" ) as handle:
            return json.load( handle )
    except json.JSONDecodeError as exc:
        raise ValueError( f"load_reference: invalid JSON in {path}: {exc}" ) from exc


def save_reference( ref: dict, path: Path ) -> Path:
    """Write ``ref`` as JSON (2-space indent) and return ``path``.

    Parent directories are created ONLY when the target lies inside the
    committed ``references/`` tree (``references/_staging/`` included);
    anything else raises :class:`ValueError` naming the offending path.
    """
    resolved = path.resolve()
    root = REFERENCES_DIR.resolve()
    if not resolved.is_relative_to( root ):
        raise ValueError(
            f"save_reference: refusing to write outside {root}: {resolved}"
        )
    resolved.parent.mkdir( parents=True, exist_ok=True )
    with resolved.open( "w", encoding="utf-8" ) as handle:
        json.dump( ref, handle, indent=2 )
        handle.write( "\n" )
    return resolved


def staging_refs_only( case_tag: str ) -> bool:
    """True iff a reference exists under ``_staging/`` but no committed one.

    Todo 10's unfrozen-SKIP rule depends on this predicate: such a case must
    report SKIP ("unfrozen reference"), never silently validate against the
    staging file.
    """
    return (
        reference_path( case_tag, staging=True ).is_file()
        and not reference_path( case_tag ).is_file()
    )
