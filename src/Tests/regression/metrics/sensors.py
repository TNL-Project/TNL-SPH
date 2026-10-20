"""Sensor-file comparator and metric functions of the regression tier.

Comparator convention (ONE convention, used consistently here, in the seeded
reference JSONs, and by :data:`regression.metrics.REGISTRY` functions):

    For ``results/<sensor_file>`` the comparator consumes the reference
    metrics named ``<sensor-stem>_avg`` and ``<sensor-stem>_std`` (stem =
    filename without suffix), comparing the per-column mean and standard
    deviation of the sensor data (first column = time, dropped) against
    ``{"value", "atol", "rtol"}`` entries via
    ``numpy.allclose(sim, value, atol=atol, rtol=rtol)`` with ``rtol``
    defaulting to ``0.0``.  Each of the two entries is one sub-check, so a
    sensor file contributes ``(passed, 2)``.

``damBreak2D_WCSPH_DBC`` keeps semantic parity with the legacy
``evaluateExamplesMetrics.damBreak2D_WCSPH_DBC``: 4 sub-checks
(``sensorsWaterLevel_avg/_std``, ``sensorsPressure_avg/_std``) with the
tolerances ported verbatim (atol 0.1 / 0.01 / 100 / 100; legacy numpy
``allclose`` defaults ``rtol=1e-05`` whose contribution is orders of
magnitude below the stated atols, and the pinned comparator default is
``rtol=0.0`` -- the seeded JSON therefore stores ``"rtol": 0.0``).

Collect mode (todo 10's ``--generate-references`` contract): registry
functions accept ``collect=True`` and return ``{metric_name: value}`` with
JSON-serializable values instead of comparing.  The metric names are the
CANONICAL ones the comparator consumes (``<sensor-stem>_avg`` /
``<sensor-stem>_std``), so a generated staging reference is
comparator-consumable once frozen.  (The legacy short aliases
``wl_*``/``ps_*`` are gone -- they produced staging refs the comparator
could never read.)  Stubs return ``(0, 0)`` / ``{}``.

Reference selection rule (per-tag resolution): a registry function serving
several case tags loads the PER-TAG committed reference
``references.reference_path( case_tag )`` when that file exists and falls
back to its hardcoded BASE reference otherwise -- entries frozen against
the base keep the pre-variant default behavior.

Diagnostic lines are printed per sub-check in the legacy per-check style,
naming the metric -- including on failure (selfcheck relies on it).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import references

# The base reference of the 2D damBreak evaluation (seeded from the legacy
# hardcoded arrays of evaluateExamplesMetrics.damBreak2D_WCSPH_DBC); the
# fallback for every served case tag without its own committed reference.
_DAMBREAK2D_WCSPH_DBC_TAG = "WCSPH-DBC/damBreak2D_WCSPH-DBC"

# Sensor filenames of the 2D damBreak evaluation; the canonical metric
# names are ``<filename-stem>_avg`` / ``<filename-stem>_std``.
_DAMBREAK2D_SENSOR_FILES = (
    "sensorsWaterLevel.dat",
    "sensorsPressure.dat",
)


def _sensor_mean_std(
    case_dir: Path, sensor_file: str, *, rows_drop_first_col: bool = True
) -> tuple[ np.ndarray, np.ndarray ]:
    """Per-column mean and std of ``case_dir/results/<sensor_file>``."""
    sensor_path = case_dir / "results" / sensor_file
    if not sensor_path.is_file():
        raise FileNotFoundError(
            f"sensor_avg_std_check: sensor file not found: {sensor_path}"
        )
    data = np.atleast_2d( np.genfromtxt( sensor_path, delimiter=" " ) )
    if rows_drop_first_col:
        data = data[ :, 1: ]
    return np.mean( data, axis=0 ), np.std( data, axis=0 )


def sensor_avg_std_check(
    case_dir: Path,
    sensor_file: str,
    ref: dict,
    *,
    rows_drop_first_col: bool = True,
) -> tuple[ int, int ]:
    """Compare per-column mean+std of ``results/<sensor_file>`` against ``ref``.

    Consumes the ``<sensor-stem>_avg`` / ``<sensor-stem>_std`` entries of
    ``ref["metrics"]``; each entry is one sub-check.  Prints one legacy-style
    diagnostic line per sub-check naming the metric.  Returns
    ``(passed, total)`` with ``total == 2``.
    """
    stem = Path( sensor_file ).stem
    label = ref.get( "case-tag" ) if isinstance( ref, dict ) else None
    if not label:
        label = Path( case_dir ).name
    metrics = ref.get( "metrics", {} )

    mean_sim, std_sim = _sensor_mean_std(
        case_dir, sensor_file, rows_drop_first_col=rows_drop_first_col
    )

    passed = 0
    total = 0
    for suffix, sim in ( ( "avg", mean_sim ), ( "std", std_sim ) ):
        metric_name = f"{stem}_{suffix}"
        total += 1
        entry = metrics.get( metric_name )
        if entry is None:
            raise ValueError(
                f"sensor_avg_std_check: metric '{metric_name}' for sensor"
                f" '{sensor_file}' not found in reference"
                f" (available: {sorted( metrics )})"
            )
        ok = np.allclose(
            sim,
            entry[ "value" ],
            atol=entry[ "atol" ],
            rtol=entry.get( "rtol", 0.0 ),
        )
        if ok:
            passed += 1
            print( f"{label}: {metric_name} comparison - successful." )
        else:
            max_diff = float(
                np.max( np.abs( np.asarray( sim ) - np.asarray( entry[ "value" ] ) ) )
            )
            print( f"{label}: {metric_name} comparison - failed." )
            print(
                f"    sim={np.asarray( sim ).tolist()} ref={entry[ 'value' ]}"
                f" atol={entry[ 'atol' ]} rtol={entry.get( 'rtol', 0.0 )}"
                f" max|diff|={max_diff}"
            )
    return passed, total


def damBreak2D_WCSPH_DBC(
    case_dir: Path, *, collect: bool = False, case_tag: str | None = None
):
    """Evaluate the 2D DBC damBreak case against its reference JSON.

    Reference selection rule: ``tag = case_tag or <base tag>``; the per-tag
    committed reference ``references.reference_path( tag )`` is loaded when
    that file exists, otherwise the BASE reference
    ``references/WCSPH-DBC/damBreak2D_WCSPH-DBC.json`` (seeded from the
    legacy arrays) -- every entry frozen against the base keeps the
    pre-variant default behavior.

    Collect mode ignores references entirely and returns the CANONICAL
    metric names (``sensorsWaterLevel_avg/_std``,
    ``sensorsPressure_avg/_std``) -- the same names the comparator looks
    up -- so generated staging references are freeze-ready.
    """
    if collect:
        values = {}
        for sensor_file in _DAMBREAK2D_SENSOR_FILES:
            stem = Path( sensor_file ).stem
            mean_sim, std_sim = _sensor_mean_std( case_dir, sensor_file )
            values[ f"{stem}_avg" ] = mean_sim.tolist()
            values[ f"{stem}_std" ] = std_sim.tolist()
        return values
    path = references.reference_path( case_tag or _DAMBREAK2D_WCSPH_DBC_TAG )
    if not path.is_file():
        path = references.reference_path( _DAMBREAK2D_WCSPH_DBC_TAG )
    ref = references.load_reference( path )
    passed = 0
    total = 0
    for sensor_file in _DAMBREAK2D_SENSOR_FILES:
        sub_passed, sub_total = sensor_avg_std_check( case_dir, sensor_file, ref )
        passed += sub_passed
        total += sub_total
    return passed, total


def damBreak2D_WCSPH_MDBC( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def damBreak3D_WCSPH_DBC( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def damBreak3D_WCSPH_MDBC( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def damBreak2D_WCSPH_BI( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def damBreak3D_WCSPH_BI( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def poiseuilleFlowWithOpenBoundary2D_WCSPH( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def poiseuilleFlowWithPeriodicBoundary2D_WCSPH( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def poiseuilleFlowWithOpenBoundary2D_WCSPH_BI( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def poiseuilleFlowWithPeriodicBoundary2D_WCSPH_BI( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Stub (legacy semantics): no evaluation until references are frozen."""
    if collect:
        return {}
    return 0, 0


def empty( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ):
    """Default no-op evaluation (legacy ``evaluateExamplesMetrics.empty``)."""
    if collect:
        return {}
    return 0, 0
