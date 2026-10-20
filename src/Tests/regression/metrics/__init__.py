"""Metrics registry of the TNL-SPH regression tier.

``REGISTRY`` maps a configuration's ``"evaluation-function"`` name
(todo 7's configuration modules carry these as strings 1:1 onto the legacy
``evaluateExamplesMetrics`` names, ``empty`` included as the default) to the
function performing the evaluation.

Registry contract (todo 6's runner imports this LAZILY inside ``evaluate()``
and depends only on the names; todo 10's generator depends on collect mode)::

    fn( case_dir: Path, *, collect: bool = False, case_tag: str | None = None ) -> tuple[ int, int ] | dict

- compare mode (default): evaluate the results of the case in ``case_dir``
  against the committed reference store (see :mod:`.references`) and return
  ``(passed, total)`` sub-check counts, printing one legacy-style diagnostic
  line per sub-check.  Stubs return ``(0, 0)`` ("no evaluation yet").
- ``case_tag`` is the configuration entry's case tag, threaded by the
  caller (``suite.runner_core.evaluate``); a function serving several tags
  resolves its PER-TAG committed reference with it and falls back to its
  hardcoded base reference when the per-tag file is absent.  ``None``
  selects the base reference; collect mode ignores it.
- collect mode (``collect=True``, used by todo 10's
  ``--generate-references`` staging): return ``{metric_name: value}`` with
  the CANONICAL metric names the comparator consumes
  (``<sensor-stem>_avg`` / ``<sensor-stem>_std``) and JSON-serializable
  values, never compare -- generated staging refs are freeze-ready.
  Stubs return ``{}``.

Selftest::

    PYTHONPATH=src/Tests python3 -m regression.metrics [--selftest]

prints the sorted REGISTRY (11 entries), runs the tag-sanitizer assertions,
seeds a TEMP reference from the legacy arrays into ``references/_staging/``
(cleaned up afterwards) and demonstrates the comparator on synthetic sensor
files: a pass case ``(4, 4)`` including the variant-tag base fallback, and
a beyond-atol perturbed case ``(3, 4)`` naming the failing metric.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

REGISTRY: dict[ str, Callable ] = {}

from . import sensors  # noqa: E402  (populates REGISTRY below)

REGISTRY.update(
    {
        "damBreak2D_WCSPH_DBC": sensors.damBreak2D_WCSPH_DBC,
        "damBreak2D_WCSPH_MDBC": sensors.damBreak2D_WCSPH_MDBC,
        "damBreak3D_WCSPH_DBC": sensors.damBreak3D_WCSPH_DBC,
        "damBreak3D_WCSPH_MDBC": sensors.damBreak3D_WCSPH_MDBC,
        "damBreak2D_WCSPH_BI": sensors.damBreak2D_WCSPH_BI,
        "damBreak3D_WCSPH_BI": sensors.damBreak3D_WCSPH_BI,
        "poiseuilleFlowWithOpenBoundary2D_WCSPH": sensors.poiseuilleFlowWithOpenBoundary2D_WCSPH,
        "poiseuilleFlowWithPeriodicBoundary2D_WCSPH": sensors.poiseuilleFlowWithPeriodicBoundary2D_WCSPH,
        "poiseuilleFlowWithOpenBoundary2D_WCSPH_BI": sensors.poiseuilleFlowWithOpenBoundary2D_WCSPH_BI,
        "poiseuilleFlowWithPeriodicBoundary2D_WCSPH_BI": sensors.poiseuilleFlowWithPeriodicBoundary2D_WCSPH_BI,
        "empty": sensors.empty,
    }
)

_EXPECTED_NAMES = frozenset(
    {
        "damBreak2D_WCSPH_DBC",
        "damBreak2D_WCSPH_MDBC",
        "damBreak3D_WCSPH_DBC",
        "damBreak3D_WCSPH_MDBC",
        "damBreak2D_WCSPH_BI",
        "damBreak3D_WCSPH_BI",
        "poiseuilleFlowWithOpenBoundary2D_WCSPH",
        "poiseuilleFlowWithPeriodicBoundary2D_WCSPH",
        "poiseuilleFlowWithOpenBoundary2D_WCSPH_BI",
        "poiseuilleFlowWithPeriodicBoundary2D_WCSPH_BI",
        "empty",
    }
)

# The one gated reference of the plan (seeded from
# evaluateExamplesMetrics.py:14-17 verbatim; see sensors.py).
_COMMITTED_TAG = "WCSPH-DBC/damBreak2D_WCSPH-DBC"

# Legacy reference arrays, evaluateExamplesMetrics.py:13-16 verbatim.
_LEGACY_METRICS = {
    "sensorsWaterLevel_avg": {
        "value": [ 0.23116514, 0.06056554, 0.03638269, 0.02139411 ],
        "atol": 0.1,
        "rtol": 0.0,
    },
    "sensorsWaterLevel_std": {
        "value": [ 0.04803559, 0.03703002, 0.02929507, 0.02061927 ],
        "atol": 0.01,
        "rtol": 0.0,
    },
    "sensorsPressure_avg": {
        "value": [ 1237.25545429, 914.68332571, 682.51443457, 187.89957863 ],
        "atol": 100,
        "rtol": 0.0,
    },
    "sensorsPressure_std": {
        "value": [ 1920.12293911, 1379.96908324, 984.69733283, 521.99612149 ],
        "atol": 100,
        "rtol": 0.0,
    },
}


def _legacy_temp_ref() -> dict:
    """TEMP reference dict seeded from the legacy arrays (selftest helper)."""
    import datetime

    return {
        "case-tag": "SELFTEST/selftest:temp",
        "source": "selftest temp seed from evaluateExamplesMetrics.py:13-17",
        "date": datetime.date.today().isoformat(),
        "verified": False,
        "metrics": _LEGACY_METRICS,
    }


def _write_sensor_file( path: Path, avgs, stds, epsilon: float, n_rows: int = 8 ) -> None:
    """Synthetic sensor dat: columns ``avg +- (std + epsilon)`` alternating.

    Per-column mean equals ``avg`` exactly (the +- branches cancel) and the
    population std equals ``std + epsilon``; with ``epsilon << atol`` both
    comparator sub-checks lie well within tolerance.  First column is time.
    """
    rows = []
    for k in range( n_rows ):
        sign = 1.0 if k % 2 == 0 else -1.0
        values = [ a + sign * ( s + epsilon ) for a, s in zip( avgs, stds ) ]
        rows.append( " ".join( [ f"{0.1 * k:.17g}" ] + [ f"{v:.17g}" for v in values ] ) )
    path.write_text( "\n".join( rows ) + "\n", encoding="utf-8" )


def selftest() -> int:
    """Run the module selftest; returns 0 only if every behavior is shown."""
    import contextlib
    import io
    import json
    import tempfile

    import numpy as np

    from . import references, sensors

    # (a) registry inventory ------------------------------------------------
    names = sorted( REGISTRY )
    print( f"REGISTRY ({len( names )} entries):" )
    for name in names:
        print( f"  {name}" )
    assert len( REGISTRY ) == 11, f"REGISTRY must have 11 entries, has {len( REGISTRY )}"
    assert set( REGISTRY ) == _EXPECTED_NAMES, "REGISTRY names deviate from the legacy 10 + empty"
    print( "registry names + count 11: OK" )

    # (b) sanitize_tag + reference_path asserts ------------------------------
    assert (
        references.sanitize_tag( "WCSPH-DBC/damBreak2D_WCSPH-DBC:MGVT" )
        == "WCSPH-DBC/damBreak2D_WCSPH-DBC__MGVT"
    )
    assert (
        references.sanitize_tag( "WCSPH-DBC/damBreak2D_WCSPH-DBC" )
        == "WCSPH-DBC/damBreak2D_WCSPH-DBC"
    )
    assert references.reference_path( "WCSPH-DBC/damBreak2D_WCSPH-DBC:MGVT" ) == (
        references.REFERENCES_DIR / "WCSPH-DBC" / "damBreak2D_WCSPH-DBC__MGVT.json"
    )
    assert references.reference_path( "WCSPH-DBC/damBreak2D_WCSPH-DBC:MGVT", staging=True ) == (
        references.REFERENCES_DIR / "_staging" / "WCSPH-DBC" / "damBreak2D_WCSPH-DBC__MGVT.json"
    )
    print( "sanitize_tag + reference_path assertions: OK" )

    committed = references.load_reference( references.reference_path( _COMMITTED_TAG ) )

    # TEMP seed from the legacy arrays into _staging/ (allowed by the guard),
    # round-trip equality with the committed seeded file, then clean up.
    staging_tag = "SELFTEST/selftest:temp"
    staging_path = references.reference_path( staging_tag, staging=True )
    assert not references.staging_refs_only( staging_tag )
    references.save_reference( _legacy_temp_ref(), staging_path )
    try:
        assert references.staging_refs_only( staging_tag ), (
            "staging_refs_only must be True when only a staging ref exists"
        )
        assert not references.staging_refs_only( _COMMITTED_TAG ), (
            "staging_refs_only must be False for the committed seeded tag"
        )
        assert references.load_reference( staging_path )[ "metrics" ] == committed[ "metrics" ], (
            "temp seed from legacy arrays must reproduce the committed seeded metrics"
        )
        print( f"temp seed round-trip ({staging_path}): OK" )
    finally:
        staging_path.unlink( missing_ok=True )
        parent = staging_path.parent
        while parent != references.REFERENCES_DIR and parent.is_dir():
            try:
                parent.rmdir()
            except OSError:
                break
            parent = parent.parent
    assert not references.staging_refs_only( staging_tag )

    # save_reference must refuse targets outside references/.
    with tempfile.TemporaryDirectory() as tmp:
        outside = Path( tmp ) / "nope.json"
        try:
            references.save_reference( _legacy_temp_ref(), outside )
        except ValueError as exc:
            assert str( outside ) in str( exc )
            print( f"save_reference outside references/ raised as required: {exc}" )
        else:
            raise AssertionError( "save_reference accepted a path outside references/" )
        assert not outside.exists()

    # (c) comparator on synthetic sensor files -------------------------------
    metrics = committed[ "metrics" ]
    wl_avg, wl_std = ( metrics[ "sensorsWaterLevel_avg" ], metrics[ "sensorsWaterLevel_std" ] )
    ps_avg, ps_std = ( metrics[ "sensorsPressure_avg" ], metrics[ "sensorsPressure_std" ] )
    # epsilon = 1e-3 of the std atol per sensor (tiny vs every atol in play).
    eps_wl, eps_ps = 1e-3 * wl_std[ "atol" ], 1e-3 * ps_std[ "atol" ]
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = Path( tmp )
        results_dir = case_dir / "results"
        results_dir.mkdir()
        wl_file = results_dir / "sensorsWaterLevel.dat"
        ps_file = results_dir / "sensorsPressure.dat"

        print(
            "synthetic pass-case inputs: columns avg +- (std + eps) alternating"
            f" over 8 rows, eps_wl={eps_wl}, eps_ps={eps_ps}"
        )
        _write_sensor_file( wl_file, wl_avg[ "value" ], wl_std[ "value" ], eps_wl )
        _write_sensor_file( ps_file, ps_avg[ "value" ], ps_std[ "value" ], eps_ps )
        print( "sensorsWaterLevel.dat (pass case):" )
        print( wl_file.read_text( encoding="utf-8" ), end="" )
        print( "sensorsPressure.dat (pass case):" )
        print( ps_file.read_text( encoding="utf-8" ), end="" )

        result = REGISTRY[ "damBreak2D_WCSPH_DBC" ]( case_dir )
        print( f"pass-case comparator result: {result} (expect (4, 4))" )
        assert result == ( 4, 4 ), f"pass case must return (4, 4), got {result}"

        # Variant tag with NO committed per-tag reference -> base fallback.
        result_variant = REGISTRY[ "damBreak2D_WCSPH_DBC" ](
            case_dir, case_tag="WCSPH-DBC/damBreak2D_WCSPH-DBC:no-such-variant"
        )
        print(
            f"variant-tag base-fallback result: {result_variant} (expect (4, 4))"
        )
        assert result_variant == ( 4, 4 ), (
            f"variant tag without a per-tag reference must fall back to the base "
            f"reference and return (4, 4), got {result_variant}"
        )

        # Perturb column 0 of the water-level file by +5*atol(avg) -> mean out.
        shift = 5.0 * wl_avg[ "atol" ]
        _write_sensor_file(
            wl_file,
            [ wl_avg[ "value" ][ 0 ] + shift ] + wl_avg[ "value" ][ 1: ],
            wl_std[ "value" ],
            eps_wl,
        )
        print(
            f"perturbed-case input: sensorsWaterLevel.dat column 0 mean shifted by +{shift}"
        )
        print( wl_file.read_text( encoding="utf-8" ), end="" )
        capture = io.StringIO()
        with contextlib.redirect_stdout( capture ):
            result_perturbed = REGISTRY[ "damBreak2D_WCSPH_DBC" ]( case_dir )
        diagnostics = capture.getvalue()
        print( diagnostics, end="" )
        print( f"perturbed comparator result: {result_perturbed} (expect (3, 4))" )
        assert result_perturbed == ( 3, 4 ), f"perturbed case must return (3, 4), got {result_perturbed}"
        assert "sensorsWaterLevel_avg" in diagnostics and "failed" in diagnostics, (
            "diagnostics must name the failing metric"
        )
        print( "failing metric named in diagnostics: OK" )

        # collect mode returns the canonical metric names the comparator
        # consumes, with JSON-serializable values.
        collected = REGISTRY[ "damBreak2D_WCSPH_DBC" ]( case_dir, collect=True )
        assert set( collected ) == {
            "sensorsWaterLevel_avg",
            "sensorsWaterLevel_std",
            "sensorsPressure_avg",
            "sensorsPressure_std",
        }
        json.dumps( collected )
        print( f"collect-mode keys {sorted( collected )} JSON-serializable: OK" )

        # stubs keep exact legacy semantics.
        for name in sorted( REGISTRY ):
            if name == "damBreak2D_WCSPH_DBC":
                continue
            assert REGISTRY[ name ]( case_dir ) == ( 0, 0 ), name
            assert REGISTRY[ name ]( case_dir, collect=True ) == {}, name
        print( "9 stubs + empty return (0, 0) / {}: OK" )

        # unknown sensor file -> clean errors (missing file, missing metric).
        try:
            sensors.sensor_avg_std_check( case_dir, "noSuchSensor.dat", committed )
        except FileNotFoundError as exc:
            assert "noSuchSensor.dat" in str( exc )
            print( f"missing sensor file raised: {exc}" )
        else:
            raise AssertionError( "missing sensor file must raise" )
        try:
            sensors.sensor_avg_std_check( case_dir, "sensorsWaterLevel.dat", { "metrics": {} } )
        except ValueError as exc:
            assert "sensorsWaterLevel_avg" in str( exc )
            print( f"missing metric in reference raised: {exc}" )
        else:
            raise AssertionError( "missing metric in reference must raise" )

    print( "regression.metrics selftest: ALL OK - pass case (4, 4), perturbed case (3, 4)" )
    return 0
