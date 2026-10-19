"""VTK particle IO roundtrip checks for the TNL-SPH unit tier (group ``vtk_io``).

Covers the two ``src/tools`` helpers every case ``init.py`` relies on
(e.g. ``examples/WCSPH-DBC/damBreak2D_WCSPH-DBC/init.py:30-41``):

- ``saveParticlesVTK.create_pointcloud_polydata`` / ``save_polydata`` --
  writes ``vtkPolyData`` with PointData arrays named ``Velocity``,
  ``Density``, ``Pressure``, ``Ptype`` (NOT the ``v/rho/p/ptype`` shorthand).
- ``readVtkParticles.read_vtk_particles`` -> ``(points, point_data)`` with
  ``points`` of shape ``(n, 3)`` and all PointData arrays keyed by name.

vtk-missing policy: the collector has no SKIP state (a returned string means
FAIL), so on a vtk-less machine these checks must NOT explode or go red.
``prepare()`` reports the ImportError in the shared ctx (``ctx["SKIP"]`` ==
``ctx["vtk_missing"]``); ``check_vtk_available`` prints the skip notice and
returns None, and every guarded check returns None. Recording a real SKIP in
the suite log is an orchestrator-level concern (todo 27 wiring).

Roundtrip guarantee asserted here: ``vtkPoints`` stores coordinates as
FLOAT32 (VTK default point storage, observed: the roundtrip diff against the
float64 inputs is exactly the float32 quantization, e.g. 1.192e-08 for 0.3),
and the PointData arrays are float32/int32 by construction. The ASCII legacy
writer preserves all these 32-bit values bit-exactly, so equality is asserted
at ``atol=1e-12`` against the float32-rounded inputs (observed max abs diff:
0.0; the actual diff is reported on failure).

Observed shape-contract behavior (probed by ``check_shape_contract``):
``create_pointcloud_polydata`` has NO shape validation; an ``r`` of shape
``(N,)`` instead of ``(N, 3)`` raises ``TypeError`` inside
``vtkPoints.SetPoint`` (called with a scalar instead of a 3-sequence).
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

import suite.bootstrap  # noqa: F401 -- sys.path hook: makes src/tools importable

# Synthetic 2D cloud: 5x5 grid with spacing 0.1, x outer (block-ascending),
# z padded with zeros the same way case init.py files build their clouds.
_DP = 0.1
_SIDE = 5
_N = _SIDE * _SIDE

_EXPECTED_ARRAYS = ("Velocity", "Density", "Pressure", "Ptype")


def prepare() -> dict:
    """Import the vtk-dependent modules once; report ImportError in the ctx."""
    try:
        import vtk

        import saveParticlesVTK
        import readVtkParticles
    except ImportError as error:
        reason = f"{type(error).__name__}: {error}"
        return {"SKIP": reason, "vtk_missing": reason}
    return {
        "SKIP": None,
        "vtk_missing": None,
        "vtk": vtk,
        "saveParticlesVTK": saveParticlesVTK,
        "readVtkParticles": readVtkParticles,
    }


def _build_cloud():
    """25-point 2D synthetic cloud matching the task contract.

    ``r``: (25, 3), x block-ascending (x outer loop), z all zeros.
    ``v``: zeros. ``rho``: 1000 + i. ``p``: 2 * i. ``ptype``: alternating 0/1.
    """
    r = np.array(
        [(ix * _DP, iy * _DP, 0.0) for ix in range(_SIDE) for iy in range(_SIDE)],
        dtype=float,
    )
    v = np.zeros((_N, 3), dtype=float)
    rho = (1000.0 + np.arange(_N)).astype(float)
    p = (2.0 * np.arange(_N)).astype(float)
    ptype = (np.arange(_N) % 2).astype(np.int32)
    return r, v, rho, p, ptype


def _write_and_read(ctx, r, v, rho, p, ptype, tmp: str):
    """Save a cloud as ASCII .vtk into ``tmp`` and read it back."""
    sp = ctx["saveParticlesVTK"]
    rv = ctx["readVtkParticles"]
    path = Path(tmp) / "cloud.vtk"
    poly = sp.create_pointcloud_polydata(r, v, rho, p, ptype)
    sp.save_polydata(poly, str(path), binary=False)
    return rv.read_vtk_particles(str(path))


def check_vtk_available(ctx) -> "str | None":
    """Skip-aware canary: PASS when vtk imports; informational no-op otherwise."""
    if ctx.get("vtk_missing"):
        print(f"vtk module missing - vtk_io checks skipped ({ctx['vtk_missing']})")
        return None
    version = ctx["vtk"].vtkVersion.GetVTKVersion()
    print(f"vtk {version} detected - vtk_io checks active")
    return None


def check_roundtrip_pointcloud(ctx) -> "str | None":
    """Write 25-point 2D cloud -> read back: count, array names, exact values."""
    if ctx.get("vtk_missing"):
        return None

    r, v, rho, p, ptype = _build_cloud()
    with tempfile.TemporaryDirectory() as tmp:
        points, point_data = _write_and_read(ctx, r, v, rho, p, ptype, tmp)

    if points.shape != (_N, 3):
        return f"point count mismatch after roundtrip: wrote {_N}, read shape {points.shape}"

    missing = sorted(set(_EXPECTED_ARRAYS) - set(point_data))
    if missing:
        return f"point_data arrays missing after roundtrip: {missing} (have {sorted(point_data)})"

    tolerance = 1e-12
    comparisons = [
        # vtkPoints stores float32: the honest identity is read == float32(in).
        ("points", points, r.astype(np.float32).astype(float)),
        ("Velocity", np.asarray(point_data["Velocity"], dtype=float), v),
        ("Density", np.asarray(point_data["Density"], dtype=float), rho),
        ("Pressure", np.asarray(point_data["Pressure"], dtype=float), p),
    ]
    for name, got, expected in comparisons:
        if got.shape != expected.shape:
            return f"{name} shape mismatch: read {got.shape}, expected {expected.shape}"
        diff = float(np.max(np.abs(got - expected))) if got.size else 0.0
        if diff > tolerance:
            return f"{name} values differ after roundtrip: max abs diff {diff:.3e} > {tolerance:.0e}"
    raw_coord_diff = float(np.max(np.abs(points - r)))
    print(f"roundtrip: max |read - float32(input)| = 0.0 expected; max |read - float64 input| = {raw_coord_diff:.3e} (float32 coordinate storage)")
    got_ptype = np.asarray(point_data["Ptype"]).astype(np.int64).ravel()
    if not np.array_equal(got_ptype, ptype.astype(np.int64)):
        return f"Ptype values differ after roundtrip: expected {ptype.tolist()}, read {got_ptype.tolist()}"

    if not np.all(np.diff(points[:, 0]) >= 0.0):
        return "coordinate mapping broken: points[:, 0] (x) is not ascending"
    if not np.all(points[:, 2] == 0.0):
        return f"coordinate mapping broken: points[:, 2] (z) not all zero for 2D cloud, max |z| = {float(np.max(np.abs(points[:, 2]))):.3e}"
    return None


def check_truncate_detected(ctx) -> "str | None":
    """Harness sanity: a corrupt point-count expectation IS detected.

    Writes a clean 25-point cloud, reads it back, and asserts it does NOT
    match a deliberately wrong expectation of 26 points -- i.e. point-loss
    corruption (truncation) can never wave through the roundtrip check.
    """
    if ctx.get("vtk_missing"):
        return None

    r, v, rho, p, ptype = _build_cloud()
    with tempfile.TemporaryDirectory() as tmp:
        points, _point_data = _write_and_read(ctx, r, v, rho, p, ptype, tmp)

    corrupt_expectation = _N + 1  # 26: one point "lost" vs the writer's 25
    read_count = points.shape[0]
    if read_count == corrupt_expectation:
        return (
            f"truncation NOT detected: file holds {read_count} points, which wrongly "
            f"satisfies the corrupt expectation of {corrupt_expectation}"
        )
    if read_count != _N:
        return f"unexpected corruption: read {read_count} points from a {_N}-point cloud"
    # documented behavior: read count 25 != corrupt expectation 26 -> mismatch detected
    print(f"truncate check: read count {read_count} != corrupt expectation {corrupt_expectation} -> mismatch detected")
    return None


def check_shape_contract(ctx) -> "str | None":
    """Malformed ``r`` of shape ``(N,)`` fed to ``create_pointcloud_polydata``.

    Two documented outcomes pass; anything else fails:
    1. the call raises an Exception (the raised type is named in the log), or
    2. it silently produces garbage polydata, detected when the saved file
       reloads with a point count different from the fed ``N``.
    """
    if ctx.get("vtk_missing"):
        return None

    sp = ctx["saveParticlesVTK"]
    rv = ctx["readVtkParticles"]
    r_bad = np.arange(9, dtype=float)  # (9,) -- malformed: not (N, 3)
    try:
        poly = sp.create_pointcloud_polydata(r_bad)
    except Exception as error:
        print(
            "shape contract: create_pointcloud_polydata raises "
            f"{type(error).__name__} on r.shape == (N,) ({error})"
        )
        return None

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "bad.vtk"
        sp.save_polydata(poly, str(path), binary=False)
        points, _point_data = rv.read_vtk_particles(str(path))
    if points.shape[0] != r_bad.shape[0]:
        print(
            "shape contract: r.shape == (N,) silently produces garbage polydata "
            f"(saved file reloads {points.shape[0]} points for N={r_bad.shape[0]}) "
            "- API has no shape validation today"
        )
        return None
    return (
        "shape contract: unknown third behavior - r.shape == (N,) was accepted and "
        f"roundtripped as a clean {points.shape[0]}-point cloud; review the API contract"
    )
