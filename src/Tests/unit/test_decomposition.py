"""Unit checks for decomposition math (group ``decomposition``).

Covers ``src/tools/decomposition.py::build_subdomain_grids`` (uniform 2D
tiling, 2D+MR rectangular refinement, 3D box), ``src/tools/splitToSubdomains.py``
(``Subdomains.divideIntoSubdomains`` x-splits and the per-subdomain config
record, incl. the known z-dimension defect at ``splitToSubdomains.py:17``)
and ``src/tools/generateCaseTool.py::getDomainLimits`` (margin/size math).

Fixture numbers are dyadic rationals (search radius 0.125, bounds in
eighths) so every ``floor``/``ceil`` inside the module under test is exact.

``build_subdomain_grids`` builds ONE grid per ``SubdomainDef`` - tilings are
expressed by passing several ``factor=1.0`` defs partitioning the axis; the
chain-contiguity rule below is the coverage invariant the plan requires
(a seeded +1dp fixture shift must flip it to FAIL by name).

Contract: ``check_*`` returns ``None`` on PASS or a message naming the
broken rule on FAIL (the runner turns raises into FAILs with a traceback).
"""

from __future__ import annotations

import contextlib
import math
import os
import re
import tempfile


def prepare():
    """Import the modules under test via the suite bootstrap."""
    import suite.bootstrap  # noqa: F401  (sys.path hook for src/tools)
    import numpy as np
    import decomposition
    import splitToSubdomains
    import generateCaseTool
    import saveParticlesVTK  # noqa: F401  (imported by generateSubdomain internally)

    return {
        "decomposition": decomposition,
        "splitToSubdomains": splitToSubdomains,
        "generateCaseTool": generateCaseTool,
        "np": np,
    }


@contextlib.contextmanager
def _cwd(path):
    """chdir into ``path``, always restoring the previous CWD."""
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


# Deterministic dyadic fixtures: every bound below is an integer multiple
# of _H0, so resolve_axis's floor() never lands mid-cell.
_H0 = 0.125
_DP = 0.125
_SETUP_2D = {
    "search_radius": _H0,
    "dp": _DP,
    "domain_origin_x": 0.0,
    "domain_origin_y": 0.0,
    "domain_size_x": 1.0,
    "domain_size_y": 0.5,
}
_SETUP_3D = {
    **_SETUP_2D,
    "domain_origin_z": 0.0,
    "domain_size_z": 0.375,
}


def _uniform_defs(decomposition):
    """Three ``factor=1.0`` defs partitioning x: [0,.375][.375,.625][.625,1]."""
    return [
        decomposition.SubdomainDef(factor=1.0, x_min=0.0, x_max=0.375),
        decomposition.SubdomainDef(factor=1.0, x_min=0.375, x_max=0.625),
        decomposition.SubdomainDef(factor=1.0, x_min=0.625, x_max=1.0),
    ]


# ---------------------------------------------------------------------------


def check_decomposition_uniform_tiling(ctx):
    """2D uniform tiling: tiles chain-cover x, propagate h0/dp per grid."""
    dec, np = ctx["decomposition"], ctx["np"]
    grids = dec.build_subdomain_grids(_uniform_defs(dec), _SETUP_2D)
    if len(grids) != 3:
        return f"uniform case: expected 3 grids from 3 defs, got {len(grids)}"

    tiles = sorted(grids, key=lambda g: g.phys_x_min)
    tol = 1e-12

    # Chain-contiguity + full cover along x.
    if abs(tiles[0].phys_x_min - _SETUP_2D["domain_origin_x"]) > tol:
        return f"x-cover: first tile starts at {tiles[0].phys_x_min}, not origin 0.0"
    last_max = tiles[-1].phys_x_max
    want_max = _SETUP_2D["domain_origin_x"] + _SETUP_2D["domain_size_x"]
    if abs(last_max - want_max) > tol:
        return f"x-cover: last tile ends at {last_max}, not {want_max}"
    for i in range(len(tiles) - 1):
        if abs(tiles[i].phys_x_max - tiles[i + 1].phys_x_min) > tol:
            return (
                f"x-tiling gap/overlap between tile {i} and {i + 1}: "
                f"tile[{i}].phys_x_max={tiles[i].phys_x_max} != "
                f"tile[{i + 1}].phys_x_min={tiles[i + 1].phys_x_min}"
            )
        if tiles[i].origin_glob_x + tiles[i].dims_x != tiles[i + 1].origin_glob_x:
            return (
                f"x-index chaining broken between tile {i} and {i + 1}: "
                f"{tiles[i].origin_glob_x}+{tiles[i].dims_x} != "
                f"{tiles[i + 1].origin_glob_x}"
            )

    # Monotone global origins; uniform resolution + dp propagation; 2D z None.
    origins = [g.origin_glob_x for g in tiles]
    if origins != sorted(origins):
        return f"origin_glob_x not monotone: {origins}"
    y0, y1 = _SETUP_2D["domain_origin_y"], _SETUP_2D["domain_origin_y"] + _SETUP_2D["domain_size_y"]
    for i, g in enumerate(tiles):
        if not np.isclose(g.search_radius, _H0):
            return f"tile {i}: search_radius {g.search_radius} != {_H0}"
        if not np.isclose(g.dp, _DP):
            return f"tile {i}: dp {g.dp} != {_DP}"
        if not (np.isclose(g.phys_y_min, y0) and np.isclose(g.phys_y_max, y1)):
            return f"tile {i}: y-span [{g.phys_y_min}, {g.phys_y_max}] != full [{y0}, {y1}]"
        if g.phys_z_min is not None or g.phys_z_max is not None:
            return f"tile {i}: 2D grid must carry phys_z=None, got {[g.phys_z_min, g.phys_z_max]}"
    return None


def check_decomposition_refinement_region(ctx):
    """2D+MR: refined-region grid carries finer, coarse-aligned dp/h0."""
    dec, np = ctx["decomposition"], ctx["np"]
    rdef = dec.RectangularRefinementDef(
        fine_factor=0.5,
        fine_x_min=0.25, fine_x_max=0.75,
        fine_y_min=0.125, fine_y_max=0.375,
    )
    coarse, fine = dec.build_rectangular_subdomain_grids(rdef, _SETUP_2D)
    tol = 1e-12

    # Per-grid arithmetic: dp = factor*setup["dp"], h = factor*setup["search_radius"].
    if not np.isclose(coarse.dp, _DP) or not np.isclose(coarse.search_radius, _H0):
        return f"coarse propagation broken: dp={coarse.dp}, h={coarse.search_radius}"
    if not np.isclose(fine.dp, 0.5 * coarse.dp):
        return f"refined dp {fine.dp} != 0.5 * coarse dp {coarse.dp}"
    if not np.isclose(fine.search_radius, 0.5 * coarse.search_radius):
        return f"refined h {fine.search_radius} != 0.5 * coarse h {coarse.search_radius}"

    # Coarse covers the whole domain; fine lies inside the declared box (and
    # hence inside coarse), with coarse-grid-aligned origins/dims (fl = 2).
    for ax, lo, hi in (("x", 0.0, 1.0), ("y", 0.0, 0.5)):
        if abs(getattr(coarse, f"phys_{ax}_min") - lo) > tol or abs(getattr(coarse, f"phys_{ax}_max") - hi) > tol:
            return f"coarse tile does not cover {ax}: [{getattr(coarse, f'phys_{ax}_min')}, {getattr(coarse, f'phys_{ax}_max')}]"
    for ax, lo, hi in (("x", 0.25, 0.75), ("y", 0.125, 0.375)):
        fmin, fmax = getattr(fine, f"phys_{ax}_min"), getattr(fine, f"phys_{ax}_max")
        if fmin < lo - tol or fmax > hi + tol:
            return f"fine tile escapes declared {ax}-region: [{fmin}, {fmax}] vs [{lo}, {hi}]"
        if getattr(fine, f"origin_glob_{ax}") % 2 != 0 or getattr(fine, f"dims_{ax}") % 2 != 0:
            return (
                f"fine grid not coarse-aligned on {ax}: "
                f"origin_glob={getattr(fine, f'origin_glob_{ax}')}, dims={getattr(fine, f'dims_{ax}')} (fl=2)"
            )
    # Firmly derivable dyadic alignment: fine x-span snaps exactly to the box.
    if abs(fine.phys_x_min - 0.25) > tol or abs(fine.phys_x_max - 0.75) > tol:
        return f"fine x-span [{fine.phys_x_min}, {fine.phys_x_max}] != [0.25, 0.75]"
    return None


def check_decomposition_3d_box(ctx):
    """3D box: z-axis propagated with phys/dims/origin; x chain-covers."""
    dec, np = ctx["decomposition"], ctx["np"]
    grids = dec.build_subdomain_grids(
        [
            dec.SubdomainDef(factor=1.0, x_min=0.0, x_max=0.25),
            dec.SubdomainDef(factor=1.0, x_min=0.25, x_max=0.5),
        ],
        _SETUP_3D,
    )
    if len(grids) != 2:
        return f"3D case: expected 2 grids, got {len(grids)}"
    tiles = sorted(grids, key=lambda g: g.phys_x_min)
    tol = 1e-12

    if abs(tiles[0].phys_x_max - tiles[1].phys_x_min) > tol:
        return (
            f"3D x-tiling gap/overlap: {tiles[0].phys_x_max} != {tiles[1].phys_x_min}"
        )
    for i, g in enumerate(tiles):
        if g.phys_z_min is None or g.phys_z_max is None:
            return f"3D grid {i}: phys_z must be resolved, got None"
        if abs(g.phys_z_min - 0.0) > tol or abs(g.phys_z_max - 0.375) > tol:
            return f"3D grid {i}: z-span [{g.phys_z_min}, {g.phys_z_max}] != [0.0, 0.375]"
        if g.dims_z != 3 or g.origin_glob_z != 0:
            return f"3D grid {i}: dims_z={g.dims_z}, origin_glob_z={g.origin_glob_z}; want 3/0"
        if not np.isclose(g.dp, _DP) or not np.isclose(g.search_radius, _H0):
            return f"3D grid {i}: propagation broken dp={g.dp}, h={g.search_radius}"
    return None


def _sorted_cloud(np):
    """30 rows [x,y,z] sorted by x with balanced mid-cell-safe splits.

    With 30 points / 3 subdomains the module splits at rows 10 and 20 via
    ``ceil(x / searchRadius)``; the values at those rows (0.95, 1.95) sit
    mid-cell, so float noise cannot flip the ceil.  Even pitch (10, 20)
    is required: the module's middle-subdomain width formula
    ``gridSplits[i-1] - gridIndexOrigins[i-1]`` equals the geometric width
    only for balanced splits, which is when ``sum(gridSizes) == gridSize_x``.
    """
    xs = []
    for i in range(30):
        if i < 10:
            xs.append(0.9 * (i + 1) / 11)          # 0.081..0.818
        elif i == 10:
            xs.append(0.95)                        # split -> ceil(9.5) = 10
        elif i < 20:
            xs.append(0.95 + (i - 10) * 0.1)       # 1.05..1.85
        elif i == 20:
            xs.append(1.95)                        # split -> ceil(19.5) = 20
        else:
            xs.append(1.95 + (i - 20) * 0.09)      # 2.04..2.76
    return np.array([[x, 0.01, 0.0] for x in xs], dtype=float)


def check_split_subdomains_tiling(ctx):
    """divideIntoSubdomains: x-splits tile [0,30] cells, origins monotone."""
    split, np = ctx["splitToSubdomains"], ctx["np"]
    sub = split.Subdomains(3)
    sub.divideIntoSubdomains(
        0.1, _sorted_cloud(np), [30, 10, 1], [0.0, 0.0, 0.0]
    )

    if sum(sub.gridSizes) != 30:
        return f"gridSizes {sub.gridSizes} sum to {sum(sub.gridSizes)}, not gridSize_x=30"
    if any(b < 0 or b > 30 for b in sub.gridSizes):
        return f"gridSizes out of [0,30]: {sub.gridSizes}"
    origins = sub.gridOrigins
    if any(origins[i] >= origins[i + 1] for i in range(len(origins) - 1)):
        return f"gridOrigins not strictly increasing: {origins}"
    if origins[0] != 0.0 or any(o > 0.0 + 30 * 0.1 + 1e-12 for o in origins):
        return f"gridOrigins out of physical [0.0, 3.0]: {origins}"
    idx = sub.gridIndexOrigins
    if idx[0] != 0 or any(idx[i] > idx[i + 1] for i in range(len(idx) - 1)):
        return f"gridIndexOrigins not monotone from 0: {idx}"
    if idx[-1] >= 30:
        return f"gridIndexOrigins out of [0,30): {idx}"
    # Physical/index consistency of every tile boundary.
    for i in range(3):
        if abs(origins[i] - (0.0 + idx[i] * 0.1)) > 1e-12:
            return f"tile {i}: origin {origins[i]} != begin + index*radius {idx[i] * 0.1}"
    return None


def check_split_subdomains_grid_record(ctx):
    """Config record: emitted gridSize triple must end in gridZsize, not gridYsize."""
    split, np = ctx["splitToSubdomains"], ctx["np"]
    sub = split.Subdomains(2)
    points_fluid = np.array(
        [[0.15, 0.01, 0.0], [2.55, 0.01, 0.0], [2.70, 0.01, 0.0]], dtype=float
    )
    points_box = np.array([[0.20, 0.0, 0.0], [2.60, 0.0, 0.0]], dtype=float)
    sub.divideIntoSubdomains(0.1, points_fluid, [30, 10, 1], [0.0, 0.0, 0.0])

    with tempfile.TemporaryDirectory() as tmp:
        os.mkdir(os.path.join(tmp, "sources"))
        with _cwd(tmp):
            info = sub.generateSubdomain(1, points_fluid, [], points_box, [])

    match = re.search(r"\.gridSize = \{\s*([^,]+?),\s*([^,]+?),\s*([^}]+?)\s*\};", info)
    if match is None:
        return "config record: no `gridSize = {X, Y, Z}` triple found"
    cx, cy, cz = (m.strip() for m in match.groups())
    if not cx.endswith("gridXsize"):
        return f"gridSize X component {cx!r} does not reference gridXsize"
    if not cy.endswith("gridYsize"):
        return f"gridSize Y component {cy!r} does not reference gridYsize"
    if not cz.endswith("gridZsize"):
        return (
            f"gridSize Z component is {cz!r} but must reference gridZsize "
            f"(subdomain record reuses gridYsize - splitToSubdomains.py:17 defect)"
        )
    # The record's scalar gridZsize must equal the z grid size handed in.
    zmatch = re.search(r"\.gridZsize = (\d+);", info)
    if zmatch is None or int(zmatch.group(1)) != 1:
        return f"gridZsize scalar missing or != 1: {zmatch.group(1) if zmatch else None}"
    return None


def check_domain_limits_margins(ctx):
    """getDomainLimits: begin <= min-r, end >= max+r, sizes == ceil-span/r."""
    gc, np = ctx["generateCaseTool"], ctx["np"]
    r = 0.1
    points_fluid = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [0.4, 0.1, 0.0]], dtype=float
    )
    points_box = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.5, 0.0]], dtype=float
    )
    gridBegin, gridEnd, gridSize = gc.getDomainLimits(
        points_fluid, points_box, r, 1.0, 1.0, 1.0
    )

    mins = [0.0, 0.0, 0.0]
    maxs = [1.0, 0.5, 0.0]
    for ax in range(3):
        # 1.005*(min-r) is below min-r whenever min-r <= 0 (all our fixtures).
        if gridBegin[ax] > mins[ax] - r + 1e-12:
            return f"axis {ax}: begin {gridBegin[ax]} not <= min-r {mins[ax] - r}"
        if gridEnd[ax] < maxs[ax] + r - 1e-12:
            return f"axis {ax}: end {gridEnd[ax]} not >= max+r {maxs[ax] + r}"
        want = math.ceil((gridEnd[ax] - gridBegin[ax]) / r - 1e-9)
        if gridSize[ax] != want:
            return f"axis {ax}: gridSize {gridSize[ax]} != ceil(span/r) {want}"
    return None
