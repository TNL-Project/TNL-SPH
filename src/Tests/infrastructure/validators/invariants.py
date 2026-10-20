"""
invariants.py

Derivable geometric/structural invariants of an ``init_report.json``
against the expectations built by
``infrastructure.validators.expectations``.

Every rule is a small named function ``(report, expectations) -> Finding | None``
collected in ``RULES`` at the bottom.  Rule names are the dotted strings
listed below - they are printed verbatim by the CLI and asserted by the
synthetic test suite, so they must not be renamed casually.

Comparison tolerance: report floats come from the solver's ``float``
(32-bit) storage while config files carry full double precision, so float
comparisons use ``REL_TOL = 1e-6`` relative to the expected value (about
8x the float32 machine epsilon).  Integer quantities (counts, dimensions,
coordinates) are compared exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional

REL_TOL = 1e-6

# Canonical rule names (printed in findings, asserted by the test suite).
SCHEMA_VERSION_IS_1 = "schema.version_is_1"
DOMAIN_GRID_DIMS = "domain.grid_dims_matches_config"
DOMAIN_GRID_ORIGIN = "domain.grid_origin_matches_config"
SETS_COUNTS = "sets.counts_match_vtk"
SETS_HALO = "sets.halo_consistency"
PATCHES_OPEN = "patches.open_boundary_geometry"
PATCHES_PERIODIC = "patches.periodic_present_if_configured"
BUFFERS_MR = "buffers.mr_frame_in_domain"
DECOMP_TILING = "decomposition.subdomain_tiling"


@dataclass
class Finding:
    rule: str
    message: str
    severity: str = "FAIL"


# ── Shared helpers ───────────────────────────────────────────────────────────

def _close(actual: float, expected: float, rel_tol: float = REL_TOL) -> bool:
    return abs(float(actual) - float(expected)) <= rel_tol * max(1.0, abs(float(expected)))


def _vec_close(actual, expected, rel_tol: float = REL_TOL) -> bool:
    return all(_close(a, e, rel_tol) for a, e in zip(actual, expected))


def _active_axes(expectations: dict) -> list[int]:
    """Axes with a nonzero domain extent in cells (2 out of 3 in 2D runs)."""
    dims = expectations.get("expected_grid_dims") or [0, 0, 0]
    return [k for k in range(3) if dims[k] > 0]


def _all_sets(report: dict):
    """Yield (section, index, entry) for every particle-set-shaped report entry."""
    for section in ("fluids", "boundaries", "multiresolution_buffers"):
        for i, entry in enumerate(report.get(section) or []):
            yield section, i, entry


# ── H3: schema gate ──────────────────────────────────────────────────────────

def rule_schema_version_is_1(report: dict, expectations: dict) -> Optional[Finding]:
    version = report.get("schema_version")
    if version != 1:
        return Finding(SCHEMA_VERSION_IS_1,
                       f"schema_version={version!r} - regenerate/init-report schema mismatch")
    return None


# ── Domain geometry ──────────────────────────────────────────────────────────

def rule_domain_grid_dims_matches_config(report: dict, expectations: dict) -> Optional[Finding]:
    if expectations.get("number_of_subdomains", 1) != 1:
        # distributed per-rank grids differ from the global ones by design;
        # the tiling rule owns the distributed geometry check.
        return None
    expected = list(expectations["expected_grid_dims"])
    axes = _active_axes(expectations)
    for section, i, entry in _all_sets(report):
        dims = entry.get("grid_dimensions")
        if dims is None:
            return Finding(DOMAIN_GRID_DIMS,
                           f"{section}[{i}] is missing grid_dimensions")
        bad = [k for k in axes if int(dims[k]) != expected[k]]
        if bad:
            return Finding(
                DOMAIN_GRID_DIMS,
                f"{section}[{i}].grid_dimensions={dims} != config-derived "
                f"{expected} (axes {bad}, origin={expectations['domain_origin']}, "
                f"size={expectations['domain_size']}, "
                f"searchRadius={expectations['search_radius']})")
    return None


def rule_domain_grid_origin_matches_config(report: dict, expectations: dict) -> Optional[Finding]:
    if expectations.get("number_of_subdomains", 1) != 1:
        return None
    expected = expectations["domain_origin"]
    for section, i, entry in _all_sets(report):
        origin = entry.get("grid_origin")
        if origin is None:
            return Finding(DOMAIN_GRID_ORIGIN,
                           f"{section}[{i}] is missing grid_origin")
        if not _vec_close(origin, expected):
            return Finding(
                DOMAIN_GRID_ORIGIN,
                f"{section}[{i}].grid_origin={origin} != config domainOrigin "
                f"{expected} (rel tol {REL_TOL})")
    return None


# ── Particle-set counts and halo ─────────────────────────────────────────────

def rule_sets_counts_match_vtk(report: dict, expectations: dict) -> Optional[Finding]:
    fluids = report.get("fluids") or []
    boundaries = report.get("boundaries") or []

    for section, i, entry in _all_sets(report):
        count = entry.get("count")
        allocated = entry.get("allocated_count")
        if count is None or allocated is None:
            return Finding(SETS_COUNTS,
                           f"{section}[{i}] is missing count/allocated_count")
        if int(allocated) < int(count):
            return Finding(
                SETS_COUNTS,
                f"{section}[{i}].allocated_count={allocated} < count={count}")

    counts = expectations.get("vtk_counts") or {}
    expected_fluid = counts.get("fluid")
    expected_boundary = counts.get("boundary")
    source = "VTK"
    if expected_fluid is None and expectations.get("fluid_n") is not None:
        expected_fluid = expectations["fluid_n"]
        expected_boundary = expectations.get("boundary_n")
        source = "config counts"
    if expected_fluid is None:
        return None  # nothing derivable to compare against
    if expectations.get("number_of_subdomains", 1) != 1:
        # distributed per-rank counts are subsets of any whole-domain VTK;
        # allocated>=count above still applied.
        return None

    # per-file comparison is meaningful only for single-set reports
    if len(fluids) == 1 and int(fluids[0]["count"]) != int(expected_fluid):
        return Finding(SETS_COUNTS,
                       f"fluids[0].count={fluids[0]['count']} != {source} "
                       f"fluid count {expected_fluid}")
    if len(boundaries) == 1 and expected_boundary is not None \
            and int(boundaries[0]["count"]) != int(expected_boundary):
        return Finding(SETS_COUNTS,
                       f"boundaries[0].count={boundaries[0]['count']} != {source} "
                       f"boundary count {expected_boundary}")
    return None


def rule_sets_halo_consistency(report: dict, expectations: dict) -> Optional[Finding]:
    axes = _active_axes(expectations)
    inactive = [k for k in range(3) if k not in axes]
    for section, i, entry in _all_sets(report):
        dims = entry.get("grid_dimensions")
        dwo = entry.get("grid_dimensions_with_overlap")
        overlap = entry.get("overlap_width")
        if dims is None or dwo is None or overlap is None:
            return Finding(
                SETS_HALO,
                f"{section}[{i}] is missing grid_dimensions / "
                f"grid_dimensions_with_overlap / overlap_width")
        for k in axes:
            if int(dwo[k]) - int(dims[k]) != 2 * int(overlap):
                return Finding(
                    SETS_HALO,
                    f"{section}[{i}] axis {k}: grid_dimensions_with_overlap - "
                    f"grid_dimensions = {int(dwo[k]) - int(dims[k])} != "
                    f"2*overlap_width ({2 * int(overlap)})")
        for k in inactive:
            if int(dwo[k]) != int(dims[k]):
                return Finding(
                    SETS_HALO,
                    f"{section}[{i}] inactive axis {k}: grid_dimensions_with_overlap "
                    f"{dwo[k]} != grid_dimensions {dims[k]}")
    return None


# ── Boundary patches ─────────────────────────────────────────────────────────

def rule_patches_open_boundary_geometry(report: dict, expectations: dict) -> Optional[Finding]:
    patches = report.get("open_boundary_patches") or []
    expected = expectations.get("open_boundary")

    if not expected:
        if patches:
            return Finding(
                PATCHES_OPEN,
                f"report has {len(patches)} open-boundary patch(es) but the "
                f"config declares no open-boundary loader")
        if expectations.get("open_boundary_declared"):
            return Finding(
                PATCHES_OPEN,
                f"config declares an open-boundary loader but it could not be "
                f"read: {expectations.get('open_boundary_error')}")
        return None

    if len(patches) != len(expected):
        return Finding(
            PATCHES_OPEN,
            f"report has {len(patches)} open-boundary patch(es), loader config "
            f"has {len(expected)}")

    origin = expectations["domain_origin"]
    size = expectations["domain_size"]

    def _in_domain(p) -> bool:
        return all(
            origin[k] - REL_TOL * max(1.0, abs(origin[k]))
            <= float(p[k])
            <= size[k] + REL_TOL * max(1.0, abs(size[k]))
            for k in range(3)
        )

    for exp, patch in zip(expected, patches):
        cfg = patch.get("config") or {}
        ident = exp.get("identifier") or f"buffer-{exp['index']}"
        for key, axis in (("position_1", None), ("position_2", None)):
            reported = [cfg.get(f"{key}_x"), cfg.get(f"{key}_y"), cfg.get(f"{key}_z")]
            if any(v is None for v in reported):
                return Finding(PATCHES_OPEN,
                               f"open patch '{ident}' is missing {key}_* entries")
            if not _vec_close(reported, exp[key]):
                return Finding(
                    PATCHES_OPEN,
                    f"open patch '{ident}' {key}={reported} != loader config "
                    f"{exp[key]} (rel tol {REL_TOL})")
            if not _in_domain(reported):
                return Finding(
                    PATCHES_OPEN,
                    f"open patch '{ident}' {key}={reported} lies outside the "
                    f"domain [{origin} .. {size}]")

        orientation = [cfg.get("orientation_x", 0.0),
                       cfg.get("orientation_y", 0.0),
                       cfg.get("orientation_z", 0.0)]
        width = [cfg.get("bufferWidth_x", 0.0),
                 cfg.get("bufferWidth_y", 0.0),
                 cfg.get("bufferWidth_z", 0.0)]
        dominant = max(range(3), key=lambda k: abs(float(orientation[k])))
        for k in range(3):
            nonzero = abs(float(width[k])) > REL_TOL
            if k == dominant and not nonzero:
                return Finding(
                    PATCHES_OPEN,
                    f"open patch '{ident}' bufferWidth is zero on the dominant "
                    f"orientation axis {'xyz'[dominant]} (width={width})")
            if k != dominant and nonzero:
                return Finding(
                    PATCHES_OPEN,
                    f"open patch '{ident}' bufferWidth[{ 'xyz'[k] }]={width[k]} is "
                    f"nonzero off the dominant orientation axis "
                    f"{'xyz'[dominant]}")
    return None


def rule_patches_periodic_present_if_configured(report: dict, expectations: dict) -> Optional[Finding]:
    expected_count = int(expectations.get("periodic_patches_expected", 0))
    patches = report.get("periodic_patches") or []
    if len(patches) != expected_count:
        return Finding(
            PATCHES_PERIODIC,
            f"report has {len(patches)} periodic patch(es), config expects "
            f"{expected_count}")
    return None


# ── Multiresolution buffers ──────────────────────────────────────────────────

def rule_buffers_mr_frame_in_domain(report: dict, expectations: dict) -> Optional[Finding]:
    buffers = report.get("multiresolution_buffers") or []
    if not buffers:
        return None
    global_dims = expectations["expected_grid_dims"]
    levels = []
    for i, entry in enumerate(buffers):
        level = entry.get("level")
        if level is None:
            return Finding(BUFFERS_MR, f"multiresolution_buffers[{i}] is missing level")
        levels.append(int(level))
        if int(level) >= 2 and not float(entry.get("resolution_factor", 0.0)) >= 2.0:
            return Finding(
                BUFFERS_MR,
                f"multiresolution_buffers[{i}] (level {level}) has "
                f"resolution_factor={entry.get('resolution_factor')} < 2 "
                f"expected for finer levels")

        # frame boxes are expressed in the buffer's OWN resolution:
        # level L cells are searchRadius / 2**(L-1), so the global domain
        # spans dims * 2**(L-1) own cells; frames are halos and may poke
        # past the domain edge by overlap_width own cells.
        scale = 2 ** (int(level) - 1)
        pad = int(entry.get("overlap_width", 0))
        for side in ("front", "back"):
            origin_key = f"frame_{side}_origin_global_coordinates"
            dims_key = f"frame_{side}_dimensions"
            frame_origin = entry.get(origin_key)
            frame_dims = entry.get(dims_key)
            if frame_origin is None or frame_dims is None:
                return Finding(BUFFERS_MR,
                               f"multiresolution_buffers[{i}] is missing "
                               f"{origin_key}/{dims_key}")
            for k in range(3):
                span = global_dims[k] * scale
                lo = int(frame_origin[k])
                hi = lo + int(frame_dims[k])
                if span == 0:
                    if lo != 0 or hi != 0:
                        return Finding(
                            BUFFERS_MR,
                            f"multiresolution_buffers[{i}] {side} frame axis {k}: "
                            f"[{lo}, {hi}] on a degenerate (2D-z) axis")
                elif lo < -pad or hi > span + pad:
                    return Finding(
                        BUFFERS_MR,
                        f"multiresolution_buffers[{i}] (level {level}) {side} "
                        f"frame axis {k}: [{lo}, {hi}] outside the domain span "
                        f"[-{pad}, {span + pad}] (global dims {global_dims}, "
                        f"resolution scale {scale}, overlap pad {pad})")
    if len(set(levels)) != len(levels):
        return Finding(BUFFERS_MR,
                       f"multiresolution buffer levels are not distinct: {levels}")
    return None


# ── Distributed decomposition ────────────────────────────────────────────────

def rule_decomposition_subdomain_tiling(report: dict, expectations: dict) -> Optional[Finding]:
    """
    Distributed path (runs only when ``expectations['distributed']`` is set,
    i.e. the CLI got ``--distributed N`` and loaded the per-rank reports
    into ``expectations['rank_reports']``).

    Mirrors splitToSubdomains.py:44-75 semantics: rank interiors
    (``origin_global_coordinates .. +grid_dimensions``) lie inside the
    config-derived global dims and union-cover them on every active axis.
    """
    if not expectations.get("distributed"):
        return None
    rank_reports = expectations.get("rank_reports") or []
    n_expected = int(expectations.get("number_of_subdomains", 1))
    if len(rank_reports) != n_expected:
        return Finding(
            DECOMP_TILING,
            f"{len(rank_reports)} rank report(s) loaded != "
            f"numberOfSubdomains={n_expected} from config")

    decomposition = report.get("decomposition") or {}
    subdomains = decomposition.get("subdomains") or []
    if subdomains and len(subdomains) != n_expected:
        return Finding(
            DECOMP_TILING,
            f"decomposition.subdomains has {len(subdomains)} entries != "
            f"numberOfSubdomains={n_expected}")

    global_dims = expectations["expected_grid_dims"]
    axes = _active_axes(expectations)
    intervals = []
    for r, rank in enumerate(rank_reports):
        fluids = rank.get("fluids") or []
        if not fluids:
            return Finding(DECOMP_TILING,
                           f"rank {r} report has no fluid set to locate")
        grid = fluids[0]
        origin_c = grid.get("grid_origin_global_coordinates")
        dims = grid.get("grid_dimensions")
        overlap = int(grid.get("overlap_width", 0))
        if origin_c is None or dims is None:
            return Finding(
                DECOMP_TILING,
                f"rank {r} fluid grid is missing "
                f"grid_origin_global_coordinates/grid_dimensions")
        for k in axes:
            lo, hi = int(origin_c[k]), int(origin_c[k]) + int(dims[k])
            if lo < 0 or hi > global_dims[k]:
                return Finding(
                    DECOMP_TILING,
                    f"rank {r} interior interval axis {k}: [{lo}, {hi}) outside "
                    f"the global dims [0, {global_dims[k]}) "
                    f"(overlap_width={overlap})")
            intervals.append((k, lo, hi, r))

    for k in axes:
        covered = [False] * global_dims[k]
        for axis, lo, hi, _r in intervals:
            if axis != k:
                continue
            for cell in range(max(0, lo), min(global_dims[k], hi)):
                covered[cell] = True
        missing = [idx for idx, ok in enumerate(covered) if not ok]
        if missing:
            return Finding(
                DECOMP_TILING,
                f"axis {k}: rank interiors do not union-cover the global dims "
                f"[0, {global_dims[k]}) - first uncovered cell {missing[0]} "
                f"({len(missing)} uncovered)")
    return None


# ── Rule registry ────────────────────────────────────────────────────────────

RULES: list[Callable[[dict, dict], Optional[Finding]]] = [
    rule_schema_version_is_1,
    rule_domain_grid_dims_matches_config,
    rule_domain_grid_origin_matches_config,
    rule_sets_counts_match_vtk,
    rule_sets_halo_consistency,
    rule_patches_open_boundary_geometry,
    rule_patches_periodic_present_if_configured,
    rule_buffers_mr_frame_in_domain,
    rule_decomposition_subdomain_tiling,
]

RULE_NAMES: dict[Callable, str] = {
    rule_schema_version_is_1: SCHEMA_VERSION_IS_1,
    rule_domain_grid_dims_matches_config: DOMAIN_GRID_DIMS,
    rule_domain_grid_origin_matches_config: DOMAIN_GRID_ORIGIN,
    rule_sets_counts_match_vtk: SETS_COUNTS,
    rule_sets_halo_consistency: SETS_HALO,
    rule_patches_open_boundary_geometry: PATCHES_OPEN,
    rule_patches_periodic_present_if_configured: PATCHES_PERIODIC,
    rule_buffers_mr_frame_in_domain: BUFFERS_MR,
    rule_decomposition_subdomain_tiling: DECOMP_TILING,
}
