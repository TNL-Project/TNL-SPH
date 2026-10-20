"""
tests.py

Synthetic self-test suite for the init-report validators.

Run:
    PYTHONPATH=src/Tests python3 -m infrastructure.validators.tests

Covers:
    (a) clean synthetic reports pass ALL rules - three scenarios
        (single-subdomain, nested multiresolution, distributed rank pair);
    (b) >=6 seeded one-field corruptions, each flipping ONLY its own rule;
    (c) schema_version 0 triggers schema.version_is_1;
    (d) deep_compare behavioral parity with the legacy engine in
        run_sph_tests.py on two sample documents.

Exit code 0 iff every assertion holds.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

from . import compare
from .invariants import (
    BUFFERS_MR,
    DECOMP_TILING,
    DOMAIN_GRID_ORIGIN,
    PATCHES_OPEN,
    PATCHES_PERIODIC,
    RULES,
    SCHEMA_VERSION_IS_1,
    SETS_COUNTS,
    SETS_HALO,
)

# ── Synthetic toy case (2D) ──────────────────────────────────────────────────

_ORIGIN = [0.0, 0.0, 0.0]
_SIZE = [1.0, 0.5, 0.0]
_SEARCH_RADIUS = 0.01
_DIMS = [100, 50, 0]  # ceil((size - origin) / search_radius) per axis
_FLUID_N = 1600
_BOUNDARY_N = 400


def _set_entry(count, allocated, dims, overlap, origin=None, origin_coords=None):
    origin = list(origin) if origin is not None else list(_ORIGIN)
    origin_coords = (list(origin_coords) if origin_coords is not None
                     else [0, 0, 0])
    return {
        "count": count,
        "allocated_count": allocated,
        "grid_origin": origin,
        "grid_origin_with_overlap": [
            origin[k] - overlap * _SEARCH_RADIUS if k < 2 else 0.0
            for k in range(3)],
        "grid_dimensions": list(dims),
        "grid_dimensions_with_overlap": [
            dims[k] + 2 * overlap if k < 2 else 0 for k in range(3)],
        "grid_origin_global_coordinates": origin_coords,
        "overlap_width": overlap,
    }


def _open_patch_entry():
    return {
        "count": 40,
        "allocated_count": 80,
        "grid_origin": [0.0, 0.1, 0.0],
        "grid_origin_with_overlap": [0.0, 0.1, 0.0],
        "grid_dimensions": [5, 20, 0],
        "grid_dimensions_with_overlap": [5, 20, 0],
        "grid_origin_global_coordinates": [0, 10, 0],
        "overlap_width": 0,
        "config": {
            "identifier": "inlet",
            "type": "inlet",
            "orientation_x": 1.0, "orientation_y": 0.0, "orientation_z": 0.0,
            "position_1_x": 0.0, "position_1_y": 0.1, "position_1_z": 0.0,
            "position_2_x": 0.0, "position_2_y": 0.3, "position_2_z": 0.0,
            "bufferWidth_x": 0.05, "bufferWidth_y": 0.0, "bufferWidth_z": 0.0,
            "paired_periodic_buffer": -1,
        },
    }


def _periodic_entry(ident, paired):
    return {
        "zone_cell_count": 100,
        "zone_particle_count": 0,
        "config": {
            "identifier": ident,
            "orientation_x": 1.0, "orientation_y": 0.0, "orientation_z": 0.0,
            "position_1_x": 0.0, "position_1_y": 0.1, "position_1_z": 0.0,
            "position_2_x": 0.0, "position_2_y": 0.3, "position_2_z": 0.0,
            "bufferWidth_x": 0.05, "bufferWidth_y": 0.0, "bufferWidth_z": 0.0,
            "paired_periodic_buffer": paired,
        },
    }


def _mr_buffer_entry(level, own_dims, own_origin, front_origin, front_dims,
                     back_origin, back_dims, resolution_factor):
    entry = _set_entry(count=0, allocated=1000, dims=own_dims, overlap=1,
                       origin=own_origin)
    entry.update({
        "level": level,
        "resolution_factor": resolution_factor,
        "frame_width": 2,
        "frame_front_origin_global_coordinates": list(front_origin),
        "frame_front_dimensions": list(front_dims),
        "frame_back_origin_global_coordinates": list(back_origin),
        "frame_back_dimensions": list(back_dims),
    })
    return entry


def make_base_report():
    return {
        "schema_version": 1,
        "config": {
            "domainOrigin": list(_ORIGIN),
            "domainSize": list(_SIZE),
            "searchRadius": _SEARCH_RADIUS,
            "dp": 0.005,
            "h": 0.007,
            "numberOfSubdomains": 1,
        },
        "fluids": [_set_entry(_FLUID_N, _FLUID_N, _DIMS, overlap=2)],
        "boundaries": [_set_entry(_BOUNDARY_N, _BOUNDARY_N, _DIMS, overlap=2)],
        "open_boundary_patches": [_open_patch_entry()],
        "periodic_patches": [_periodic_entry("periodicityLeft", 2),
                             _periodic_entry("periodicityRight", 1)],
        "multiresolution_buffers": [],
        "decomposition": {},
    }


def make_mr_report():
    report = make_base_report()
    report["config"]["numberOfSubdomains"] = 2
    report["open_boundary_patches"] = []
    report["periodic_patches"] = []
    report["multiresolution_buffers"] = [
        # level-1 band owned by the coarse (whole-domain) grid
        _mr_buffer_entry(level=1, own_dims=[100, 50, 0], own_origin=_ORIGIN,
                         front_origin=[48, 10, 0], front_dims=[12, 10, 0],
                         back_origin=[49, 11, 0], back_dims=[10, 8, 0],
                         resolution_factor=0.5),
        # level-2 band owned by the fine patch; frame coords in fine cells
        _mr_buffer_entry(level=2, own_dims=[24, 20, 0],
                         own_origin=[0.5, 0.1, 0.0],
                         front_origin=[96, 20, 0], front_dims=[24, 20, 0],
                         back_origin=[95, 19, 0], back_dims=[26, 22, 0],
                         resolution_factor=2.0),
    ]
    report["decomposition"] = {
        "global": {"origin": list(_ORIGIN), "dimensions": list(_DIMS)},
        "subdomains": [
            {"index": 0, "origin_global_coordinates": [0, 0, 0],
             "region_dimensions": [100.0, 50.0, 0.0]},
            {"index": 1, "origin_global_coordinates": [100, 20, 0],
             "region_dimensions": [12.0, 10.0, 0.0]},
        ],
    }
    return report


def make_rank_pair():
    rank0 = make_base_report()
    rank0["config"]["numberOfSubdomains"] = 2
    rank0["open_boundary_patches"] = []
    rank0["periodic_patches"] = []
    rank0["fluids"] = [_set_entry(800, 800, [50, 50, 0], overlap=2,
                                  origin_coords=[0, 0, 0])]
    rank0["boundaries"] = [_set_entry(200, 200, [50, 50, 0], overlap=2,
                                      origin_coords=[0, 0, 0])]
    rank0["decomposition"] = {
        "global": {"origin": list(_ORIGIN), "dimensions": list(_DIMS)},
        "subdomains": [
            {"index": 0, "origin_global_coordinates": [0, 0, 0],
             "region_dimensions": [50.0, 50.0, 0.0]},
            {"index": 1, "origin_global_coordinates": [50, 0, 0],
             "region_dimensions": [50.0, 50.0, 0.0]},
        ],
    }
    rank1 = copy.deepcopy(rank0)
    rank1["fluids"] = [_set_entry(800, 800, [50, 50, 0], overlap=2,
                                  origin=[0.5, 0.0, 0.0],
                                  origin_coords=[50, 0, 0])]
    rank1["boundaries"] = [_set_entry(200, 200, [50, 50, 0], overlap=2,
                                      origin=[0.5, 0.0, 0.0],
                                      origin_coords=[50, 0, 0])]
    return rank0, rank1


def make_expectations(*, scenario="single"):
    expectations = {
        "config_path": Path("<synthetic>"),
        "case_dir": Path("<synthetic>"),
        "domain_origin": list(_ORIGIN),
        "domain_size": list(_SIZE),
        "search_radius": _SEARCH_RADIUS,
        "dp": 0.005,
        "h": 0.007,
        "number_of_subdomains": 1,
        "fluid_n": _FLUID_N,
        "boundary_n": _BOUNDARY_N,
        "fluid_vtk": "sources/toy_fluid.vtk",
        "boundary_vtk": "sources/toy_boundary.vtk",
        "output_directory": "results",
        "expected_grid_dims": list(_DIMS),
        "open_boundary": [{
            "index": 1, "type": "inlet", "identifier": "inlet",
            "particles": "sources/inlet.vtk", "number_of_particles": 40,
            "orientation": [1.0, 0.0, 0.0],
            "position_1": [0.0, 0.1, 0.0],
            "position_2": [0.0, 0.3, 0.0],
            "buffer_width": [0.05, 0.0, 0.0],
            "paired_periodic_buffer": -1,
            "shift_vector": [0.0, 0.0, 0.0],
        }],
        "open_boundary_declared": True,
        "open_boundary_error": None,
        "periodic_boundary": None,
        "periodic_boundary_declared": False,
        "periodic_boundary_error": None,
        "periodic_patches_expected": 2,
        "subdomains": None,
        "subdomains_declared": False,
        "subdomains_error": None,
        "vtk_counts": {"fluid": _FLUID_N, "boundary": _BOUNDARY_N,
                       "fluid_file": "<toy>", "boundary_file": "<toy>",
                       "error": None},
        "rank_reports": None,
        "distributed": False,
    }
    if scenario == "mr":
        expectations["number_of_subdomains"] = 2
        expectations["open_boundary"] = None
        expectations["open_boundary_declared"] = False
        expectations["periodic_patches_expected"] = 0
    elif scenario == "distributed":
        expectations["number_of_subdomains"] = 2
        expectations["open_boundary"] = None
        expectations["open_boundary_declared"] = False
        expectations["periodic_patches_expected"] = 0
        rank0, rank1 = make_rank_pair()
        expectations["rank_reports"] = [rank0, rank1]
        expectations["distributed"] = True
    return expectations


# ── Assertion helpers ────────────────────────────────────────────────────────

_FAILURES: list[str] = []


def _check(condition: bool, label: str):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {label}")
    if not condition:
        _FAILURES.append(label)


def _run_all_rules(report, expectations):
    """Return {rule_name: Finding} for every rule that fired."""
    from .invariants import RULE_NAMES
    fired = {}
    for rule in RULES:
        finding = rule(report, expectations)
        if finding is not None:
            fired[RULE_NAMES[rule]] = finding
    return fired


# ── (a) clean scenarios ─────────────────────────────────────────────────────

def _test_clean_scenarios():
    print("(a) clean synthetic reports pass all rules")
    scenarios = [
        ("single-subdomain", make_base_report, make_expectations()),
        ("multiresolution-nested", make_mr_report,
         make_expectations(scenario="mr")),
    ]
    rank0, _rank1 = make_rank_pair()
    scenarios.append(("distributed-rank0", lambda: rank0,
                      make_expectations(scenario="distributed")))

    for name, report_fn, expectations in scenarios:
        fired = _run_all_rules(report_fn(), expectations)
        _check(not fired, f"clean {name} report: no rule fired"
               + (f" (got {list(fired)})" if fired else ""))


# ── (b)+(c) seeded one-field corruptors ──────────────────────────────────────

def _corrupt_open_patch_position(report):
    report["open_boundary_patches"][0]["config"]["position_1_x"] = 5.0


def _corrupt_periodic_removed(report):
    report["periodic_patches"] = []


def _corrupt_fluid_count(report):
    report["fluids"][0]["count"] += 1


def _corrupt_boundary_overlap(report):
    report["boundaries"][0]["overlap_width"] += 1


def _corrupt_grid_origin(report):
    report["fluids"][0]["grid_origin"][0] = 0.5


def _corrupt_schema_0(report):
    report["schema_version"] = 0


def _corrupt_mr_level(report):
    report["multiresolution_buffers"][1]["level"] = 1


def _corrupt_rank_origin(report):
    report["rank_reports"][1]["fluids"][0][
        "grid_origin_global_coordinates"][0] = 51


def _test_corruptors():
    print("(b) seeded one-field corruptions flip exactly their own rule")
    cases = [
        ("corrupt grid_origin", make_base_report, _corrupt_grid_origin,
         make_expectations(), DOMAIN_GRID_ORIGIN),
        ("corrupt fluids[0].count", make_base_report, _corrupt_fluid_count,
         make_expectations(), SETS_COUNTS),
        ("corrupt boundaries[0].overlap_width", make_base_report,
         _corrupt_boundary_overlap, make_expectations(), SETS_HALO),
        ("corrupt open patch position_1_x", make_base_report,
         _corrupt_open_patch_position, make_expectations(), PATCHES_OPEN),
        ("corrupt periodic patches removed", make_base_report,
         _corrupt_periodic_removed, make_expectations(), PATCHES_PERIODIC),
        ("corrupt multiresolution_buffers[1].level", make_mr_report,
         _corrupt_mr_level, make_expectations(scenario="mr"), BUFFERS_MR),
    ]
    print("(c) schema_version 0 triggers the schema rule")
    cases.append(("corrupt schema_version -> 0", make_base_report,
                  _corrupt_schema_0, make_expectations(),
                  SCHEMA_VERSION_IS_1))

    for label, report_fn, corrupt, expectations, expected_rule in cases:
        report = report_fn()
        corrupt(report)
        fired = _run_all_rules(report, expectations)
        _check(list(fired) == [expected_rule],
               f"{label}: only {expected_rule} fired"
               + (f" (got {list(fired)})" if fired else ""))
        if expected_rule == SCHEMA_VERSION_IS_1 and fired:
            _check("regenerate" in fired[expected_rule].message,
                   "schema finding message mentions regeneration")

    # distributed tiling lives on expectations["rank_reports"], not the
    # primary report - corrupt through a wrapper
    expectations = make_expectations(scenario="distributed")
    rank0, _rank1 = make_rank_pair()
    corrupt_report = {"rank_reports": expectations["rank_reports"]}
    _corrupt_rank_origin(corrupt_report)
    fired = _run_all_rules(rank0, expectations)
    _check(list(fired) == [DECOMP_TILING],
           f"corrupt rank-1 interior origin: only {DECOMP_TILING} fired"
           + (f" (got {list(fired)})" if fired else ""))


# ── (d) deep_compare parity with the legacy engine ───────────────────────────

def _legacy_diffs(expected, actual):
    infra_dir = Path(__file__).resolve().parents[1]
    if str(infra_dir) not in sys.path:
        sys.path.insert(0, str(infra_dir))
    import run_sph_tests  # noqa: PLC0415 - legacy engine, imported lazily
    return run_sph_tests.deep_compare(expected, actual)


def _diff_signature(diffs):
    return [(d.path, repr(d.expected), repr(d.actual)) for d in diffs]


def _test_deep_compare_parity():
    print("(d) deep_compare parity with run_sph_tests' engine")
    samples = [
        ({"a": 1, "b": [1.0, 2.0, 3.0], "c": {"d": "x", "e": 5}},
         {"a": 2, "b": [1.0, 2.0], "c": {"d": "x"}, "extra": True}),
        ({"x": 1.0000000005, "y": [0.5, "s"], "z": None},
         {"x": 1.0, "y": [0.5, 7], "z": 0}),
    ]
    for i, (expected, actual) in enumerate(samples):
        ours = compare.deep_compare(expected, actual)
        legacy = _legacy_diffs(expected, actual)
        _check(_diff_signature(ours) == _diff_signature(legacy),
               f"sample {i + 1}: identical diff lists "
               f"({len(ours)} ours vs {len(legacy)} legacy)")
    ours0 = compare.deep_compare(samples[0][0], samples[0][1])
    _check(len(ours0) == 3,
           f"sample 1 diff count is 3 (a, b.length, c.e) - got {len(ours0)}")


# ── Entry point ──────────────────────────────────────────────────────────────

def main() -> int:
    print("infrastructure.validators synthetic test suite")
    _test_clean_scenarios()
    _test_corruptors()
    _test_deep_compare_parity()
    if _FAILURES:
        print(f"\n{len(_FAILURES)} assertion(s) FAILED")
        return 1
    print("\nall assertions passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
