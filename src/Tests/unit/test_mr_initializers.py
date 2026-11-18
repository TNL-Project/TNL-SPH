"""Unit checks for the 2D and 3D multiresolution dam-break initializers."""

from __future__ import annotations

import contextlib
import importlib.util
import os
import tempfile
from pathlib import Path

import suite.bootstrap  # noqa: F401


def _load_initializer(name, relative_path):
    """Load one case script without importing its sibling initializer."""
    from suite.paths import ROOT

    path = Path(ROOT) / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load initializer {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare():
    """Load both initializers and the real VTK reader."""
    import decomposition
    import numpy as np
    import readVtkParticles

    uniform = _load_initializer(
        "mr_uniform_initializer",
        "examples/WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution/init_multi_resolution.py",
    )
    lts = _load_initializer(
        "mr_lts_initializer",
        "examples/WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution_local-timestepping/init_multi_resolution_lts.py",
    )
    uniform_3d = _load_initializer(
        "mr_uniform_initializer_3d",
        "examples/WCSPH-BI/damBreak3D_WCSPH-BI_multiresolution/init_multi_resolution.py",
    )
    dummy_3d = _load_initializer(
        "mr_dummy_initializer_3d",
        "src/Tests/infrastructure/fixtures/dummyMultiresolutionSimulation3D/init_3dmr_configuration.py",
    )
    return {"dec": decomposition, "np": np, "reader": readVtkParticles,
            "uniform": uniform, "lts": lts, "uniform_3d": uniform_3d, "dummy_3d": dummy_3d}


@contextlib.contextmanager
def _cwd(path):
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _grid(dec, factor, x_min, x_max, y_min, y_max):
    return dec.SubdomainGrid(
        factor=factor,
        search_radius=0.1 * factor,
        dp=0.1 * factor,
        phys_x_min=x_min,
        phys_x_max=x_max,
        phys_y_min=y_min,
        phys_y_max=y_max,
    )


def _grids(dec):
    return [
        _grid(dec, 1.0, 0.0, 1.0, 0.0, 1.0),
        _grid(dec, 0.5, 0.5, 1.0, 0.0, 0.5),
        _grid(dec, 0.25, 0.75, 1.0, 0.0, 0.25),
    ]


def _setup():
    return {
        "fluid_length": 1.0,
        "fluid_height": 1.0,
        "box_length": 1.0,
        "box_height": 0.5,
        "n_boundary_layers": 1,
        "density": 1000.0,
        "speed_of_sound": 34.3,
    }


def _check_intermediate_exclusion(ctx, module):
    """When tiling three levels, L1 must not retain points inside L2."""
    if not hasattr(module, "generate_fluid_particles_nested"):
        return "initializer has no nested-chain fluid generator"
    with tempfile.TemporaryDirectory() as tmp:
        Path(tmp, "sources").mkdir()
        with _cwd(tmp):
            module.generate_fluid_particles_nested(_grids(ctx["dec"]), _setup())
            clouds = [
                ctx["reader"].read_vtk_particles(
                    f"sources/subdomain-{level}-dambreak_fluid.vtk"
                )[0]
                for level in range(3)
            ]
    inside_l1 = (
        (clouds[0][:, 0] > 0.5) & (clouds[0][:, 0] <= 1.0)
        & (clouds[0][:, 1] > 0.0) & (clouds[0][:, 1] <= 0.5)
    )
    if ctx["np"].any(inside_l1):
        return f"L0 overlaps L1 at {int(inside_l1.sum())} particles"
    inside_l2 = (
        (clouds[1][:, 0] > 0.75) & (clouds[1][:, 0] <= 1.0)
        & (clouds[1][:, 1] > 0.0) & (clouds[1][:, 1] <= 0.25)
    )
    if ctx["np"].any(inside_l2):
        return f"L1 overlaps L2 at {int(inside_l2.sum())} particles"
    if len(clouds[1]) == 0:
        return "intermediate level unexpectedly has no particles outside L2"
    if len(clouds[2]) == 0:
        return "finest level unexpectedly has no particles inside L2"
    return None


def check_uniform_intermediate_level_exclusion(ctx):
    """Given three uniform-step levels, adjacent boxes tile without overlap."""
    return _check_intermediate_exclusion(ctx, ctx["uniform"])


def check_lts_intermediate_level_exclusion(ctx):
    """Given three LTS levels, L2 is excluded from L1 as well as L1 from L0."""
    return _check_intermediate_exclusion(ctx, ctx["lts"])


def check_problem_box_serializes_far_corner(ctx):
    """Both initializers emit maximum corners, not lengths, as domainSize*."""
    for name in ("uniform", "lts"):
        setup = {"search_radius": 0.1, "box_length": 1.0, "box_height": 0.5}
        ctx[name].define_problem_bounding_box(setup)
        expected = {
            "domain_origin_x": -0.25,
            "domain_origin_y": -0.25,
            "domain_size_x": 1.3,
            "domain_size_y": 0.8,
        }
        for key, want in expected.items():
            got = setup[key]
            if abs(got - want) > 1e-12:
                return f"{name} {key}: expected far corner {want}, got {got}"
    return None


def check_right_wall_has_common_physical_position(ctx):
    """The right wall is level-independent and does not use a leaked rx value."""
    for name in ("uniform", "lts"):
        wall_positions = []
        for grid in _grids(ctx["dec"]):
            result = ctx[name]._wall_lattice(grid, _setup(), lambda _x, _y: True)
            box_x, _box_y, _gx, _gy, normal_x, _ny = result
            right_wall = {x for x, nx in zip(box_x, normal_x) if nx == -1.0}
            wall_positions.append(right_wall)
        if any(wall != {1.0} for wall in wall_positions):
            return f"{name} level-dependent right wall positions: {wall_positions}"
    return None


def check_wall_lattice_has_unique_positions(ctx):
    """Given corrected right-wall coordinates, corners add no duplicate sites."""
    for name in ("uniform", "lts"):
        for level, grid in enumerate(_grids(ctx["dec"])):
            result = ctx[name]._wall_lattice(grid, _setup(), lambda _x, _y: True)
            sites = list(zip(result[0], result[1]))
            duplicate_sites = sorted({site for site in sites if sites.count(site) > 1})
            if duplicate_sites:
                return f"{name} L{level} duplicate wall sites: {duplicate_sites}"
    return None


def check_dummy3d_reference_domain(ctx):
    """Given the legacy unit box, far-corner output preserves its grid."""
    setup = {"search_radius": 0.14, "box_x": 1.0, "box_y": 1.0, "box_z": 1.0, "dp": 0.04}
    ctx["dummy_3d"].define_problem_bounding_box(setup)
    for axis in "xyz":
        if abs(setup[f"domain_origin_{axis}"] + 0.35) > 1e-12:
            return f"3D dummy domain_origin_{axis} changed the reference origin"
        if abs(setup[f"domain_size_{axis}"] - 1.0) > 1e-12:
            return (f"3D dummy domain_size_{axis}: expected far corner 1.0, "
                    f"got {setup[f'domain_size_{axis}']}")

    bounds = ctx["dummy_3d"].resolve_refinement_bounds(
        {f"fine_{axis}_{side}": 0.25 if side == "min" else None
         for axis in "xyz" for side in ("min", "max")}, setup)
    for axis in "xyz":
        if abs(bounds[f"fine_{axis}_max"] - 1.0) > 1e-12:
            return f"3D dummy fine_{axis}_max is no longer centered on the unit box"

    zone = ctx["dec"].RectangularRefinementDef(
        fine_factor=0.5,
        fine_x_min=bounds["fine_x_min"], fine_x_max=bounds["fine_x_max"],
        fine_y_min=bounds["fine_y_min"], fine_y_max=bounds["fine_y_max"],
        fine_z_min=bounds["fine_z_min"], fine_z_max=bounds["fine_z_max"],
    )
    coarse, fine = ctx["dec"].build_rectangular_subdomain_grids(zone, setup)
    for axis in "xyz":
        if getattr(coarse, f"dims_{axis}") != 10:
            return f"3D dummy coarse dims_{axis}: expected 10, got {getattr(coarse, f'dims_{axis}')}"
        if getattr(coarse, f"phys_{axis}_min") > getattr(fine, f"phys_{axis}_min"):
            return f"3D dummy fine {axis} minimum lies outside coarse bounds"
        if getattr(coarse, f"phys_{axis}_max") < getattr(fine, f"phys_{axis}_max"):
            return f"3D dummy fine {axis} maximum lies outside coarse bounds"
    return None


def check_dam_break3d_domain_tiling(ctx):
    """Given the default tank, its 1:2 zones use far corners and nest."""
    setup = {
        "search_radius": 0.14,
        "box_length": 3.3,
        "box_height": 1.1,
        "box_width": 3.0,
        "dp": 0.04,
    }
    ctx["uniform_3d"].define_problem_bounding_box(setup)
    for axis, box in zip("xyz", (3.3, 1.1, 3.0)):
        expected_corner = box + 3 * setup["search_radius"]
        if abs(setup[f"domain_size_{axis}"] - expected_corner) > 1e-12:
            return (f"3D domain_size_{axis}: expected far corner {expected_corner}, "
                    f"got {setup[f'domain_size_{axis}']}")
        if abs(setup[f"domain_origin_{axis}"] + 2.5 * setup["search_radius"]) > 1e-12:
            return f"3D domain_origin_{axis} is not -2.5 search radii"

    zone = ctx["dec"].RectangularRefinementDef(
        fine_factor=0.5,
        fine_x_min=0.4, fine_x_max=1.2,
        fine_y_min=0.2, fine_y_max=0.82,
        fine_z_min=setup["domain_origin_z"], fine_z_max=0.5,
    )
    coarse, fine = ctx["dec"].build_rectangular_subdomain_grids(zone, setup)
    for axis, expected_coarse, expected_fine in zip("xyz", (30, 14, 27), (12, 10, 12)):
        if getattr(coarse, f"dims_{axis}") != expected_coarse:
            return f"3D coarse dims_{axis}: expected {expected_coarse}, got {getattr(coarse, f'dims_{axis}')}"
        if getattr(fine, f"dims_{axis}") != expected_fine:
            return f"3D fine dims_{axis}: expected {expected_fine}, got {getattr(fine, f'dims_{axis}')}"
        tolerance = 1e-9 * setup["search_radius"]
        if getattr(fine, f"phys_{axis}_min") < getattr(coarse, f"phys_{axis}_min") - tolerance:
            return f"3D fine {axis} minimum lies outside the coarse grid"
        if getattr(fine, f"phys_{axis}_max") > getattr(coarse, f"phys_{axis}_max") + tolerance:
            return f"3D fine {axis} maximum lies outside the coarse grid"
    return None
