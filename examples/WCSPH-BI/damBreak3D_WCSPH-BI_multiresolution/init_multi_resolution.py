#! /usr/bin/env python3
"""
Multiresolution SPH dam-break initializer — 3-D version.

Fluid particles are generated from scratch (meshgrid).
Boundary particles are generated via DualSPHysics genCase (called once per
resolution level) and then filtered into subdomains by physical coordinates.

Subdomains
----------
subdomain-0 : coarse  (dp_coarse = dp)
subdomain-1 : fine    (dp_fine   = fine_factor * dp,  fine_factor < 1)

CLI-specified fine region  --fine-{x,y,z}-{min,max}
"""

import math
import os
import sys
import argparse
import subprocess
import numpy as np
from pprint import pprint

sys.path.append('../../../src/tools')
import saveParticlesVTK
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from pathlib import Path

import decomposition as dec
import initialCondifionFunctions as ic
import writeInitConfigFile as cf

# initialize directories
example_dir = Path(__file__).parent
project_dir = (example_dir / ".." / ".." / ".." ).resolve()

# ---------------------------------------------------------------------------
# genCase helpers
# ---------------------------------------------------------------------------

def run_gencase(dp: float) -> None:
    """Call the geometry-generation shell script for the given particle spacing."""
    subprocess.check_call(
        ['./generateGeometryWithDualSPHysicsGenCase.sh', str(dp)],
        cwd='./template/generateGeometryWithDualSPHysicsGenCase/'
    )


def read_gencase_vtk(vtk_path: str):
    """Read a genCase VTK file and return the wrapped polydata object."""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_path)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    return dsa.WrapDataObject(reader.GetOutput())

# ---------------------------------------------------------------------------
# Fluid particle generation  (meshgrid — no genCase dependency)
# ---------------------------------------------------------------------------

def _box_mask_3d(pts, x_min, x_max, y_min, y_max, z_min, z_max):
    """Boolean mask: points strictly inside the given axis-aligned box."""
    return (
        (pts[:, 0] > x_min) & (pts[:, 0] <= x_max) &
        (pts[:, 1] > y_min) & (pts[:, 1] <= y_max) &
        (pts[:, 2] > z_min) & (pts[:, 2] <= z_max)
    )

def generate_fluid_particles_for_subdomain(
        idx:         int,
        g,                       # SubdomainGrid (from dec)
        vtk_path:    str,        # genCase boundary VTK at g.dp resolution
        setup:       dict,
        exclude_box: dict = None  # {"x_min","x_max","y_min","y_max","z_min","z_max"}
) -> None:

    rho0 = setup["density"]
    dp = setup["dp"]
    """
    Read a genCase boundary VTK (at the correct dp for this subdomain),
    keep only the particles that fall inside the subdomain physical box
    (and outside the exclude_box), then write the filtered cloud.
    """
    wpd = read_gencase_vtk(vtk_path)
    pts = np.array(wpd.Points, dtype=float)

    mask = _box_mask_3d(
        pts,
        g.phys_x_min, g.phys_x_max,
        g.phys_y_min, g.phys_y_max,
        g.phys_z_min, g.phys_z_max,
    )

    if exclude_box is not None:
        in_excl = _box_mask_3d(
            pts,
            exclude_box["x_min"], exclude_box["x_max"],
            exclude_box["y_min"], exclude_box["y_max"],
            exclude_box["z_min"], exclude_box["z_max"],
        )
        mask = mask & ~in_excl

    r       = pts[mask]
    n       = len(r)
    v       = np.zeros((n, 3))
    rho     = rho0 * np.ones(n)
    p       = np.zeros(n)
    pt      = np.ones(n)


    cloud = saveParticlesVTK.create_pointcloud_polydata(r, v, rho, p, pt)
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_fluid.vtk")

    mass = rho0 * dp**3
    print(f"[subdomain-{idx}] fluid particles : {n}")
    print(f"  dp = {dp:.5f}  Epot0 = {mass * 9.81 * r[ :, 2 ].sum():.4f} J")
    print(f"  box  x:[{g.phys_x_min:.4f},{g.phys_x_max:.4f}]"
          f"  y:[{g.phys_y_min:.4f},{g.phys_y_max:.4f}]"
          f"  z:[{g.phys_z_min:.4f},{g.phys_z_max:.4f}]")
    g.fluid_n = n


# ---------------------------------------------------------------------------
# Boundary particle generation  (genCase → filter by subdomain)
# ---------------------------------------------------------------------------


def process_boundary_particles_for_subdomain(
        idx:         int,
        g,                       # SubdomainGrid
        vtk_path:    str,        # genCase boundary VTK at g.dp resolution
        setup:       dict,
        exclude_box: dict = None
) -> None:
    """
    Read a genCase boundary VTK (at the correct dp for this subdomain),
    keep only the particles that fall inside the subdomain physical box
    (and outside the exclude_box), then write the filtered cloud.
    """
    wpd = read_gencase_vtk(vtk_path)
    pts = np.array(wpd.Points, dtype=float)

    mask = _box_mask_3d(
        pts,
        g.phys_x_min, g.phys_x_max,
        g.phys_y_min, g.phys_y_max,
        g.phys_z_min, g.phys_z_max,
    )

    if exclude_box is not None:
        in_excl = _box_mask_3d(
            pts,
            exclude_box["x_min"], exclude_box["x_max"],
            exclude_box["y_min"], exclude_box["y_max"],
            exclude_box["z_min"], exclude_box["z_max"],
        )
        mask = mask & ~in_excl

    r       = pts[mask]
    n       = len(r)
    v       = np.zeros((n, 3))
    rho     = setup["density"] * np.ones(n)
    p       = np.zeros(n)
    pt      = np.ones(n)
    elem_sz = (g.dp ** 2) * np.ones(n)

    # Normals — genCase boundary VTK should carry a 'Normal' array
    raw_normals = np.array(wpd.PointData["Normal"], dtype=float)[mask]
    zero_count  = 0
    normals     = np.zeros_like(raw_normals)
    for i, nrm in enumerate(raw_normals):
        mag = np.linalg.norm(nrm)
        if mag < 1e-12:
            print(f"  [subdomain-{idx}] WARNING: zero normal at particle {i}, pos {r[i]}")
            zero_count += 1
        else:
            normals[i] = nrm / mag
    if zero_count:
        print(f"  [subdomain-{idx}] {zero_count} undefined normals.")

    cloud = saveParticlesVTK.create_pointcloud_polydata(
        r, v, rho, p, pt,
        elementSize=elem_sz, normals=normals
    )
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_boundary.vtk")

    print(f"[subdomain-{idx}] boundary particles: {n}")
    g.boundary_n = n


# ---------------------------------------------------------------------------
# Top-level orchestration helpers
# ---------------------------------------------------------------------------

def _fine_box(fine_grid) -> dict:
    return {
        "x_min": fine_grid.phys_x_min, "x_max": fine_grid.phys_x_max,
        "y_min": fine_grid.phys_y_min, "y_max": fine_grid.phys_y_max,
        "z_min": fine_grid.phys_z_min, "z_max": fine_grid.phys_z_max,
    }


def generate_fluid_particles_rectangular(
        coarse_grid, fine_grid, setup: dict,
        coarse_vtk: str, fine_vtk: str
) -> None:
    """
    coarse_vtk : genCase boundary VTK at dp_coarse resolution
    fine_vtk   : genCase boundary VTK at dp_fine   resolution
    """
    generate_fluid_particles_for_subdomain(
            0, coarse_grid, coarse_vtk, setup,
            exclude_box=_fine_box(fine_grid)
    )

    generate_fluid_particles_for_subdomain(
            1, fine_grid, fine_vtk, setup
    )

def generate_boundary_particles_rectangular(
        coarse_grid, fine_grid, setup: dict,
        coarse_vtk: str, fine_vtk: str
) -> None:
    """
    coarse_vtk : genCase boundary VTK at dp_coarse resolution
    fine_vtk   : genCase boundary VTK at dp_fine   resolution
    """
    process_boundary_particles_for_subdomain(
        0, coarse_grid, coarse_vtk, setup,
        exclude_box=_fine_box(fine_grid)
    )
    process_boundary_particles_for_subdomain(
        1, fine_grid, fine_vtk, setup
    )


# ---------------------------------------------------------------------------
# Domain bounding box
# ---------------------------------------------------------------------------

def define_problem_bounding_box(setup: dict) -> None:
    sr = setup["search_radius"]
    setup.update({
        "domain_origin_x": -2.5 * sr,
        "domain_origin_y": -2.5 * sr,
        "domain_origin_z": -2.5 * sr,
        "domain_size_x":   setup["box_length"] + 5.0 * sr,
        "domain_size_y":   setup["box_height"] + 5.0 * sr,
        "domain_size_z":   setup["box_width"]  + 5.0 * sr,
    })


# ---------------------------------------------------------------------------
# Extended config writer (3-D, rectangular MRB)
# ---------------------------------------------------------------------------

def write_distributed_domain_params_rectangular(
        coarse_grid, fine_grid, setup: dict
) -> None:
    """
    Writes standard per-subdomain grid entries then appends the [fine-region]
    block (origin, size, active faces) that MultiresolutionBoundary::
    initZonesRectangular reads.
    """
    grids = [coarse_grid, fine_grid]
    dec.write_distributed_domain_params(grids, setup)

    axes = ["x", "y", "z"]
    sr   = setup["search_radius"]

    with open("sources/config-distributed-domain.ini", "a") as f:
        f.write("# Fine region specification (for rectangular MRB)\n")

        for ax in axes:
            phys_min = getattr(fine_grid, f"phys_{ax}_min")
            phys_max = getattr(fine_grid, f"phys_{ax}_max")
            f.write(f"fine-region-origin-{ax} = {phys_min:.7f}\n")
            f.write(f"fine-region-size-{ax}   = {phys_max - phys_min:.7f}\n")

        for ax in axes:
            domain_min = setup[f"domain_origin_{ax}"]
            domain_max = domain_min + setup[f"domain_size_{ax}"]
            phys_min   = getattr(fine_grid, f"phys_{ax}_min")
            phys_max   = getattr(fine_grid, f"phys_{ax}_max")

            face_min_active = abs(phys_min - domain_min) > sr * 1e-3
            face_max_active = abs(phys_max - domain_max) > sr * 1e-3
            f.write(f"fine-region-face-min-{ax}-active = {int(face_min_active)}\n")
            f.write(f"fine-region-face-max-{ax}-active = {int(face_max_active)}\n")

        f.write("\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Multiresolution SPH 3-D dam-break setup generator"
    )

    g = ap.add_argument_group("resolution")
    g.add_argument("--dp",      type=float, default=0.04,  help="coarse particle spacing")
    g.add_argument("--h-coef",  type=float, default=1.75,  help="smoothing length coefficient")

    g = ap.add_argument_group("domain geometry  (SPHERIC test case 2 defaults)")
    g.add_argument("--box-length", type=float, default=3.3,  help="tank length (x)")
    g.add_argument("--box-height", type=float, default=1.1,  help="tank height (y)")
    g.add_argument("--box-width",  type=float, default=3.0, help="tank width  (z)")
    g.add_argument("--fluid-length", type=float, default=0.40,  help="fluid column length (x)")
    g.add_argument("--fluid-height", type=float, default=0.30,  help="fluid column height (y)")
    g.add_argument("--fluid-width",  type=float, default=0.161, help="fluid column width  (z, usually = box-width)")

    g = ap.add_argument_group("fine zone  (physical coords)")
    g.add_argument("--fine-factor", type=float, default=0.5, help="dp_fine = fine_factor * dp  (< 1 → refinement)")
    g.add_argument("--fine-x-min", type=float, default=0.4)
    g.add_argument("--fine-x-max", type=float, default=1.2)
    g.add_argument("--fine-y-min", type=float, default=0.2)
    g.add_argument("--fine-y-max", type=float, default=0.82)
    g.add_argument("--fine-z-min", type=float, default=0.0)
    g.add_argument("--fine-z-max", type=float, default=0.5, help="Set to box-width to span the full tank width")

    g = ap.add_argument_group("physics")
    g.add_argument("--density",          type=float, default=1000.0)
    g.add_argument("--speed-of-sound",   type=float, default=45.17)
    g.add_argument("--cfl",              type=float, default=0.125)
    g.add_argument("--alpha",            type=float, default=0.02)
    g.add_argument("--dynamic-viscosity",type=float, default=0.001)

    g = ap.add_argument_group("formulation")
    g.add_argument("--bc-type",          type=str, default="BIConservative_numeric")
    g.add_argument("--diffusive-term",   type=str, default="MolteniDiffusiveTerm")
    g.add_argument("--viscous-term",     type=str, default="ArtificialViscosity")
    g.add_argument("--time-integration", type=str, default="VerletScheme", help="time integration scheme")

    g = ap.add_argument_group("genCase paths  (relative to this script)")
    g.add_argument("--coarse-fluid-vtk", type=str, default=None,
                   help="Pre-existing genCase fluid VTK at dp_coarse. "
                        "If omitted, genCase is called automatically.")
    g.add_argument("--fine-fluid-vtk",   type=str, default=None,
                   help="Pre-existing genCase fluid VTK at dp_fine. "
                        "If omitted, genCase is called automatically.")
    g.add_argument("--coarse-bound-vtk", type=str, default=None,
                   help="Pre-existing genCase boundary VTK at dp_coarse. "
                        "If omitted, genCase is called automatically.")
    g.add_argument("--fine-bound-vtk",   type=str, default=None,
                   help="Pre-existing genCase boundary VTK at dp_fine. "
                        "If omitted, genCase is called automatically.")
    g.add_argument("--no-gencase", action="store_true",
                   help="Skip genCase calls entirely (both VTK paths must be supplied).")

    return ap.parse_args()


def build_setup(args) -> dict:
    dp = args.dp
    h  = args.h_coef * dp
    return {
        # geometry
        "box_length":            args.box_length,
        "box_height":            args.box_height,
        "box_width":             args.box_width,
        "fluid_length":          args.fluid_length,
        "fluid_height":          args.fluid_height,
        "fluid_width":           args.fluid_width,
        # resolution
        "dp":                    dp,
        "h_coef":                args.h_coef,
        "n_boundary_layers":     1,
        "boundary_element_size": dp**2,
        # physics
        "density":               args.density,
        "speed_of_sound":        args.speed_of_sound,
        "cfl":                   args.cfl,
        "particle_mass":         args.density * dp**3,
        "smoothing_length":      h,
        "search_radius":         2.0 * h,
        "time_step":             args.cfl * h / args.speed_of_sound,
        "alpha":                 args.alpha,
        "dynamic_viscosity":     args.dynamic_viscosity,
        # formulation
        "bc_type":               args.bc_type,
        "diffusive_term":        args.diffusive_term,
        "viscous_term":          args.viscous_term,
        "time_integration" :     args.time_integration,
    }


def set_paraview_states_paths():
    import re
    template_dir = ( example_dir / "template" ).resolve()
    for pvsmfile in sorted(template_dir.iterdir()):
        if not pvsmfile.is_file() or ".pvsm" not in pvsmfile.name:
            continue
        print(pvsmfile)
        content = pvsmfile.read_text()
        content = content.replace("{project_dir}", str(project_dir))
        pvsmfile.write_text(content)

# ---------------------------------------------------------------------------
# Derived genCase VTK paths  (mirrors your single-res 3D naming convention)
# ---------------------------------------------------------------------------

def gencase_boundary_vtk(dp: float) -> str:
    return f"./sources/genCaseGeometries/dambreak_bound_dp{dp}.vtk"

def gencase_fluid_vtk(dp: float) -> str:
    return f"./sources/genCaseGeometries/dambreak_fluid_dp{dp}.vtk"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    args  = parse_args()
    setup = build_setup(args)

    dp_coarse = setup["dp"]
    dp_fine   = args.fine_factor * dp_coarse

    for path in ("./results", "./sources"):
        os.makedirs(path, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: Domain bounding box
    # ------------------------------------------------------------------
    define_problem_bounding_box(setup)

    # ------------------------------------------------------------------
    # Step 2: Fine region definition
    # ------------------------------------------------------------------
    rdef = dec.RectangularRefinementDef(
        fine_factor = args.fine_factor,
        fine_x_min  = args.fine_x_min,
        fine_x_max  = args.fine_x_max,
        fine_y_min  = args.fine_y_min,
        fine_y_max  = args.fine_y_max,
        fine_z_min  = setup["domain_origin_z"],
        fine_z_max  = args.fine_z_max,
    )
    print(rdef)

    # ------------------------------------------------------------------
    # Step 3: Build subdomain grids
    # ------------------------------------------------------------------
    coarse_grid, fine_grid = dec.build_rectangular_subdomain_grids(rdef, setup)

    print("\nCoarse grid:"); pprint(vars(coarse_grid))
    print("\nFine grid:");   pprint(vars(fine_grid))

    # ------------------------------------------------------------------
    # Step 4: Run genCase for boundaries (if needed)
    # ------------------------------------------------------------------
    coarse_fluid_vtk = args.coarse_fluid_vtk or gencase_fluid_vtk(dp_coarse)
    fine_fluid_vtk   = args.fine_fluid_vtk   or gencase_fluid_vtk(dp_fine)

    coarse_vtk = args.coarse_bound_vtk or gencase_boundary_vtk(dp_coarse)
    fine_vtk   = args.fine_bound_vtk   or gencase_boundary_vtk(dp_fine)

    if not args.no_gencase:
        if not os.path.exists(coarse_vtk):
            print(f"\n--- Running genCase at dp_coarse = {dp_coarse} ---")
            run_gencase(dp_coarse)
        else:
            print(f"[genCase] Using existing coarse boundary VTK: {coarse_vtk}")

        if not os.path.exists(fine_vtk):
            print(f"\n--- Running genCase at dp_fine = {dp_fine} ---")
            run_gencase(dp_fine)
        else:
            print(f"[genCase] Using existing fine boundary VTK:   {fine_vtk}")
    else:
        # Sanity check: files must exist when --no-gencase is set
        for path, label in [(coarse_vtk, "coarse"), (fine_vtk, "fine")]:
            if not os.path.exists(path):
                raise FileNotFoundError(
                    f"--no-gencase set but {label} boundary VTK not found: {path}"
                )

    # ------------------------------------------------------------------
    # Step 5: Generate fluid particles (meshgrid, from scratch)
    # ------------------------------------------------------------------
    print("\n--- Generating fluid particles ---")
    generate_fluid_particles_rectangular(
            coarse_grid, fine_grid, setup,
            coarse_vtk=coarse_fluid_vtk,
            fine_vtk=fine_fluid_vtk
    )

    # ------------------------------------------------------------------
    # Step 6: Generate boundary particles (genCase VTK → filter)
    # ------------------------------------------------------------------
    print("\n--- Generating boundary particles ---")
    generate_boundary_particles_rectangular(
        coarse_grid, fine_grid, setup,
        coarse_vtk=coarse_vtk,
        fine_vtk=fine_vtk,
    )

    # ------------------------------------------------------------------
    # Step 7: Global particle counts
    # ------------------------------------------------------------------
    setup["fluid_n"]    = coarse_grid.fluid_n + fine_grid.fluid_n
    setup["boundary_n"] = coarse_grid.boundary_n + fine_grid.boundary_n

    print(f"\nTotal fluid    particles : {setup['fluid_n']}")
    print(f"Total boundary particles : {setup['boundary_n']}")
    pprint(setup)

    # ------------------------------------------------------------------
    # Step 8: Write config files
    # ------------------------------------------------------------------
    cf.write_simulation_params(setup)
    write_distributed_domain_params_rectangular(coarse_grid, fine_grid, setup)
    cf.write_measuretool_params(setup)
    set_paraview_states_paths()
