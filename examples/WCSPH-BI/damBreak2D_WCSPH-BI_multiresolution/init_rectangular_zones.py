#! /usr/bin/env python3
"""
Multiresolution SPH dam-break initializer.

Subdomains are defined as axis-aligned parts in physical coordinates.
Each subdomain carries its own refinement factor (dp_subdomain = factor * dx_L0).
"""

import math
import os
import sys
import argparse
import numpy as np
from dataclasses import dataclass
from typing import List
from pprint import pprint

sys.path.append('../../../src/tools')
import saveParticlesVTK
import init_generate_standard as single_resolution

#  -----------
import decomposition as dec
import initialCondifionFunctions as ic
import writeInitConfigFile as cf


# ---------------------------------------------------------------------------
# Particle generators
# ---------------------------------------------------------------------------

def generate_fluid_particles(
        idx:         int,
        g:           SubdomainGrid,
        setup:       dict,
        exclude_box: dict = None   # optional {"x_min", "x_max", "y_min", "y_max"}
) -> None:
    dp           = g.dp
    fluid_length = setup["fluid_length"]
    fluid_height = setup["fluid_height"]
    rho0         = setup["density"]
    c            = setup["speed_of_sound"]

    nx = round(fluid_length / dp)
    ny = round(fluid_height / dp)

    xs = dp * (np.arange(nx) + 1)
    ys = dp * (np.arange(ny) + 1)
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    XX, YY = XX.ravel(), YY.ravel()

    # Keep only particles inside this subdomain's physical box
    mask = (
        (XX > g.phys_x_min) & (XX <= g.phys_x_max) &
        (YY > g.phys_y_min) & (YY <= g.phys_y_max)
    )

    # Optionally exclude a rectangular region (e.g. fine region from coarse domain)
    if exclude_box is not None:
        in_excluded = (
            (XX > exclude_box["x_min"]) & (XX <= exclude_box["x_max"]) &
            (YY > exclude_box["y_min"]) & (YY <= exclude_box["y_max"])
        )
        mask = mask & ~in_excluded

    rx, ry = XX[mask], YY[mask]
    n      = len(rx)

    rho = np.array([ic.hydrostatic_density(y, fluid_height, rho0, c) for y in ry])
    r   = np.column_stack([rx, ry, np.zeros(n)])
    v, p, pt = np.zeros((n, 3)), np.zeros(n), np.zeros(n)

    cloud = saveParticlesVTK.create_pointcloud_polydata(r, v, rho, p, pt)
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_fluid.vtk")

    mass = rho0 * dp**2
    print(f"[subdomain-{idx}] fluid particles: {n},  Epot0 = {mass * 9.81 * ry.sum():.4f} J")
    print(f"  physical box  x:[{g.phys_x_min:.4f}, {g.phys_x_max:.4f}]"
          f"  y:[{g.phys_y_min:.4f}, {g.phys_y_max:.4f}]")
    g.fluid_n = n


def generate_boundary_particles(
        idx:         int,
        g:           SubdomainGrid,
        setup:       dict,
        exclude_box: dict = None
) -> None:
    dp   = g.dp
    rho0 = setup["density"]
    nl   = setup["n_boundary_layers"]

    box_length_n = round(setup["box_length"] / dp)
    box_height_n = round(setup["box_height"] / dp)

    box_rx, box_ry     = [], []
    ghost_rx, ghost_ry = [], []
    normal_x, normal_y = [], []

    def in_subdomain(rx, ry):
        return (g.phys_x_min < rx <= g.phys_x_max and
                g.phys_y_min < ry <= g.phys_y_max)

    def in_excluded(rx, ry):
        if exclude_box is None:
            return False
        return (exclude_box["x_min"] < rx <= exclude_box["x_max"] and
                exclude_box["y_min"] < ry <= exclude_box["y_max"])

    def add(bx, by, gx, gy, nx, ny):
        if in_subdomain(bx, by) and not in_excluded(bx, by):
            box_rx.append(bx);   box_ry.append(by)
            ghost_rx.append(gx); ghost_ry.append(gy)
            normal_x.append(nx); normal_y.append(ny)

    # Left wall
    for layer in range(nl):
        for j in range(box_height_n - 1):
            ry = (j + 1) * dp
            add(-layer * dp, ry, (layer + 1) * dp, ry, 1.0, 0.0)

    # Bottom wall
    for layer in range(nl):
        for i in range(box_length_n - nl + 1):
            rx = (i + 1) * dp
            add(rx, -layer * dp, rx, (layer + 1) * dp, 0.0, 1.0)

    #x_last = box_rx[-1] + dp if box_rx else dp * box_length_n
    x_last = rx + dp if box_rx else dp * box_length_n

    # Right wall
    for layer in range(nl):
        for j in range(box_height_n - 1):
            ry = (j + 1) * dp
            add(x_last + layer * dp, ry, x_last - (layer + 1) * dp, ry, -1.0, 0.0)

    # Corners
    def corner(cx, cy, sx, sy):
        for layer in range(nl):
            for k in range(nl):
                bx = cx + k     * dp * sx
                by = cy + layer * dp * sy
                gx = cx - (k + 1)     * dp * sx
                gy = cy - (layer + 1) * dp * sy
                dr = np.array([gx - bx, gy - by])
                n  = dr / np.linalg.norm(dr)
                if in_subdomain(bx, by) and not in_excluded(bx, by):
                    box_rx.append(bx);   box_ry.append(by)
                    ghost_rx.append(gx); ghost_ry.append(gy)
                    normal_x.append(n[0]); normal_y.append(n[1])

    corner(0.0,    0.0, -1, -1)
    corner(x_last, 0.0, +1, -1)

    n = len(box_rx)
    r       = np.column_stack([box_rx,   box_ry,   np.zeros(n)])
    ghosts  = np.column_stack([ghost_rx, ghost_ry, np.zeros(n)])
    normals = np.column_stack([normal_x, normal_y, np.zeros(n)])
    elem_sz = dp * np.ones(n)
    v       = np.zeros((n, 3))
    rho     = rho0 * np.ones(n)
    p, pt   = np.zeros(n), np.ones(n)

    cloud = saveParticlesVTK.create_pointcloud_polydata(
        r, v, rho, p, pt,
        ghostNodes=ghosts, elementSize=elem_sz, normals=normals)
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_boundary.vtk")

    print(f"[subdomain-{idx}] boundary particles: {n}")
    g.boundary_n = n

def _fine_box(fine_grid: SubdomainGrid) -> dict:
    return {
        "x_min": fine_grid.phys_x_min, "x_max": fine_grid.phys_x_max,
        "y_min": fine_grid.phys_y_min, "y_max": fine_grid.phys_y_max,
    }

def generate_fluid_particles_rectangular(
        coarse_grid: SubdomainGrid,
        fine_grid:   SubdomainGrid,
        setup:       dict
) -> None:
    generate_fluid_particles(1, fine_grid,   setup)
    generate_fluid_particles(0, coarse_grid, setup, exclude_box=_fine_box(fine_grid))

def generate_boundary_particles_rectangular(
        coarse_grid: SubdomainGrid,
        fine_grid:   SubdomainGrid,
        setup:       dict
) -> None:
    generate_boundary_particles(0, coarse_grid, setup, exclude_box=_fine_box(fine_grid))
    generate_boundary_particles(1, fine_grid, setup)

def save_grid(grids: List[SubdomainGrid], setup: dict) -> None:
    import domainGrid
    h0 = setup["search_radius"]
    for i, g in enumerate(grids):
        ox = setup["domain_origin_x"] + h0 * g.factor * g.origin_glob_x
        oy = setup["domain_origin_y"] + h0 * g.factor * g.origin_glob_y
        n_cells = g.dims_x * g.dims_y
        domainGrid.domainGrid(
            g.dims_x, g.dims_y, 1,
            ox, oy, 0,
            np.zeros(n_cells),   # deprecated gridSector
            g.search_radius,
            f"sources/subdomain-{i}dambreak_grid.vtk"
        )



def write_distributed_domain_params_rectangular(
        coarse_grid: SubdomainGrid,
        fine_grid:   SubdomainGrid,
        setup:       dict
) -> None:
    """
    Extended config writer for the rectangular nested refinement case.
    Writes the standard per-subdomain entries plus a [fine-region] block
    that MultiresolutionBoundary::initZonesRectangular reads to build
    the BufferSide array and frame zone.
    """
    grids = [coarse_grid, fine_grid]

    # Reuse the standard writer for per-subdomain grid entries
    dec.write_distributed_domain_params(grids, setup)

    # Append the fine-region block — read by initZonesRectangular
    is_3d = "domain_size_z" in setup
    axes  = ["x", "y", "z"] if is_3d else ["x", "y"]

    with open("sources/config-distributed-domain.ini", "a") as f:
        f.write("# Fine region specification (for rectangular MRB)\n")

        # Physical origin and size of the fine region
        for ax in axes:
            phys_min = getattr(fine_grid, f"phys_{ax}_min")
            phys_max = getattr(fine_grid, f"phys_{ax}_max")
            f.write(f"fine-region-origin-{ax} = {phys_min:.7f}\n")
            f.write(f"fine-region-size-{ax}   = {phys_max - phys_min:.7f}\n")

        # Which faces are active interfaces (not touching domain boundary)
        # A face is a domain boundary if phys_min/max coincides with domain origin/end
        sr = setup["search_radius"]
        for ax in axes:
            domain_min = setup[f"domain_origin_{ax}"]
            domain_max = setup[f"domain_origin_{ax}"] + setup[f"domain_size_{ax}"]
            phys_min   = getattr(fine_grid, f"phys_{ax}_min")
            phys_max   = getattr(fine_grid, f"phys_{ax}_max")

            # Face at min end: active if fine region does not touch domain boundary
            face_min_active = abs(phys_min - domain_min) > sr * 1e-3
            face_max_active = abs(phys_max - domain_max) > sr * 1e-3
            f.write(f"fine-region-face-min-{ax}-active = {int(face_min_active)}\n")
            f.write(f"fine-region-face-max-{ax}-active = {int(face_max_active)}\n")

        f.write("\n")

# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Multiresolution SPH dam-break setup generator")

    g = ap.add_argument_group("resolution")
    g.add_argument("--dp", type=float, default=0.002, help="coarse particle spacing")
    g.add_argument("--h-coef", type=float, default=1.75, help="smoothing length coefficient")

    g = ap.add_argument_group("domain geometry")
    g.add_argument("--box-length", type=float, default=1.61, help="tank length")
    g.add_argument("--box-height", type=float, default=2.0, help="tank height")
    g.add_argument("--fluid-length", type=float, default=0.6, help="initial fluid column length")
    g.add_argument("--fluid-height", type=float, default=0.3, help="initial fluid column height")

    g = ap.add_argument_group("physics")
    g.add_argument("--density", type=float, default=1000,  help="reference density [kg/m³]")
    g.add_argument("--speed-of-sound", type=float, default=34.3,  help="numerical speed of sound [m/s]")
    g.add_argument("--cfl", type=float, default=0.125) #FIXME: Use 0.125 instead of 0.25 due to missing subcycling
    g.add_argument("--alpha", type=float, default=0.02, help="artificial viscosity coefficient")
    g.add_argument("--dynamic-viscosity", type=float, default=0.001, help="dynamic viscosity [Pa·s]")

    g = ap.add_argument_group("formulation")
    g.add_argument("--bc-type", type=str, default="DBC")
    g.add_argument("--diffusive-term", type=str, default="MolteniDiffusiveTerm")
    g.add_argument("--viscous-term", type=str, default="ArtificialViscosity")

    return ap.parse_args()

def build_setup(args) -> dict:
    dp = args.dp
    h  = args.h_coef * dp
    return {
        "box_height":       args.box_height,
        "box_length":       args.box_length,
        "fluid_height":     args.fluid_height,
        "fluid_length":     args.fluid_length,
        "dp":               dp,
        "h_coef":           args.h_coef,
        "n_boundary_layers":1,
        "density":          args.density,
        "speed_of_sound":   args.speed_of_sound,
        "cfl":              args.cfl,
        "particle_mass":    args.density * dp**2,
        "smoothing_length": h,
        "search_radius":    2.0 * h,
        "time_step":        args.cfl * h / args.speed_of_sound,
        "alpha":            args.alpha,
        "dynamic_viscosity":args.dynamic_viscosity,
        "bc_type":          args.bc_type,
        "diffusive_term":   args.diffusive_term,
        "viscous_term":     args.viscous_term,
    }

def define_problem_bounding_box(setup):
    sr = setup["search_radius"]
    domain = {
        "domain_origin_x" : -2.5 * sr,
        "domain_origin_y" : -2.5 * sr,
        "domain_size_x" : setup["box_length"] + 2.5 * sr,
        "domain_size_y" : setup["box_height"] + 2.5 * sr,
    }
    setup.update( domain )

if __name__ == "__main__":
    args  = parse_args()
    setup = build_setup(args)

    for path in ("./results", "./sources"):
        os.makedirs(path, exist_ok=True)

    # Step 1: Domain bounding box
    define_problem_bounding_box(setup)

    # Step 2: Define rectangular refinement region
    #   fine region: x = [1.2, box_length], y = [0, 0.3]
    rdef = dec.RectangularRefinementDef(
        fine_factor = 0.5,
        fine_x_min  = 1.2,
        #fine_x_max  = args.box_length,   # touches right domain boundary → no MRB face there
        fine_x_max  = setup["domain_origin_x"] + setup["domain_size_x"] + 2.5 * setup["search_radius"],   # touches right domain boundary → no MRB face there
        fine_y_min  = setup["domain_origin_y"],  # touches bottom boundary → no MRB face there
        fine_y_max  = 0.3,
    )
    print(rdef)

    # Step 3: Build grids
    coarse_grid, fine_grid = dec.build_rectangular_subdomain_grids(rdef, setup)

    print("\nCoarse grid:"); pprint(coarse_grid)
    print("\nFine grid:");   pprint(fine_grid)

    # Step 4: Generate particles
    generate_fluid_particles_rectangular(coarse_grid, fine_grid, setup)
    generate_boundary_particles_rectangular(coarse_grid, fine_grid, setup)

    # Step 5: Global counts
    setup["fluid_n"]    = coarse_grid.fluid_n + fine_grid.fluid_n
    setup["boundary_n"] = coarse_grid.boundary_n

    # Step 6: Write outputs
    save_grid([coarse_grid, fine_grid], setup)
    cf.write_simulation_params(setup)
    write_distributed_domain_params_rectangular(coarse_grid, fine_grid, setup)
    cf.write_measuretool_params(setup)
