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

def generate_fluid_particles(idx: int, g: SubdomainGrid, setup: dict) -> None:
    dp              = g.dp
    fluid_length    = setup["fluid_length"]
    fluid_height    = setup["fluid_height"]
    rho0            = setup["density"]
    c               = setup["speed_of_sound"]

    nx = round(fluid_length / dp)
    ny = round(fluid_height / dp)

    xs = dp * (np.arange(nx) + 1)          # candidate x positions
    ys = dp * (np.arange(ny) + 1)          # candidate y positions
    XX, YY = np.meshgrid(xs, ys, indexing="ij")
    XX, YY = XX.ravel(), YY.ravel()

    # Keep only particles inside this subdomain's physical box
    mask = (
        (XX > g.phys_x_min) & (XX <= g.phys_x_max) &
        (YY > g.phys_y_min) & (YY <= g.phys_y_max)
    )
    rx, ry = XX[mask], YY[mask]
    n      = len(rx)

    rho = np.array([ic.hydrostatic_density(y, fluid_height, rho0, c) for y in ry])
    r   = np.column_stack([rx, ry, np.zeros(n)])
    v   = np.zeros((n, 3))
    p   = np.zeros(n)
    pt  = np.zeros(n)

    cloud = saveParticlesVTK.create_pointcloud_polydata(r, v, rho, p, pt)
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_fluid.vtk")

    mass  = rho0 * dp**2
    print(f"[subdomain-{idx}] fluid particles: {n},  Epot0 = {mass * 9.81 * ry.sum():.4f} J")
    print(f"  physical box  x:[{g.phys_x_min:.4f}, {g.phys_x_max:.4f}]"
          f"  y:[{g.phys_y_min:.4f}, {g.phys_y_max:.4f}]")

    g.fluid_n = n


def generate_boundary_particles(idx: int, g: SubdomainGrid, setup: dict) -> None:
    dp   = g.dp
    rho0 = setup["density"]
    nl   = setup["n_boundary_layers"]

    box_length_n = round(setup["box_length"] / dp)
    box_height_n = round(setup["box_height"] / dp)

    box_rx, box_ry   = [], []
    ghost_rx, ghost_ry = [], []
    normal_x, normal_y = [], []

    def in_box(rx, ry):
        return (g.phys_x_min < rx <= g.phys_x_max and
                g.phys_y_min < ry <= g.phys_y_max)

    def add(bx, by, gx, gy, nx, ny):
        if in_box(bx, by):
            box_rx.append(bx);   box_ry.append(by)
            ghost_rx.append(gx); ghost_ry.append(gy)
            normal_x.append(nx); normal_y.append(ny)

    # Left wall
    for layer in range(nl):
        for j in range(box_height_n - 1):
            ry = (j + 1) * dp
            add(-layer * dp,          ry,
                 (layer + 1) * dp,    ry,
                 1.0, 0.0)

    # Bottom wall
    for layer in range(nl):
        for i in range(box_length_n - nl + 1):
            rx = (i + 1) * dp
            add(rx,  -layer * dp,
                rx,   (layer + 1) * dp,
                0.0,  1.0)

    x_last = box_rx[-1] + dp   # right-wall anchor

    # Right wall
    for layer in range(nl):
        for j in range(box_height_n - 1):
            ry = (j + 1) * dp
            add(x_last + layer * dp,          ry,
                x_last - (layer + 1) * dp,    ry,
                -1.0, 0.0)

    #: y_last = box_ry[-1] + dp

    #: # TODO: Top wall
    #: for layer in range(nl):
    #:     for i in range(box_length_n - nl + 1):
    #:         rx = (i + 1) * dp
    #:         add(rx,   y_last,
    #:             rx,   (layer + 1) * dp,
    #:             0.0,  -1.0)

    # Corners (90°)
    def corner(cx, cy, sx, sy):
        for layer in range(nl):
            for k in range(nl):
                bx = cx + k     * dp * sx
                by = cy + layer * dp * sy
                gx = cx - (k + 1)     * dp * sx
                gy = cy - (layer + 1) * dp * sy
                dr = np.array([gx - bx, gy - by])
                n  = dr / np.linalg.norm(dr)
                if in_box(bx, by):
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
    p       = np.zeros(n)
    pt      = np.ones(n)

    cloud = saveParticlesVTK.create_pointcloud_polydata(
        r, v, rho, p, pt,
        ghostNodes=ghosts, elementSize=elem_sz, normals=normals)
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_boundary.vtk")

    print(f"[subdomain-{idx}] boundary particles: {n}")
    g.boundary_n = n

def save_grid(grids: List[SubdomainGrid], setup: dict) -> None:
    import domainGrid
    h0 = setup["search_radius"]
    for i, g in enumerate(grids):
        ox = setup["domain_origin_x"] + h0 * g.origin_glob_x
        oy = setup["domain_origin_y"] + h0 * g.origin_glob_y
        n_cells = g.dims_x * g.dims_y
        domainGrid.domainGrid(
            g.dims_x, g.dims_y, 1,
            ox, oy, 0,
            np.zeros(n_cells),   # deprecated gridSector
            g.search_radius,
            f"sources/subdomain-{i}dambreak_grid.vtk"
        )


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

    if args.fluid_length > args.box_length or args.fluid_height > args.box_length:
        sys.exit("Error: fluid block must fit inside the tank.")

    for path in ("./results", "./sources"):
        os.makedirs(path, exist_ok=True)

    # Step 1: Based on the box size, compute domain bounding box
    define_problem_bounding_box(setup)

    # Step 2: Define subdomains (physical coords + refinement factor)
    subdomain_defs = [
        dec.SubdomainDef(factor=1.0, x_max=1.0), # coarse
        dec.SubdomainDef(factor=0.5, x_min=1.0), # fine
    ]

    # Step 3: Get parameters of subdomain grids
    grids = dec.build_subdomain_grids(subdomain_defs, setup)

    print("\nComplete setup:")
    pprint(setup)
    print("\nSubdomain grids:")
    for i, g in enumerate(grids):
        pprint({f"subdomain-{i}": g})

    # Step 4: Generate particles for each subdomain
    for i, g in enumerate(grids):
        generate_fluid_particles(i, g, setup)
        generate_boundary_particles(i, g, setup)

    # Step 5: Accumulate global particles counts
    setup["fluid_n"]    = sum(g.fluid_n    for g in grids)
    setup["boundary_n"] = sum(g.boundary_n for g in grids)

    # Step 6: Write outputs
    save_grid(grids, setup)
    cf.write_simulation_params(setup)
    dec.write_distributed_domain_params(grids, setup)
    cf.write_measuretool_params(setup)
