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


def _wall_lattice(g: dec.SubdomainGrid, setup: dict, keep_fn) -> tuple:
    """Tank wall lattice (4 walls + corners) at the subdomain's own
    resolution; keep_fn(bx, by) selects which box particles are kept."""
    dp   = g.dp
    nl   = setup["n_boundary_layers"]

    box_length_n = round(setup["box_length"] / dp)
    box_height_n = round(setup["box_height"] / dp)

    box_rx, box_ry     = [], []
    ghost_rx, ghost_ry = [], []
    normal_x, normal_y = [], []

    def add(bx, by, gx, gy, nx, ny):
        if keep_fn(bx, by):
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

    # Top wall
    y_top = box_height_n * dp
    for layer in range(nl):
        for i in range(box_length_n - nl + 1):
            rx = (i + 1) * dp
            add(rx, y_top + layer * dp, rx, y_top - (layer + 1) * dp, 0.0, -1.0)

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
                add(bx, by, gx, gy, n[0], n[1])

    corner(0.0,    0.0, -1, -1)
    corner(x_last, 0.0, +1, -1)
    corner(0.0,    y_top, -1, +1)
    corner(x_last, y_top, +1, +1)

    return box_rx, box_ry, ghost_rx, ghost_ry, normal_x, normal_y


def generate_boundary_particles(
        idx:         int,
        g:           dec.SubdomainGrid,
        setup:       dict,
        exclude_box: dict = None
) -> None:
    rho0 = setup["density"]

    def keep(bx, by):
        in_subdomain = (g.phys_x_min < bx <= g.phys_x_max and
                        g.phys_y_min < by <= g.phys_y_max)
        in_excluded = (exclude_box is not None and
                       exclude_box["x_min"] < bx <= exclude_box["x_max"] and
                       exclude_box["y_min"] < by <= exclude_box["y_max"])
        return in_subdomain and not in_excluded

    box_rx, box_ry, ghost_rx, ghost_ry, normal_x, normal_y = _wall_lattice(g, setup, keep)

    n = len(box_rx)
    r       = np.column_stack([box_rx,   box_ry,   np.zeros(n)])
    ghosts  = np.column_stack([ghost_rx, ghost_ry, np.zeros(n)])
    normals = np.column_stack([normal_x, normal_y, np.zeros(n)])
    elem_sz = g.dp * np.ones(n)
    v       = np.zeros((n, 3))
    rho     = rho0 * np.ones(n)
    p, pt   = np.zeros(n), np.ones(n)

    cloud = saveParticlesVTK.create_pointcloud_polydata(
        r, v, rho, p, pt,
        ghostNodes=ghosts, elementSize=elem_sz, normals=normals)
    saveParticlesVTK.save_polydata(cloud, f"sources/subdomain-{idx}-dambreak_boundary.vtk")

    print(f"[subdomain-{idx}] boundary particles: {n}")
    g.boundary_n = n


def _frame_and_band(own_g: dec.SubdomainGrid, nb_g: dec.SubdomainGrid, setup: dict) -> tuple:
    """Interface frame boxes in own cell units, mirroring
    MultiresolutionBoundary::initZones: the front frame is the own grid box
    for an outer interface, or the neighbor box mapped by the resolution
    ratio for an inner interface; the back frame expands/shrinks the front
    by one overlap cell accordingly."""
    sr_o = own_g.search_radius
    sr_n = nb_g.search_radius
    own_origin = (setup["domain_origin_x"] + own_g.origin_glob_x * sr_o,
                  setup["domain_origin_y"] + own_g.origin_glob_y * sr_o)
    nb_origin  = (setup["domain_origin_x"] + nb_g.origin_glob_x * sr_n,
                  setup["domain_origin_y"] + nb_g.origin_glob_y * sr_n)
    own_max = (own_origin[0] + own_g.dims_x * sr_o, own_origin[1] + own_g.dims_y * sr_o)
    nb_max  = (nb_origin[0]  + nb_g.dims_x * sr_n,  nb_origin[1]  + nb_g.dims_y * sr_n)

    inner = all(own_origin[d] <= nb_origin[d] <= own_max[d] for d in (0, 1))
    outer = all(nb_origin[d] <= own_origin[d] <= nb_max[d] for d in (0, 1))
    assert inner != outer

    resfac = sr_n / sr_o
    if inner:
        front_origin = (round(nb_g.origin_glob_x * resfac), round(nb_g.origin_glob_y * resfac))
        front_dims   = (round(nb_g.dims_x * resfac), round(nb_g.dims_y * resfac))
        orient = -1
    else:
        front_origin = (own_g.origin_glob_x, own_g.origin_glob_y)
        front_dims   = (own_g.dims_x, own_g.dims_y)
        orient = 1
    overlap = 1
    back_origin = tuple(front_origin[d] - overlap * orient for d in (0, 1))
    back_dims   = tuple(front_dims[d] + 2 * overlap * orient for d in (0, 1))
    return back_origin, back_dims, front_origin, front_dims


def generate_ghost_band_boundaries(grids: List[dec.SubdomainGrid], setup: dict) -> None:
    """Own-resolution wall lattice clipped to the interface band (frameBack
    XOR frameFront cell ring) for every directed interface, written as one
    VTK set per interface and consumed by the solver through the
    boundary-ghost-buffer-<p> entries.  The enumeration order must match
    the solver's own-major order over the chain links."""
    links = [(i, i + 1) for i in range(len(grids) - 1)]
    directed = [d for (a, b) in links for d in ((a, b), (b, a))]
    interfaces = [iface for own in range(len(grids)) for iface in directed if iface[0] == own]

    ox = setup["domain_origin_x"]
    oy = setup["domain_origin_y"]

    for p, (own, nb) in enumerate(interfaces):
        own_g, nb_g = grids[own], grids[nb]
        back_origin, back_dims, front_origin, front_dims = _frame_and_band(own_g, nb_g, setup)

        def in_band(bx, by, _bo=back_origin, _bd=back_dims,
                    _fo=front_origin, _fd=front_dims, _sr=own_g.search_radius):
            gc = (math.floor((bx - ox) / _sr), math.floor((by - oy) / _sr))
            in_back  = all(_bo[d] <= gc[d] < _bo[d] + _bd[d] for d in (0, 1))
            in_front = all(_fo[d] <= gc[d] < _fo[d] + _fd[d] for d in (0, 1))
            return in_back != in_front

        box_rx, box_ry, ghost_rx, ghost_ry, normal_x, normal_y = _wall_lattice(own_g, setup, in_band)

        n = len(box_rx)
        if n == 0:
            print(f"[interface {own}->{nb}] no boundary ghost particles in the band - solver will copy them")
            continue

        r       = np.column_stack([box_rx,   box_ry,   np.zeros(n)])
        ghosts  = np.column_stack([ghost_rx, ghost_ry, np.zeros(n)])
        normals = np.column_stack([normal_x, normal_y, np.zeros(n)])
        elem_sz = own_g.dp * np.ones(n)
        v       = np.zeros((n, 3))
        rho     = setup["density"] * np.ones(n)
        p_arr, pt = np.zeros(n), np.ones(n)

        cloud = saveParticlesVTK.create_pointcloud_polydata(
            r, v, rho, p_arr, pt,
            ghostNodes=ghosts, elementSize=elem_sz, normals=normals)
        saveParticlesVTK.save_polydata(cloud, f"sources/boundary_ghost_buffer_{p}.vtk")
        setup[f"ghost_buffer_{p}_n"] = n
        print(f"[interface {own}->{nb}] boundary ghost particles: {n}")

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
        subdomain_setup = {
            "search_radius": g.search_radius,
            "domain_size_x": g.dims_x * g.search_radius,
            "domain_size_y": g.dims_y * g.search_radius,
            "domain_origin_x": ox,
            "domain_origin_y": oy,
        }
        domainGrid.write_domain_grid(
            subdomain_setup,
            f"sources/subdomain-{i}-dambreak_grid.vtk"
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
    g.add_argument("--box-height", type=float, default=0.6, help="tank height")
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
        "fine_alloc_fact":  2,
    }

def define_problem_bounding_box(setup):
    sr = setup["search_radius"]
    domain = {
        "domain_origin_x" : -2.5 * sr,
        "domain_origin_y" : -2.5 * sr,
        "domain_size_x" : setup["box_length"] + 5.5 * sr,
        "domain_size_y" : setup["box_height"] + 5.5 * sr,
    }
    setup.update( domain )

if __name__ == "__main__":
    args  = parse_args()
    setup = build_setup(args)

    for path in ("./results", "./sources"):
        os.makedirs(path, exist_ok=True)

    # Step 1: Domain bounding box
    define_problem_bounding_box(setup)

    # Step 2: Define subdomains
    #   L0: whole tank (base resolution)
    #   L1: right-bottom zone, factor 0.5, touches right and bottom tank walls
    #   L2: right-bottom corner inside L1, factor 0.25, touches right and bottom tank walls
    domain_far_x = setup["domain_origin_x"] + setup["domain_size_x"] + 5.5 * setup["search_radius"]
    subdomain_defs = [
        dec.SubdomainDef(factor=1.0),
        dec.SubdomainDef(
            factor = 0.5,
            x_min  = 1.2,
            x_max  = domain_far_x,
            y_min  = setup["domain_origin_y"],
            y_max  = 0.3,
        ),
        dec.SubdomainDef(
            factor = 0.25,
            x_min  = 1.4,
            x_max  = domain_far_x,
            y_min  = setup["domain_origin_y"],
            y_max  = 0.15,
        ),
    ]

    # Step 3: Build grids
    coarse_grid, l1_grid, l2_grid = dec.build_subdomain_grids(subdomain_defs, setup)
    l2_grid.alloc_fact = 4  # the corner is fully submerged during the impacting slosh

    print("\nCoarse grid:"); pprint(coarse_grid)
    print("\nL1 grid:");      pprint(l1_grid)
    print("\nL2 grid:");      pprint(l2_grid)

    # Step 4: Generate particles
    generate_fluid_particles(1, l1_grid, setup)
    generate_fluid_particles(2, l2_grid, setup)
    generate_fluid_particles(0, coarse_grid, setup, exclude_box=_fine_box(l1_grid))

    generate_boundary_particles(0, coarse_grid, setup, exclude_box=_fine_box(l1_grid))
    generate_boundary_particles(1, l1_grid, setup, exclude_box=_fine_box(l2_grid))
    generate_boundary_particles(2, l2_grid, setup)

    generate_ghost_band_boundaries([coarse_grid, l1_grid, l2_grid], setup)

    # Step 5: Global counts and timestep driven by the finest resolution
    setup["fluid_n"] = coarse_grid.fluid_n + l1_grid.fluid_n + l2_grid.fluid_n
    setup["boundary_n"] = coarse_grid.boundary_n
    setup["number_of_subdomains"] = len(subdomain_defs)
    setup["time_step"] = setup["cfl"] * setup["smoothing_length"] * l2_grid.factor / setup["speed_of_sound"]

    # Step 6: Write outputs
    save_grid([coarse_grid, l1_grid, l2_grid], setup)
    cf.write_simulation_params(setup)
    dec.write_distributed_domain_params([coarse_grid, l1_grid, l2_grid], setup, particles_filename_pattern="dambreak")
    cf.write_measuretool_params(setup)
