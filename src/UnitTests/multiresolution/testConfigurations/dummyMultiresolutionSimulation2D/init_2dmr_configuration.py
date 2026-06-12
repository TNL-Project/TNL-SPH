#! /usr/bin/env python3
"""
Multiresolution SPH test configuration initializer.

Supports multiple named configurations defined in configurations.py.
Each configuration specifies geometrical parameters only (domain size,
refinement region bounds, particle spacing).

Fine-region bounds set to None in a configuration are resolved to the
corresponding domain boundary after the bounding box is computed:
  fine_x_min=None -> domain_origin_x
  fine_x_max=None -> domain_origin_x + domain_size_x
  fine_y_min=None -> domain_origin_y
  fine_y_max=None -> domain_origin_y + domain_size_y
"""

import os
import sys
import argparse
import numpy as np
from typing import List
from pprint import pprint

_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.normpath(os.path.join(_script_dir, '../../../../tools')))
import saveParticlesVTK

import decomposition as dec
import initialCondifionFunctions as ic
import writeInitConfigFile as cf

from configurations import CONFIGURATIONS


def resolve_refinement_bounds(config: dict, setup: dict) -> dict:
    resolved = dict(config)
    boundary_map = {
        "fine_x_min": "domain_origin_x",
        "fine_x_max": ("domain_origin_x", "domain_size_x"),
        "fine_y_min": "domain_origin_y",
        "fine_y_max": ("domain_origin_y", "domain_size_y"),
    }
    for fine_key, domain_key in boundary_map.items():
        if resolved.get(fine_key) is None:
            if isinstance(domain_key, tuple):
                origin_key, size_key = domain_key
                resolved[fine_key] = setup[origin_key] + setup[size_key]
            else:
                resolved[fine_key] = setup[domain_key]
    return resolved


def save_grid(grids: List[dec.SubdomainGrid], setup: dict, output_dir: str) -> None:
    import domainGrid
    h0 = setup["search_radius"]
    for i, g in enumerate(grids):
        ox = setup["domain_origin_x"] + h0 * g.factor * g.origin_glob_x
        oy = setup["domain_origin_y"] + h0 * g.factor * g.origin_glob_y
        n_cells = g.dims_x * g.dims_y
        domainGrid.domainGrid(
            g.dims_x, g.dims_y, 1,
            ox, oy, 0,
            np.zeros(n_cells),
            g.search_radius,
            f"{output_dir}/subdomain-{i}dambreak_grid.vtk"
        )


def write_distributed_domain_params_rectangular(
        coarse_grid: dec.SubdomainGrid,
        fine_grid:   dec.SubdomainGrid,
        setup:       dict,
        output_dir:  str
) -> None:
    grids = [coarse_grid, fine_grid]
    h0    = setup["search_radius"]
    axes  = ["x", "y"]
    fact  = 2

    domain_origin = {ax: setup[f"domain_origin_{ax}"] for ax in axes}

    with open(f"{output_dir}/config-distributed-domain.ini", "w") as f:
        f.write("# Subdomain information\n\n")
        for i, g in enumerate(grids):
            prefix = f"subdomain-{i}-"

            params = {
                "fluid-particles":      f"{output_dir}/subdomain-{i}-dambreak_fluid.vtk",
                "boundary-particles":   f"{output_dir}/subdomain-{i}-dambreak_boundary.vtk",
                "fluid_n":              g.fluid_n,
                "boundary_n":           g.boundary_n,
                "fluid_n_allocated":    fact * g.fluid_n if fact * g.fluid_n > 0 else setup["fluid_n"],
                "boundary_n_allocated": fact * g.boundary_n if fact * g.boundary_n > 0 else setup["boundary_n"],
                "refinement-factor":    g.factor,
            }

            ax_groups = {
                "origin-global-coords": {},
                "grid-dimensions":      {},
                "origin":               {},
                "size":                 {},
            }
            for ax in axes:
                origin_glob = getattr(g, f"origin_glob_{ax}")
                dims        = getattr(g, f"dims_{ax}")
                origin_phys = domain_origin[ax] + g.search_radius * origin_glob
                size_phys   = g.search_radius * dims

                ax_groups["origin-global-coords"][ax] = origin_glob
                ax_groups["grid-dimensions"][ax]      = dims
                ax_groups["origin"][ax]               = f"{origin_phys:.7f}"
                ax_groups["size"][ax]                 = f"{size_phys:.7f}"

            for keyword, ax_values in ax_groups.items():
                for ax, value in ax_values.items():
                    params[f"{keyword}-{ax}"] = value

            for k, v in params.items():
                f.write(f"{prefix}{k} = {v}\n")
            f.write("\n")


def write_simulation_params(setup: dict, output_dir: str, config_name: str):
    with open(os.path.join(_script_dir, "dummyConfig2D_template.ini"), "r") as f:
        cfg = cf.safe_replace(f.read(), cf.ini_replacements, setup)
    cfg = cfg.replace("placeholderSubdomainsConfigPath", f"{output_dir}/config-distributed-domain.ini")
    with open(f"{output_dir}/dummyConfig2D.ini", "w") as f:
        f.write(cfg)


def build_setup(config: dict) -> dict:
    dp = config["dp"]
    h  = config["h_coef"] * dp
    return {
        "box_x":            config["box_x"],
        "box_y":            config["box_y"],
        "dp":               dp,
        "smoothing_length": h,
        "search_radius":    2.0 * h,
    }


def define_problem_bounding_box(setup: dict):
    sr = setup["search_radius"]
    domain = {
        "domain_origin_x": -2.5 * sr,
        "domain_origin_y": -2.5 * sr,
        "domain_size_x":   setup["box_x"] + 2.5 * sr,
        "domain_size_y":   setup["box_y"] + 2.5 * sr,
    }
    setup.update(domain)


def list_configurations():
    print("Available configurations:\n")
    for name, cfg in CONFIGURATIONS.items():
        desc = cfg.get("description", "(no description)")
        print(f"  {name:25s}  {desc}")
    print()


def parse_args():
    ap = argparse.ArgumentParser(description="Multiresolution SPH test configuration generator")
    ap.add_argument("--config-name", type=str, default="dummy-center",
                    choices=list(CONFIGURATIONS.keys()),
                    help="named configuration to initialize")
    ap.add_argument("--list", action="store_true",
                    help="list available configurations and exit")

    g = ap.add_argument_group("overrides")
    g.add_argument("--dp", type=float, default=None, help="override particle spacing")
    g.add_argument("--h-coef", type=float, default=None, help="override smoothing length coefficient")
    g.add_argument("--box-x", type=float, default=None, help="override domain x-size")
    g.add_argument("--box-y", type=float, default=None, help="override domain y-size")

    return ap.parse_args()


def init_configuration(config_name: str, overrides: dict):
    config = dict(CONFIGURATIONS[config_name])

    cli_map = {
        "dp":     "dp",
        "h_coef": "h_coef",
        "box_x":  "box_x",
        "box_y":  "box_y",
    }
    for cli_key, cfg_key in cli_map.items():
        val = overrides.get(cli_key)
        if val is not None:
            config[cfg_key] = val

    setup = build_setup(config)
    define_problem_bounding_box(setup)

    config = resolve_refinement_bounds(config, setup)

    output_dir = f"sources/{config_name}"
    results_dir = f"results/{config_name}"
    for path in (results_dir, output_dir):
        os.makedirs(path, exist_ok=True)

    rdef = dec.RectangularRefinementDef(
        fine_factor = config["fine_factor"],
        fine_x_min  = config["fine_x_min"],
        fine_x_max  = config["fine_x_max"],
        fine_y_min  = config["fine_y_min"],
        fine_y_max  = config["fine_y_max"],
    )
    print(rdef)

    coarse_grid, fine_grid = dec.build_rectangular_subdomain_grids(rdef, setup)
    print("\nCoarse grid:"); pprint(coarse_grid)
    print("\nFine grid:");   pprint(fine_grid)

    setup["fluid_n"]    = coarse_grid.fluid_n + fine_grid.fluid_n
    setup["boundary_n"] = coarse_grid.boundary_n

    save_grid([coarse_grid, fine_grid], setup, output_dir)
    write_simulation_params(setup, output_dir, config_name)
    write_distributed_domain_params_rectangular(coarse_grid, fine_grid, setup, output_dir)

    print(f"\nConfiguration '{config_name}' initialized.")
    print(f"  Output directory: {output_dir}/")


if __name__ == "__main__":
    args = parse_args()

    if args.list:
        list_configurations()
        sys.exit(0)

    overrides = {
        "dp":     args.dp,
        "h_coef": args.h_coef,
        "box_x":  args.box_x,
        "box_y":  args.box_y,
    }

    init_configuration(args.config_name, overrides)
