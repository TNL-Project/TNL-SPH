#! /usr/bin/env python3
"""Distributed (multi-GPU) SPH dam-break 3-D initial condition generator."""

import json
import os
import sys
import argparse
import subprocess
import numpy as np
from pprint import pprint

sys.path.append('../../../src/tools')
import saveParticlesVTK
import readVtkParticles
import writeInitConfigFile as cf
import computeDomainSize
import decomposition as dec
import domainGrid

def generate_geometry(dp: float) -> None:
   subprocess.check_call(
      ['./generateGeometryWithDualSPHysicsGenCase.sh', str(dp)],
      cwd='./template/generateGeometryWithDualSPHysicsGenCase/'
   )

def process_fluid_particles(setup: dict) -> np.ndarray:
   vtk_path = f'./sources/genCaseGeometries/dambreak_fluid_dp{setup["dp"]}.vtk'
   points, point_data = readVtkParticles.read_vtk_particles(vtk_path)

   fluid_n = len(points)
   cloud = saveParticlesVTK.create_pointcloud_polydata(
      points,
      point_data['Vel'].astype(float),
      point_data['Rhop'],
      np.zeros(fluid_n),
      np.zeros(fluid_n),
   )
   saveParticlesVTK.save_polydata(cloud, "sources/dambreak_fluid.vtk")
   setup["fluid_n"] = fluid_n
   return points

def process_boundary_particles(setup: dict) -> np.ndarray:
   vtk_path = f'./sources/genCaseGeometries/dambreak_bound_dp{setup["dp"]}.vtk'
   points, _ = readVtkParticles.read_vtk_particles(vtk_path)

   box_n = len(points)
   cloud = saveParticlesVTK.create_pointcloud_polydata(
      points,
      np.zeros((box_n, 3)),
      setup['density'] * np.ones(box_n),
      np.zeros(box_n),
      np.zeros(box_n),
   )
   saveParticlesVTK.save_polydata(cloud, "sources/dambreak_boundary.vtk")

   setup["boundary_n"] = box_n
   setup["domain_origin_x"] = float(np.min(points[:, 0]))
   setup["domain_origin_y"] = float(np.min(points[:, 1]))
   setup["domain_origin_z"] = float(np.min(points[:, 2]))
   setup["domain_end_x"] = float(np.max(points[:, 0]))
   setup["domain_end_y"] = float(np.max(points[:, 1]))
   setup["domain_end_z"] = float(np.max(points[:, 2]))
   return points

def write_subdomain_vtk(grids: list, splits: dict, setup: dict) -> None:
   density = setup["density"]

   for g in grids:
      sub_fluid_r, sub_boundary_r = splits[(g.ix, g.iy, g.iz)]

      fn = len(sub_fluid_r)
      cloud = saveParticlesVTK.create_pointcloud_polydata(
         sub_fluid_r,
         np.zeros((fn, 3)),
         density * np.ones(fn),
         np.zeros(fn),
         np.zeros(fn),
      )
      saveParticlesVTK.save_polydata(
         cloud, f"sources/dambreak_subdomain-x-{g.ix}-y-{g.iy}-fluid.vtk"
      )

      bn = len(sub_boundary_r)
      cloud = saveParticlesVTK.create_pointcloud_polydata(
         sub_boundary_r,
         np.zeros((bn, 3)),
         density * np.ones(bn),
         np.zeros(bn),
         np.ones(bn),
      )
      saveParticlesVTK.save_polydata(
         cloud, f"sources/dambreak_subdomain-x-{g.ix}-y-{g.iy}-boundary.vtk"
      )

      print(f"  [subdomain-x-{g.ix}-y-{g.iy}] fluid: {fn}, boundary: {bn}")

def write_simulation_params(setup: dict, grids: list) -> None:
   with open('template/config_template.jsonc', 'r') as f:
      cfg = f.read()

   cfg = cf.safe_replace(cfg, cf.ini_replacements, setup)

   dd_data = dec.build_distributed_domain_data(grids, setup, distributed=True, particles_filename_pattern="dambreak")
   dd_json = json.dumps(dd_data, indent=8)
   dd_json = dd_json.strip()[1:-1]
   cfg = cfg.replace('placeholderDistributedDomainContent', dd_json)

   with open('sources/config.jsonc', 'w') as f:
      f.write(cfg)

def parse_args():
   ap = argparse.ArgumentParser(
      description="Distributed (multi-GPU) SPH dam-break 3-D initial condition generator"
   )

   g = ap.add_argument_group("distribution parameters")
   g.add_argument("--subdomains-x", type=int, default=3, help="number of subdomains in x direction")
   g.add_argument("--subdomains-y", type=int, default=1, help="number of subdomains in y direction")
   g.add_argument("--overlap-width", type=int, default=1, help="width of domain overlap in cells")

   g = ap.add_argument_group("resolution parameters")
   g.add_argument("--dp", type=float, default=0.02, help="initial distance between particles")
   g.add_argument("--h-coef", type=float, default=2, help="smoothing length coefficient")

   g = ap.add_argument_group("simulation parameters")
   g.add_argument("--density", type=float, default=1000, help="referential density of the fluid")
   g.add_argument("--speed-of-sound", type=float, default=45.17, help="speed of sound")
   g.add_argument("--cfl", type=float, default=0.15, help="CFL condition number")

   g = ap.add_argument_group("control initialization")
   g.add_argument('--generate-geometry', default=True, action=argparse.BooleanOptionalAction,
                   help="generate new geometry with gencase")

   return ap.parse_args()

def build_setup(args) -> dict:
   dp = args.dp
   h  = args.h_coef * dp
   return {
      "dp":                   dp,
      "h_coef":               args.h_coef,
      "density":              args.density,
      "speed_of_sound":       args.speed_of_sound,
      "cfl":                  args.cfl,
      "particle_mass":        args.density * (dp ** 3),
      "smoothing_length":     h,
      "search_radius":        2 * h,
      "time_step":            args.cfl * h / args.speed_of_sound,
      "subdomains_x":         args.subdomains_x,
      "subdomains_y":         args.subdomains_y,
      "number_of_subdomains": args.subdomains_x * args.subdomains_y,
      "overlap_width":        args.overlap_width,
      "fluid_alloc_factor":   2,
      "boundary_alloc_factor": 3,
   }

if __name__ == "__main__":
   args  = parse_args()
   setup = build_setup(args)

   for path in ("./results", "./sources"):
      os.makedirs(path, exist_ok=True)

   if args.generate_geometry:
      generate_geometry(setup["dp"])

   fluid_r = process_fluid_particles(setup)
   box_r = process_boundary_particles(setup)

   setup["allocated_fluid_n"] = setup["fluid_n"]
   setup["allocated_boundary_n"] = setup["boundary_n"]

   computeDomainSize.compute_domain_size(setup)

   grids = dec.decompose_domain(setup, fluid_r)
   splits = dec.split_particles_to_subdomains(grids, fluid_r, box_r)
   write_subdomain_vtk(grids, splits, setup)

   print("\nComplete example setup:")
   pprint(setup)
   write_simulation_params(setup, grids)
   dec.write_distributed_domain_params(grids, setup, distributed=True, particles_filename_pattern="dambreak")

   domainGrid.write_domain_grid(setup, "sources/dambreak_grid.vtk")
