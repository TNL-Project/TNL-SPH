import json
import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

@dataclass
class SubdomainDef:
   factor: float
   x_min:  float = -math.inf;  x_max: float = math.inf
   y_min:  float = -math.inf;  y_max: float = math.inf
   z_min:  float = -math.inf;  z_max: float = math.inf

@dataclass
class SubdomainGrid:
   factor:        float
   search_radius: float
   dp:            float
   phys_x_min:    float
   phys_x_max:    float
   phys_y_min:    float
   phys_y_max:    float
   phys_z_min:    float = None
   phys_z_max:    float = None
   origin_glob_x: int   = 0
   origin_glob_y: int   = 0
   origin_glob_z: int   = None
   dims_x:        int   = 0
   dims_y:        int   = 0
   dims_z:        int   = None
   fluid_n:       int   = 0
   boundary_n:    int   = 0
   index:         int   = 0
   ix:            int   = 0
   iy:            int   = 0
   iz:            int   = 0

@dataclass
class RectangularRefinementDef:
   fine_factor: float
   fine_x_min:  float
   fine_x_max:  float
   fine_y_min:  float
   fine_y_max:  float
   fine_z_min:  float = -math.inf
   fine_z_max:  float =  math.inf


# --- Multiresolution decomposition ---

def build_subdomain_grids(defs: List[SubdomainDef], setup: dict) -> List[SubdomainGrid]:
   h0     = setup["search_radius"]
   is_3d  = "domain_size_z" in setup
   axes   = ["x", "y", "z"] if is_3d else ["x", "y"]

   domain_origin = {ax: setup[f"domain_origin_{ax}"] for ax in axes}
   domain_size   = {ax: setup[f"domain_size_{ax}"]   for ax in axes}
   n_cells_glob  = {ax: int(math.ceil(domain_size[ax] / h0)) for ax in axes}

   def resolve_axis(sd: SubdomainDef, ax: str, f: float) -> tuple:
      fl   = int(round(1.0 / f))
      o    = domain_origin[ax]
      n_g  = n_cells_glob[ax]

      p0   = max(getattr(sd, f"{ax}_min"), o)
      og   = max(0, int(math.floor((p0 - o) / h0)))
      p0   = o + og * h0

      p1   = min(getattr(sd, f"{ax}_max"), o + n_g * h0)
      eg   = min(n_g, int(math.floor((p1 - o) / h0)))
      dims = (eg - og) * fl

      h    = f * h0
      return og * fl, dims, p0, p0 + dims * h

   grids: List[SubdomainGrid] = []
   for i, sd in enumerate(defs):
      f  = sd.factor
      h  = f * h0
      dp = f * setup["dp"]

      resolved = {ax: resolve_axis(sd, ax, f) for ax in axes}

      grids.append(SubdomainGrid(
         factor        = f,
         search_radius = h,
         dp            = dp,
         phys_x_min    = resolved["x"][2],
         phys_x_max    = resolved["x"][3],
         phys_y_min    = resolved["y"][2],
         phys_y_max    = resolved["y"][3],
         phys_z_min    = resolved["z"][2] if is_3d else None,
         phys_z_max    = resolved["z"][3] if is_3d else None,
         origin_glob_x = resolved["x"][0],
         origin_glob_y = resolved["y"][0],
         origin_glob_z = resolved["z"][0] if is_3d else None,
         dims_x        = resolved["x"][1],
         dims_y        = resolved["y"][1],
         dims_z        = resolved["z"][1] if is_3d else None,
         index         = i,
      ))

   return grids

def build_rectangular_subdomain_grids(
      rdef:  RectangularRefinementDef,
      setup: dict
) -> tuple:
   coarse_def = SubdomainDef(factor=1.0)
   fine_def   = SubdomainDef(
      factor = rdef.fine_factor,
      x_min  = rdef.fine_x_min,  x_max = rdef.fine_x_max,
      y_min  = rdef.fine_y_min,  y_max = rdef.fine_y_max,
      z_min  = rdef.fine_z_min,  z_max = rdef.fine_z_max,
   )
   coarse_grid, fine_grid = build_subdomain_grids([coarse_def, fine_def], setup)
   return coarse_grid, fine_grid


# --- Distributed (multi-GPU) decomposition ---
#
# Uniform-resolution 3-D grid decomposition (subdomains_x x subdomains_y x subdomains_z).
# For this particular case subdomains_z = 1, so z is not decomposed.
# Splits are computed from the fluid particle distribution (load-balanced).

def compute_grid_splits(
   fluid_r:           np.ndarray,
   fluid_n:           int,
   n_subdomains:      int,
   n_subdomains_axis: int,
   axis:              int,
   search_radius:     float,
) -> list:
   ptcs_per_sub = int(math.ceil(fluid_n / n_subdomains))
   splits = []
   for i in range(n_subdomains_axis - 1):
      idx = ptcs_per_sub * (i + 1)
      splits.append(math.ceil(fluid_r[idx, axis] / search_radius))
   return splits

def decompose_axis(
   n_subdomains:  int,
   grid_splits:   list,
   domain_origin: float,
   domain_size:   float,
   search_radius: float,
   overlap_width: int,
) -> Tuple[list, list, list, list]:
   n_cells_total = int(math.ceil(domain_size / search_radius))

   if n_subdomains == 1:
      return (
         [domain_origin],
         [domain_size],
         [n_cells_total],
         [overlap_width],
      )

   origins = []
   sizes = []
   grid_dims = []
   grid_index_origins = []

   for i in range(n_subdomains):
      if i == 0:
         dims       = grid_splits[0]
         origin     = domain_origin
         idx_origin = overlap_width
         size       = grid_splits[0] * search_radius
      elif i == n_subdomains - 1:
         dims       = n_cells_total - grid_splits[i - 1]
         origin     = domain_origin + grid_splits[i - 1] * search_radius
         idx_origin = sum(grid_dims[:i]) + overlap_width
         size       = domain_size - grid_splits[i - 1] * search_radius
      else:
         dims       = grid_splits[i] - grid_splits[i - 1]
         origin     = domain_origin + grid_splits[i - 1] * search_radius
         idx_origin = sum(grid_dims[:i]) + overlap_width
         size       = (grid_splits[i] - grid_splits[i - 1]) * search_radius

      grid_dims.append(dims)
      origins.append(origin)
      grid_index_origins.append(idx_origin)
      sizes.append(size)

   return origins, sizes, grid_dims, grid_index_origins

def decompose_domain(setup: dict, fluid_r: np.ndarray) -> List[SubdomainGrid]:
   sr = setup["search_radius"]
   ow = setup["overlap_width"]
   sx = setup["subdomains_x"]
   sy = setup["subdomains_y"]
   sz = setup.get("subdomains_z", 1)
   n  = setup["number_of_subdomains"]
   fn = setup["fluid_n"]

   splits_x = compute_grid_splits(fluid_r, fn, n, sx, 0, sr)
   splits_y = compute_grid_splits(fluid_r, fn, n, sy, 1, sr)

   origins_x, sizes_x, dims_x, idx_origins_x = decompose_axis(
      sx, splits_x, setup["domain_origin_x"], setup["domain_size_x"], sr, ow)
   origins_y, sizes_y, dims_y, idx_origins_y = decompose_axis(
      sy, splits_y, setup["domain_origin_y"], setup["domain_size_y"], sr, ow)

   is_3d = "domain_size_z" in setup
   if is_3d:
      dims_z = int(math.ceil(setup["domain_size_z"] / sr))
      phys_z_min = setup["domain_origin_z"]
      phys_z_max = setup["domain_origin_z"] + setup["domain_size_z"]
   else:
      dims_z = None
      phys_z_min = None
      phys_z_max = None

   grids = []
   for ix in range(sx):
      for iy in range(sy):
         for iz in range(sz):
            if sx == 1:
               eps = 1.005
               lower_x = eps * setup["domain_origin_x"]
               upper_x = eps * setup["domain_size_x"]
            else:
               lower_x = origins_x[ix]
               upper_x = setup["domain_size_x"] if ix == sx - 1 else origins_x[ix + 1]

            if sy == 1:
               lower_y = setup["domain_origin_y"]
               upper_y = setup["domain_size_y"]
            else:
               lower_y = origins_y[iy]
               if iy == sy - 1:
                  lower_y = origins_y[iy] - sr
                  upper_y = setup["domain_origin_y"] + (setup["domain_size_y"] + 1) * sr
               else:
                  upper_y = origins_y[iy + 1]

            grids.append(SubdomainGrid(
               factor        = 1.0,
               search_radius = sr,
               dp            = setup["dp"],
               phys_x_min    = lower_x,
               phys_x_max    = upper_x,
               phys_y_min    = lower_y,
               phys_y_max    = upper_y,
               phys_z_min    = phys_z_min,
               phys_z_max    = phys_z_max,
               origin_glob_x = idx_origins_x[ix],
               origin_glob_y = idx_origins_y[iy],
               origin_glob_z = ow if is_3d else None,
               dims_x        = dims_x[ix],
               dims_y        = dims_y[iy],
               dims_z        = dims_z,
               ix            = ix,
               iy            = iy,
               iz            = iz,
            ))

   return grids

def split_particles_to_subdomains(
   grids:      List[SubdomainGrid],
   fluid_r:    np.ndarray,
   boundary_r: np.ndarray,
) -> dict:
   fluid_rx = fluid_r[:, 0]
   fluid_ry = fluid_r[:, 1]
   bound_rx = boundary_r[:, 0]
   bound_ry = boundary_r[:, 1]

   splits = {}
   for g in grids:
      fluid_mask = (
         (fluid_rx > g.phys_x_min) & (fluid_rx <= g.phys_x_max) &
         (fluid_ry > g.phys_y_min) & (fluid_ry <= g.phys_y_max)
      )
      boundary_mask = (
         (bound_rx > g.phys_x_min) & (bound_rx <= g.phys_x_max) &
         (bound_ry > g.phys_y_min) & (bound_ry <= g.phys_y_max)
      )

      sub_fluid    = fluid_r[fluid_mask]
      sub_boundary = boundary_r[boundary_mask]

      g.fluid_n    = len(sub_fluid)
      g.boundary_n = len(sub_boundary)

      splits[(g.ix, g.iy, g.iz)] = (sub_fluid, sub_boundary)

   return splits


# --- Config data building (shared by multiresolution and distributed) ---

def build_subdomain_data(
   g:             SubdomainGrid,
   setup:        dict,
   fluid_path:    str,
   boundary_path: str,
) -> dict:
   overlap = setup.get("overlap_width", 0)
   is_3d = "domain_size_z" in setup
   axes  = ["x", "y", "z"] if is_3d else ["x", "y"]
   domain_origin = {ax: setup[f"domain_origin_{ax}"] for ax in axes}

   fluid_alloc    = setup.get("fluid_alloc_factor", setup.get("fine_alloc_fact", 2))
   boundary_alloc = setup.get("boundary_alloc_factor", setup.get("fine_alloc_fact", 2))

   fn_alloc = fluid_alloc * g.fluid_n if g.fluid_n > 0 else fluid_alloc * setup.get("fluid_n", 0)
   bn_alloc = boundary_alloc * g.boundary_n if g.boundary_n > 0 else setup.get("boundary_n", 0)

   sd: dict = {
      "fluid-particles":      fluid_path,
      "boundary-particles":   boundary_path,
      "fluid_n":              g.fluid_n,
      "boundary_n":           g.boundary_n,
      "fluid_n_allocated":    fn_alloc,
      "boundary_n_allocated": bn_alloc,
   }

   #if g.factor != 1.0:
   #   sd["refinement-factor"] = g.factor

   for ax in axes:
      origin_glob = getattr(g, f"origin_glob_{ax}")
      dims        = getattr(g, f"dims_{ax}")
      origin_phys = domain_origin[ax] + g.search_radius * (origin_glob - overlap)
      size_phys   = g.search_radius * dims

      sd[f"origin-global-coords-{ax}"] = origin_glob
      sd[f"grid-dimensions-{ax}"]      = dims
      sd[f"origin-{ax}"]               = round(origin_phys, 7)
      sd[f"size-{ax}"]                 = round(size_phys, 7)

   return sd

def build_distributed_domain_data(
   grids:       List[SubdomainGrid],
   setup:       dict,
   distributed: bool = False,
   particles_filename_pattern = ""
) -> dict:
   if particles_filename_pattern != "":
       particles_filename_pattern += "_"

   if distributed:
      data = {}
      for g in grids:
         key = f"subdomain-x-{g.ix}-y-{g.iy}"
         fluid_path = f"sources/{particles_filename_pattern}subdomain-x-{g.ix}-y-{g.iy}-fluid.vtk"
         boundary_path = f"sources/{particles_filename_pattern}ubdomain-x-{g.ix}-y-{g.iy}-boundary.vtk"
         data[key] = build_subdomain_data(g, setup, fluid_path, boundary_path)
      return data
   else:
      subdomains = {}
      for i, g in enumerate(grids):
         fluid_path = f"sources/subdomain-{i}-{particles_filename_pattern}fluid.vtk"
         boundary_path = f"sources/subdomain-{i}-{particles_filename_pattern}boundary.vtk"
         subdomains[f"subdomain-{i}"] = build_subdomain_data(g, setup, fluid_path, boundary_path)
      return {"subdomains": subdomains}

def write_distributed_domain_params(
   grids:       List[SubdomainGrid],
   setup:       dict,
   distributed: bool = False,
   particles_filename_pattern = ""
) -> None:
   data = build_distributed_domain_data(grids, setup, distributed, particles_filename_pattern)
   with open("sources/config-distributed-domain.jsonc", "w") as f:
      json.dump(data, f, indent=4, sort_keys=True)
