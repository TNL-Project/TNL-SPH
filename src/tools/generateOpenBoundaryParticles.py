"""
General-purpose open boundary (inlet/outlet) particle generator for TNL-SPH.

This tool extrudes a surface discretization read from a VTK polydata file
along the surface normal to create the buffer layers required by the open
boundary condition.  The surface is expected to carry ``Points``,
``Normals`` and ``Areas`` point-data arrays (as produced by the typical
pre-processing pipeline used in the TNL-SPH examples).

The function is a generalisation of the ``generate_open_boundary_particles``
helper that used to live inline in several case ``init.py`` scripts.  It
supports user-supplied velocity and density profiles so that non-trivial
inlet conditions (parabolic velocity, hydrostatic density, ...) can be
injected without modifying the tool.

Profiles
--------
``velocity_profile`` and ``density_profile`` accept one of:

* ``None``        - use the default (zero velocity / constant ``setup["density"]``)
* scalar / tuple  - constant value (velocity as ``(vx, vy, vz)``, density as ``rho0``)
* callable        - ``fn(r, setup, prefix) -> ndarray`` where ``r`` has shape ``(n, 3)``

When ``coord_mode == "relative"`` the coordinates passed to the profile
callables are shifted so that the minimum corner of the generated buffer is
at the origin.  When ``coord_mode == "absolute"`` (default) the raw particle
coordinates are passed.

Setup contract
--------------
Reads:
    setup["dp"], setup["search_radius"], setup["{prefix}_filename"]
    setup["density"]  (used only when density_profile is None)

Writes (in-place):
    setup["{prefix}_n"]
    setup["{prefix}_position_x"], setup["{prefix}_position_y"], setup["{prefix}_position_z"]
    setup["{prefix}_size_x"],     setup["{prefix}_size_y"],     setup["{prefix}_size_z"]
    setup["{prefix}_orientation_x"], setup["{prefix}_orientation_y"], setup["{prefix}_orientation_z"]
    setup["{prefix}_layers"]
    setup["{prefix}_width"]
"""

import numpy as np
import vtk
from vtk.numpy_interface import dataset_adapter as dsa

import saveParticlesVTK

def _resolve_velocity(velocity_profile, r, r_profile, setup, prefix, np_patch):
   """Return a (np_patch, 3) velocity array from the profile argument."""
   if velocity_profile is None:
      return np.zeros((np_patch, 3))
   if callable(velocity_profile):
      v = np.asarray(velocity_profile(r_profile, setup, prefix), dtype=float)
      if v.ndim == 1:
         if v.shape[0] == 3:
            v = np.broadcast_to(v, (np_patch, 3)).copy()
         elif v.shape[0] == np_patch:
            v = np.column_stack([v, np.zeros(np_patch), np.zeros(np_patch)])
         else:
            raise ValueError(
               f"velocity_profile returned array of shape {v.shape}; "
               f"expected (n,3), (3,) or (n,)"
            )
      return v
   # constant tuple / list
   v_const = np.asarray(velocity_profile, dtype=float).reshape(-1)
   if v_const.shape[0] != 3:
      raise ValueError(
         f"velocity_profile constant must have 3 components, got {v_const.shape[0]}"
      )
   return np.broadcast_to(v_const, (np_patch, 3)).copy()

def _resolve_density(density_profile, r_profile, setup, prefix, np_patch):
   """Return a (np_patch,) density array from the profile argument."""
   if density_profile is None:
      return setup["density"] * np.ones(np_patch)
   if callable(density_profile):
      rho = np.asarray(density_profile(r_profile, setup, prefix), dtype=float).reshape(-1)
      return rho
   # scalar constant
   return float(density_profile) * np.ones(np_patch)

def generate_open_boundary_particles(
   setup,
   prefix,
   velocity_profile=None,
   density_profile=None,
   coord_mode="absolute",
   output_filename=None,
   verbose=False,
):
   """
   Generate open boundary buffer particles by extruding a surface along its normal.

   Parameters
   ----------
   setup : dict
      Simulation setup dictionary (see module docstring for read/written keys).
   prefix : str
      Patch identifier, e.g. ``"inlet"`` or ``"outlet"``.  The VTK file is read
      from ``setup[prefix + "_filename"]``.
   velocity_profile : None | callable | tuple
      Velocity profile.  ``None`` -> zero velocity.  A 3-tuple/list -> constant
      velocity.  A callable ``fn(r, setup, prefix) -> (n,3)`` -> per-particle.
   density_profile : None | callable | float
      Density profile.  ``None`` -> ``setup["density"]``.  A scalar -> constant.
      A callable ``fn(r, setup, prefix) -> (n,)`` -> per-particle.
   coord_mode : str
      ``"absolute"`` (default) passes raw coordinates to profile callables.
      ``"relative"`` passes coordinates shifted so the buffer min-corner is origin.
   output_filename : str, optional
      Destination VTK path.  Defaults to ``f"sources/{prefix}.vtk"``.
   verbose : bool
      Print layer count / array shapes / resulting patch parameters.

   Returns
   -------
   int
      Number of generated particles (also stored as ``setup[prefix + "_n"]``).
   """
   dp = setup["dp"]
   search_radius = setup["search_radius"]
   inlet_layers = int(np.ceil(search_radius / dp))

   # --- read the surface discretization -------------------------------------
   reader = vtk.vtkPolyDataReader()
   reader.SetFileName(setup[prefix + "_filename"])
   reader.ReadAllScalarsOn()
   reader.ReadAllVectorsOn()
   reader.Update()
   polydata = reader.GetOutput()

   r_loaded = dsa.WrapDataObject(polydata).Points
   n_loaded = dsa.WrapDataObject(polydata).PointData["Normals"]
   elem_sizes_loaded = np.array(dsa.WrapDataObject(polydata).PointData["Areas"])
   normal = n_loaded[0]

   # --- stack positions for the buffer layers -------------------------------
   r = np.empty((0, r_loaded.shape[1]), dtype=r_loaded.dtype)
   n = np.empty((0, n_loaded.shape[1]), dtype=n_loaded.dtype)
   elem_size = np.empty((0,), dtype=elem_sizes_loaded.dtype)

   if verbose:
      print(f"{prefix}_normals: {normal}")
      print(f"{prefix}_layers: {inlet_layers}")
      print(
         f"r_loaded.shape: {r_loaded.shape}, "
         f"elem_sizes_loaded.shape: {elem_sizes_loaded.shape}"
      )
      print(
         f"r.shape: {r.shape}, n.shape: {n.shape} "
         f"elementSize: {elem_size.shape}"
      )

   for layer in range(inlet_layers):
      r_layer = r_loaded - layer * dp * normal
      r = np.vstack([r, r_layer])
      n = np.vstack([n, n_loaded])
      elem_size = np.concatenate([elem_size, elem_sizes_loaded])

   # add interface offset (shift half a particle into the buffer)
   r -= 0.5 * dp * normal

   np_patch = len(r)

   # --- resolve velocity / density profiles --------------------------------
   if coord_mode == "relative":
      r_profile = r - np.min(r, axis=0)
   elif coord_mode == "absolute":
      r_profile = r
   else:
      raise ValueError(
         f"coord_mode must be 'absolute' or 'relative', got '{coord_mode}'"
      )

   v = _resolve_velocity(velocity_profile, r, r_profile, setup, prefix, np_patch)
   rho = _resolve_density(density_profile, r_profile, setup, prefix, np_patch)
   p = np.zeros(np_patch)
   ptype = np.zeros(np_patch)

   if verbose:
      print(
         f"rshape: {r.shape}, vshape: {v.shape}, rhoshape: {rho.shape}, "
         f"nshape: {n.shape}, elem_size {elem_size.shape}"
      )

   # --- write the VTK -------------------------------------------------------
   if output_filename is None:
      output_filename = f"sources/{prefix}.vtk"
   polydata_out = saveParticlesVTK.create_pointcloud_polydata(
      r, v, rho, p, ptype, normals=n, elementSize=elem_size
   )
   saveParticlesVTK.save_polydata(polydata_out, output_filename)

   # --- compute patch descriptor (preserves the legacy "disgusting" logic) --
   x_min = min(r[:, 0])
   y_min = min(r[:, 1])
   z_min = min(r[:, 2])
   x_max = max(r[:, 0])
   y_max = max(r[:, 1])
   z_max = max(r[:, 2])

   # disgusting correction: along a positive-normal axis the buffer grows
   # towards smaller coordinates, so min/max must be swapped before the
   # offset is counted back.
   if normal[0] > 0:
      x_min, x_max = x_max, x_min
   if normal[1] > 0:
      y_min, y_max = y_max, y_min
   if normal[2] > 0:
      z_min, z_max = z_max, z_min

   # offset correction: count back the half-dp shift applied above
   x_min += 0.5 * dp * normal[0]
   y_min += 0.5 * dp * normal[1]
   z_min += 0.5 * dp * normal[2]

   # disgusting size correction: the patch has zero extent along its normal
   size_x = np.abs(x_max - x_min)
   size_y = np.abs(y_max - y_min)
   size_z = np.abs(z_max - z_min)
   if np.abs(normal[0]) > 1e-6:
      size_x = 0
   if np.abs(normal[1]) > 1e-6:
      size_y = 0
   if np.abs(normal[2]) > 1e-6:
      size_z = 0

   openbc_parameters = {
      f"{prefix}_n": np_patch,
      f"{prefix}_position_x": x_min,
      f"{prefix}_position_y": y_min,
      f"{prefix}_position_z": z_min,
      f"{prefix}_size_x": size_x,
      f"{prefix}_size_y": size_y,
      f"{prefix}_size_z": size_z,
      f"{prefix}_orientation_x": n[0, 0],
      f"{prefix}_orientation_y": n[0, 1],
      f"{prefix}_orientation_z": n[0, 2],
      f"{prefix}_layers": inlet_layers,
      f"{prefix}_width": dp * inlet_layers,
   }
   setup.update(openbc_parameters)

   if verbose:
      from pprint import pprint

      pprint(openbc_parameters)

   return np_patch
