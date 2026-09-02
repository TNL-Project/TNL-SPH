import json
import os
import json
import numpy as np

import numpy as np

ini_replacements = [
   # (placeholder,                          setup key,               format spec)
   ("placeholderSearchRadius",              "search_radius",         ""),
   ("placeholderDomainOrigin-x",            "domain_origin_x",       ".7f"),
   ("placeholderDomainOrigin-y",            "domain_origin_y",       ".7f"),
   ("placeholderDomainOrigin-z",            "domain_origin_z",       ".7f"),
   ("placeholderDomainSize-x",              "domain_size_x",         ".7f"),
   ("placeholderDomainSize-y",              "domain_size_y",         ".7f"),
   ("placeholderDomainSize-z",              "domain_size_z",         ".7f"),
   ("placeholderInitParticleDistance",      "dp",                    ""),
   ("placeholderSmoothingLength",           "smoothing_length",      ""),
   ("placeholderMass",                      "particle_mass",         ""),
   ("placeholderBoundaryMass",              "particle_mass",         ""),
   ("placeholderBoundaryElementSize",       "boundary_element_size", ""),
   ("placeholderSpeedOfSound",              "speed_of_sound",        ""),
   ("placeholderDensity",                   "density",               ""),
   ("placeholderTimeStep",                  "time_step",             ""),
   ("placeholderNumberOfSubdomains",        "number_of_subdomains",  ""),
   ("placeholderCFL",                       "cfl",                   ""),
   ("placeholderAlpha",                     "alpha",                 ""),
   ("placeholderDynamicVicosity",           "dynamic_viscosity",     ""),
   ("placeholderFluidParticles",            "fluid_n",               ""),
   ("placeholderAllocatedFluidParticles",   "allocated_fluid_n",     ""),
   ("placeholderBoundaryParticles",         "boundary_n",            ""),
   ("placeholderAllocatedBoundaryParticles","allocated_boundary_n",  ""),
   ("placeholderSubdomains-x",              "subdomains_x",          ""),
   ("placeholderSubdomains-y",              "subdomains_y",          ""),
   ("placeholderSubdomains-z",              "subdomains_z",          ""),
]

header_replacements = [
   # (placeholder,                          setup key,               format spec)
   ("#placeholderBoundaryConditionsType",   "bc_type",               ""),
   ("#placeholderDiffusiveTerm",            "diffusive_term",        ""),
   ("#placeholderViscosTerm",               "viscous_term",          ""),
   ("#placeholderTimeIntegration",          "time_integration",      ""),
]

mt_repmacements = [
   # (placeholder,                          setup key,           format spec)
   ("placeholderInitParticleDistance",      "dp",                    ""),
   ("placeholderSmoothingLength",           "smoothing_length",      ""),
]

open_boundary_replacements = [
   # (placeholder,                          setup key,                    format spec)
   ("placeholderInletParticles",            "inlet_n",                    ""),
   ("placeholderAllocatedInletParticles",   "allocated_inlet_n",          ""),
   ("placeholderOutletParticles",           "outlet_n",                   ""),
   ("placeholderAllocatedOutletParticles",  "allocated_outlet_n",         ""),
   ("placeholderInletVelocity_x",           "inlet_velocity_x",           ""),
   ("placeholderInletVelocity_y",           "inlet_velocity_y",           ""),
   ("placeholderInletVelocity_z",           "inlet_velocity_z",           ""),
   ("placeholderInletPosition1_x",          "inlet_position1_x",          ""),
   ("placeholderInletPosition1_y",          "inlet_position1_y",          ""),
   ("placeholderInletPosition1_z",          "inlet_position1_z",          ""),
   ("placeholderInletPosition2_x",          "inlet_position2_x",          ""),
   ("placeholderInletPosition2_y",          "inlet_position2_y",          ""),
   ("placeholderInletPosition2_z",          "inlet_position2_z",          ""),
   ("placeholderInletDensity",              "density",                    ""),
   ("placeholderInletWidth_x",              "inlet_width",                ".7f"),
   ("placeholderOutletVelocity_x",          "outlet_velocity_x",          ""),
   ("placeholderOutletVelocity_y",          "outlet_velocity_y",          ""),
   ("placeholderOutletVelocity_z",          "outlet_velocity_z",          ""),
   ("placeholderOutletPosition1_x",         "outlet_position1_x",         ""),
   ("placeholderOutletPosition1_y",         "outlet_position1_y",         ""),
   ("placeholderOutletPosition1_z",         "outlet_position1_z",         ""),
   ("placeholderOutletPosition2_x",         "outlet_position2_x",         ""),
   ("placeholderOutletPosition2_y",         "outlet_position2_y",         ""),
   ("placeholderOutletPosition2_z",         "outlet_position2_z",         ""),
   ("placeholderOutletDensity",             "density",                    ""),
   ("placeholderOutletWidth_x",             "outlet_width",               ".7f"),
]

def safe_replace(text: str, replacements: list, setup: dict) -> str:
   for placeholder, key, fmt in replacements:
      value = setup.get(key)
      if value is not None:
         text = text.replace(placeholder, format(value, fmt))
   return text

def write_simulation_params(setup: dict) -> None:
   # Build a config-writing copy with standardized keys and defaults so the
   # caller's setup dict is not polluted (keeps JSON dumps stable).
   cfg_setup = dict(setup)
   # Support fluid_points_count as an alias for fluid_n
   if "fluid_n" not in cfg_setup and "fluid_points_count" in cfg_setup:
      cfg_setup["fluid_n"] = cfg_setup["fluid_points_count"]
   # Default allocation counts to the particle counts (1x) when not given.
   cfg_setup.setdefault("allocated_fluid_n", cfg_setup.get("fluid_n", 0))
   cfg_setup.setdefault("allocated_boundary_n", cfg_setup.get("boundary_n", 0))
   cfg_setup.setdefault("boundary_element_size", 0)

   with open("template/config_template.jsonc", "r") as f:
      cfg = safe_replace(f.read(), ini_replacements, cfg_setup)

   # Handle the open-boundary-patches count placeholder.  When patches exist,
   # emit a valid JSONC reference to the external open-boundary config file.
   if "placeholderNumberOfOpenBoundaryPatches" in cfg:
      patches = cfg_setup.get("open_boundary_patches", [])
      if len(patches) > 0:
         obc_block = f'{len(patches)},\n    "open-boundary-config": "sources/config-open-boundary.jsonc"'
      else:
         obc_block = "0"
      cfg = cfg.replace("placeholderNumberOfOpenBoundaryPatches", obc_block)

   with open("sources/config.jsonc", "w") as f:
      f.write(cfg)

   if os.path.exists("template/config_template.h"):
      with open("template/config_template.h", "r") as f:
         hdr = safe_replace(f.read(), header_replacements, setup)
      with open("template/config.h", "w") as f:
         f.write(hdr)

def write_open_boundary_params(setup: dict) -> None:
   with open("template/config-open-boundary_template.jsonc", "r") as f:
       cfg = safe_replace(f.read(), open_boundary_replacements, setup)

   # Patch-interpolated placeholders (e.g. placeholder_inletPosition1_x).
   # Supports templates that use placeholder_{patch}* naming instead of the
   # fixed CamelCase placeholderInlet*/placeholderOutlet* scheme above.  Both
   # schemes can coexist in the same template (missing placeholders are no-ops).
   for patch in setup.get("open_boundary_patches", []):
      n = setup.get(f"{patch}_n")
      if n is None:
         continue
      allocated = setup.get(f"allocated_{patch}_n", 2 * n)

      cfg = cfg.replace(f"placeholder_{patch}Particles", f"{n}")
      cfg = cfg.replace(f"placeholder_{patch}AllocatedParticles", f"{allocated}")

      for axis in "xyz":
         cfg = cfg.replace(
            f"placeholder_{patch}Orientation_{axis}",
            f"{setup.get(f'{patch}_orientation_{axis}', 0)}",
         )
         cfg = cfg.replace(
            f"placeholder_{patch}Velocity_{axis}",
            f"{setup.get(f'{patch}_velocity_{axis}', 0)}",
         )
         pos = setup.get(f"{patch}_position_{axis}", 0)
         size = setup.get(f"{patch}_size_{axis}", 0)
         cfg = cfg.replace(f"placeholder_{patch}Position1_{axis}", f"{pos}")
         cfg = cfg.replace(f"placeholder_{patch}Position2_{axis}", f"{pos + size}")

      cfg = cfg.replace(f"placeholder_{patch}Density", f"{setup.get('density', 0)}")

      width = setup.get(f"{patch}_width", 0)
      cfg = cfg.replace(f"placeholder_{patch}Width_x", f"{width:.5f}")
      cfg = cfg.replace(f"placeholder_{patch}Width_y", "0.0")
      cfg = cfg.replace(f"placeholder_{patch}Width_z", "0.0")

   with open("sources/config-open-boundary.jsonc", "w") as f:
      f.write(cfg)

def write_open_boundary_params_inline(
   setup: dict,
   output_path: str = "sources/config-open-boundary.jsonc",
   particles_filename_pattern: str = "sources/{patch}.vtk",
) -> None:
   """Write the open-boundary JSONC config directly from the setup dict.

   Unlike write_open_boundary_params, this does not require a template file.
   It emits one buffer block per patch listed in setup["open_boundary_patches"],
   reading patch geometry/velocity data written by generate_open_boundary_particles.
   """
   patches = setup.get("open_boundary_patches", [])
   if not patches:
      with open(output_path, "w") as f:
         f.write("{}\n")
      return

   lines = ["{"]
   for i, patch in enumerate(patches, start=1):
      n = int(setup.get(f"{patch}_n", 0))
      allocated_n = int(setup.get(f"allocated_{patch}_n", 2 * n))

      orientation = {ax: float(setup.get(f"{patch}_orientation_{ax}", 0.0)) for ax in "xyz"}
      position = {ax: float(setup.get(f"{patch}_position_{ax}", 0.0)) for ax in "xyz"}
      size = {ax: float(setup.get(f"{patch}_size_{ax}", 0.0)) for ax in "xyz"}
      velocity = {ax: float(setup.get(f"{patch}_velocity_{ax}", 0.0)) for ax in "xyz"}
      width = float(setup.get(f"{patch}_width", 0.0))
      density = float(setup.get("density", 0.0))

      normal_axis = max(((ax, abs(orientation[ax])) for ax in "xyz"), key=lambda item: item[1])[0]
      buffer_width = {ax: 0.0 for ax in "xyz"}
      buffer_width[normal_axis] = width

      default_rho_bc = "fixed" if patch == "inlet" else "extrapolated"
      default_v_bc = "fixed" if patch == "inlet" else "extrapolated"
      rho_bc = setup.get(f"{patch}_rho_bc", default_rho_bc)
      v_bc = setup.get(f"{patch}_v_bc", default_v_bc)

      lines.append(f"    // Buffer {i} - {patch}")
      lines.append(f'    "buffer-{i}-particles": "{particles_filename_pattern.format(patch=patch)}",')
      lines.append(f'    "buffer-{i}-numberOfParticles": {n},')
      lines.append(f'    "buffer-{i}-numberOfAllocatedParticles": {allocated_n},')
      lines.append(f'    "buffer-{i}-type": "{patch}",')
      lines.append(f'    "buffer-{i}-identifier": "{patch}",')
      for ax in "xyz":
         lines.append(f'    "buffer-{i}-orientation-{ax}": {orientation[ax]:.7f},')
      for ax in "xyz":
         lines.append(f'    "buffer-{i}-position-1-{ax}": {position[ax]:.7f},')
      for ax in "xyz":
         lines.append(f'    "buffer-{i}-position-2-{ax}": {position[ax] + size[ax]:.7f},')
      for ax in "xyz":
         lines.append(f'    "buffer-{i}-bufferWidth-{ax}": {buffer_width[ax]:.7f},')
      lines.append(f'    "buffer-{i}-rho_bc": "{rho_bc}",')
      lines.append(f'    "buffer-{i}-density": {density:.7f},')
      lines.append(f'    "buffer-{i}-v_bc": "{v_bc}",')
      for ax in "xyz":
         lines.append(f'    "buffer-{i}-velocity-{ax}": {velocity[ax]:.7f},')
      if patch == "outlet":
         lines.append(f'    "buffer-{i}-extrapolationDetTreshold": 1e3,')

   if lines[-1].endswith(","):
      lines[-1] = lines[-1][:-1]
   lines.append("}")

   with open(output_path, "w") as f:
      f.write("\n".join(lines) + "\n")


def write_measuretool_params(setup: dict) -> None:
    with open("template/config-measuretool_template.jsonc", "r") as f:
        cfg = safe_replace(f.read(), mt_repmacements, setup)
    with open("sources/config-measuretool.jsonc", "w") as f:
        f.write(cfg)

def write_user_defined_params(setup: dict, udp_replacements = []) -> None:
    with open("template/config-user-defined_template.jsonc", "r") as f:
        cfg = safe_replace(f.read(), udp_replacements, setup)
    with open("sources/config-user-defined.jsonc", "w") as f:
        f.write(cfg)

def save_params_to_json(data: dict, filename: str):
   #FIX: We need to somehow deal with the
   def json_converter(obj):
      if isinstance(obj, np.generic):  # np.float32, np.int64, np.bool_, ...
         return obj.item()             # Convert to native Python type
      if isinstance(obj, np.ndarray):
         return obj.tolist()           # Convert arrays to lists
      raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")

   with open(filename, "w") as f:
      json.dump(data, f, indent=4, default=json_converter)
