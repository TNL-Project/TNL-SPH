import math
from dataclasses import dataclass

@dataclass
class SubdomainDef:
    factor: float # dp = factor * dp_L0  (< 1: finer, > 1: coarser)
    x_min:  float = -math.inf;  x_max: float = math.inf
    y_min:  float = -math.inf;  y_max: float = math.inf
    z_min:  float = -math.inf;  z_max: float = math.inf

@dataclass
class SubdomainGrid:
    factor:        float
    search_radius: float # = factor * search_radius_L0
    dp:            float # = factor * dp_L0

    # Physical bounding box (clipped to global domain)
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

#TODO: Is this necessary?
@dataclass
class RectangularRefinementDef:
    """
    Describes a nested rectangular fine region inside a coarse domain.
    The coarse domain covers the full domain; particles inside fine_* bounds
    are skipped during coarse particle generation.
    fine_factor: refinement factor for the fine region (< 1.0)
    """
    fine_factor: float
    fine_x_min:  float
    fine_x_max:  float
    fine_y_min:  float
    fine_y_max:  float
    fine_z_min:  float = -math.inf
    fine_z_max:  float =  math.inf


def build_subdomain_grids(defs: List[SubdomainDef], setup: dict) -> List[SubdomainGrid]:
    """
    Convert physical SubdomainDefs into SubdomainGrid objects.
    All integer coordinates are expressed in the coarse (L0) cell system.
    Works for both 2D (no z entries in setup) and 3D.
    """
    h0     = setup["search_radius"] # coarse cell size
    is_3d  = "domain_size_z" in setup
    axes   = ["x", "y", "z"] if is_3d else ["x", "y"]

    domain_origin = {ax: setup[f"domain_origin_{ax}"] for ax in axes}
    domain_size   = {ax: setup[f"domain_size_{ax}"]   for ax in axes}
    n_cells_glob  = {ax: int(math.ceil(domain_size[ax] / h0)) for ax in axes}

    def resolve_axis(sd: SubdomainDef, ax: str, f: float) -> tuple:
        """
        Returns (origin_glob, dims, phys_min, phys_max) for one axis.
        """
        fl   = int(round(1.0 / f))
        o    = domain_origin[ax]
        n_g  = n_cells_glob[ax]

        p0   = max(getattr(sd, f"{ax}_min"), o)
        og   = max(0, int(math.floor((p0 - o) / h0)))
        p0   = o + og * h0 # snap to L0 grid (TODO: we want to snap L grid to L-1)

        p1   = min(getattr(sd, f"{ax}_max"), o + n_g * h0)
        eg   = min(n_g, int(math.floor((p1 - o) / h0)))
        dims = (eg - og) * fl

        h    = f * h0
        return og * fl, dims, p0, p0 + dims * h

    grids: List[SubdomainGrid] = []
    for sd in defs:
        f  = sd.factor
        h  = f * h0
        dp = f * setup["dp"]

        resolved = {ax: resolve_axis(sd, ax, f) for ax in axes}
        # resolved[ax] = (origin_glob, dims, phys_min, phys_max)

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
        ))

    return grids

def build_rectangular_subdomain_grids(
        rdef:  RectangularRefinementDef,
        setup: dict
) -> tuple:   # returns (coarse_grid, fine_grid)
    """
    Build two SubdomainGrid objects for a rectangular nested refinement:
      - coarse_grid: covers the full domain (factor = 1.0)
      - fine_grid:   covers the fine rectangle (factor = rdef.fine_factor)

    The coarse grid's phys bounds are the full domain — the caller is
    responsible for skipping particles inside the fine region during
    particle generation (use the 'exclusion_box' returned alongside).
    """
    # Reuse the existing linear builder for both grids
    coarse_def = SubdomainDef(factor=1.0)   # full domain
    fine_def   = SubdomainDef(
        factor = rdef.fine_factor,
        x_min  = rdef.fine_x_min,  x_max = rdef.fine_x_max,
        y_min  = rdef.fine_y_min,  y_max = rdef.fine_y_max,
        z_min  = rdef.fine_z_min,  z_max = rdef.fine_z_max,
    )
    coarse_grid, fine_grid = build_subdomain_grids([coarse_def, fine_def], setup)
    coarse_grid.index = 0
    fine_grid.index   = 1
    return coarse_grid, fine_grid

def write_distributed_domain_params(grids: List[SubdomainGrid], setup: dict) -> None:
    h0    = setup["search_radius"]
    is_3d = "domain_size_z" in setup
    axes  = ["x", "y", "z"] if is_3d else ["x", "y"]
    fact = 2  # TODO

    domain_origin = {ax: setup[f"domain_origin_{ax}"] for ax in axes}

    with open("sources/config-distributed-domain.ini", "w") as f:
        f.write("# Subdomain information\n\n")
        for i, g in enumerate(grids):
            prefix = f"subdomain-{i}-"

            # Axis-independent entries
            params = {
                "fluid-particles":      f"sources/subdomain-{i}-dambreak_fluid.vtk",
                "boundary-particles":   f"sources/subdomain-{i}-dambreak_boundary.vtk",
                "fluid_n":              g.fluid_n,
                "boundary_n":           g.boundary_n,
                "fluid_n_allocated":    fact * g.fluid_n if fact * g.fluid_n > 0 else setup["fluid_n"],
                "boundary_n_allocated": fact * g.boundary_n if fact * g.boundary_n > 0 else setup["boundary_n"],
                "refinement-factor":    g.factor,
            }

            # Axis-dependent entries — grouped by keyword, not by axis
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

            # Flatten: all x/y/z for each keyword before moving to the next
            for keyword, ax_values in ax_groups.items():
                for ax, value in ax_values.items():
                    params[f"{keyword}-{ax}"] = value

            for k, v in params.items():
                f.write(f"{prefix}{k} = {v}\n")
            f.write("\n")
