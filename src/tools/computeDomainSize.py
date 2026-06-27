"""
Domain size computation for TNL-SPH simulation setups.

Expands the domain bounding box (``domain_origin_*`` / ``domain_end_*``) by
one search-radius cell layer and derives ``domain_size_*``.  The expansion
factor ``eps`` defaults to ``1.005``.  Axes without ``domain_origin_*`` /
``domain_end_*`` keys are skipped, so the same tool works for 2D and 3D cases.

Setup contract
--------------
Reads:
    setup["search_radius"]
    setup["domain_origin_x"], setup["domain_end_x"]   (and y, z if present)

Writes (in-place):
    setup["domain_origin_x"]  (overwritten with expanded value)
    setup["domain_size_x"]    (and y, z if present)
"""

def compute_domain_size(setup: dict, eps: float = 1.005) -> None:
   search_radius = setup["search_radius"]

   extra_parameters = {}
   for axis in ("x", "y", "z"):
      origin_key = f"domain_origin_{axis}"
      end_key = f"domain_end_{axis}"
      if origin_key not in setup or end_key not in setup:
         continue
      domain_origin = eps * (setup[origin_key] - search_radius)
      domain_end = eps * (setup[end_key] + search_radius)
      extra_parameters[origin_key] = domain_origin
      extra_parameters[f"domain_size_{axis}"] = domain_end - domain_origin

   setup.update(extra_parameters)
