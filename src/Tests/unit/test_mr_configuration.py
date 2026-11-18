"""Multiresolution domain-corner and refinement-chain configuration checks."""

from __future__ import annotations


def prepare():
    """Import the shared preprocessing tools."""
    import suite.bootstrap  # noqa: F401
    import decomposition

    return {"dec": decomposition}


_SETUP = {
    "search_radius": 0.125,
    "dp": 0.125,
    "domain_origin_x": -0.25,
    "domain_origin_y": 0.125,
    # Solver schema name retained for compatibility; values are far corners.
    "domain_size_x": 0.875,
    "domain_size_y": 0.625,
    "fluid_n": 3,
    "boundary_n": 2,
}


def check_max_corner_with_nonzero_origin(ctx):
    """Given non-zero origins, grid dimensions use far corner minus origin."""
    grid = ctx["dec"].build_subdomain_grids(
        [ctx["dec"].SubdomainDef(factor=1.0)], _SETUP
    )[0]
    expected = {
        "dims_x": 9,
        "dims_y": 4,
        "phys_x_min": -0.25,
        "phys_x_max": 0.875,
        "phys_y_min": 0.125,
        "phys_y_max": 0.625,
    }
    for name, want in expected.items():
        got = getattr(grid, name)
        if abs(got - want) > 1e-12:
            return f"{name}: expected far-corner result {want}, got {got}"
    return None


def check_ordered_one_to_two_chain(ctx):
    """Given factors 1, 1/2, 1/4, each adjacent interface remains valid."""
    dec = ctx["dec"]
    defs = [dec.SubdomainDef(factor=factor) for factor in (1.0, 0.5, 0.25)]
    try:
        grids = dec.build_subdomain_grids(defs, _SETUP)
    except ValueError as error:
        return f"ordered 1:2 chain rejected: {error}"
    if len(grids) != 3:
        return f"ordered chain produced {len(grids)} grids instead of 3"
    return None


def check_non_one_to_two_chain_rejected(ctx):
    """Given a 1:3 adjacent interface, preprocessing rejects it explicitly."""
    dec = ctx["dec"]
    defs = [dec.SubdomainDef(factor=1.0), dec.SubdomainDef(factor=1.0 / 3.0)]
    try:
        dec.build_subdomain_grids(defs, _SETUP)
    except ValueError as error:
        if "1:2" not in str(error):
            return f"ratio error does not name the supported 1:2 interface: {error}"
        return None
    return "non-1:2 adjacent refinement was accepted without a runtime guard"


def check_parser_facing_subdomain_keys(ctx):
    """Nested subdomain data flattens to the solver's camelCase key schema."""
    dec = ctx["dec"]
    grid = dec.build_subdomain_grids([dec.SubdomainDef(factor=1.0)], _SETUP)[0]
    grid.fluid_n = _SETUP["fluid_n"]
    grid.boundary_n = _SETUP["boundary_n"]
    data = dec.build_distributed_domain_data([grid], _SETUP)
    if set(data) != {"subdomains"}:
        return f"unexpected top-level serialized sections: {sorted(data)}"
    nested = data["subdomains"]["subdomain0"]

    # Equivalent of flattenJSON's parentKey + capitalized childKey convention.
    flattened = {
        "subdomain0" + key[0].upper() + key[1:]: value
        for key, value in nested.items()
    }
    expected = {
        "subdomain0RefinementFactor": 1.0,
        "subdomain0OriginGlobalCoordsX": 0,
        "subdomain0GridDimensionsX": 9,
        "subdomain0OriginX": -0.25,
        "subdomain0SizeX": 1.125,
        "subdomain0GridDimensionsY": 4,
        "subdomain0SizeY": 0.5,
    }
    for key, want in expected.items():
        got = flattened.get(key)
        if got is None or abs(got - want) > 1e-7:
            return f"{key}: expected {want}, got {got}; keys={sorted(flattened)}"
    return None


def check_parser_facing_buffer_keys(ctx):
    """MR and boundary-ghost sections retain their nested parser key prefixes."""
    dec = ctx["dec"]
    setup = {**_SETUP, "mrb_buffer_0_n": 4, "ghost_buffer_0_n": 6}
    grid = dec.build_subdomain_grids([dec.SubdomainDef(factor=1.0)], setup)[0]
    data = dec.build_distributed_domain_data([grid], setup)
    if set(data) != {"subdomains", "multiresolutionBuffers", "boundaryGhostBuffers"}:
        return f"unexpected serialized buffer sections: {sorted(data)}"

    expected = {
        "multiresolutionBuffers": (
            "multiresolutionBuffer0", "sources/multiresolution_buffer_0.vtk", 4
        ),
        "boundaryGhostBuffers": (
            "boundaryGhostBuffer0", "sources/boundary_ghost_buffer_0.vtk", 6
        ),
    }
    for section, (prefix, particles, count) in expected.items():
        nested = data[section].get(prefix)
        if nested != {"particles": particles, "n": count}:
            return f"{section}[{prefix!r}]: expected {particles}/{count}, got {nested}"
        flattened = {
            prefix + key[0].upper() + key[1:]: value for key, value in nested.items()
        }
        if flattened != {prefix + "Particles": particles, prefix + "N": count}:
            return f"wrong flattened {section} schema: {flattened}"
    return None
