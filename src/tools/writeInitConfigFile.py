def safe_replace(text: str, replacements: dict, setup: dict) -> str:
    for placeholder, key, fmt in replacements:
        value = setup.get(key)
        if value is not None:
            text = text.replace(placeholder, format(value, fmt))
    return text

def write_simulation_params(setup: dict) -> None:
    ini_replacements = [
        # (placeholder,                          setup key,           format spec)
        ("placeholderSearchRadius",              "search_radius",     ""),
        ("placeholderDomainOrigin-x",            "domain_origin_x",   ".7f"),
        ("placeholderDomainOrigin-y",            "domain_origin_y",   ".7f"),
        ("placeholderDomainSize-x",              "domain_size_x",     ".7f"),
        ("placeholderDomainSize-y",              "domain_size_y",     ".7f"),
        ("placeholderInitParticleDistance",      "dp",                ""),
        ("placeholderSmoothingLength",           "smoothing_length",  ""),
        ("placeholderMass",                      "particle_mass",     ""),
        ("placeholderSpeedOfSound",              "speed_of_sound",    ""),
        ("placeholderDensity",                   "density",           ""),
        ("placeholderTimeStep",                  "time_step",         ""),
        ("placeholderCFL",                       "cfl",               ""),
        ("placeholderAlpha",                     "alpha",             ""),
        ("placeholderDynamicVicosity",           "dynamic_viscosity", ""),
        ("placeholderFluidParticles",            "fluid_n",           ""),
        ("placeholderAllocatedFluidParticles",   "fluid_n",           ""),
        ("placeholderBoundaryParticles",         "boundary_n",        ""),
        ("placeholderAllocatedBoundaryParticles","boundary_n",        ""),
    ]

    header_replacements = [
        # (placeholder,                          setup key,           format spec)
        ("#placeholderBoundaryConditionsType",   "bc_type",           ""),
        ("#placeholderDiffusiveTerm",            "diffusive_term",    ""),
        ("#placeholderViscosTerm",               "viscous_term",      ""),
    ]

    with open("template/config_template.ini", "r") as f:
        cfg = safe_replace(f.read(), ini_replacements, setup)
    with open("sources/config.ini", "w") as f:
        f.write(cfg)

    with open("template/config_template.h", "r") as f:
        hdr = safe_replace(f.read(), header_replacements, setup)
    with open("template/config.h", "w") as f:
        f.write(hdr)

def write_measuretool_params(setup: dict) -> None:
    mt_repmacements = [
        # (placeholder,                          setup key,           format spec)
    ]

    with open("template/config-measuretool_template.ini", "r") as f:
        cfg = safe_replace(f.read(), mt_repmacements, setup)
    with open("sources/config-measuretool.ini", "w") as f:
        f.write(cfg)
