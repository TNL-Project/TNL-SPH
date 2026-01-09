def do_nothing( case_dir ):
    return 0, 0

configurations = []
# WCSPH-BI: dam break 3D 
for dp in [0.0075, 0.01, 0.02]:
    configurations.append({
        "case-tag" : f"damBreak3D_WCSPH-BI_dp{dp}-unconditionally-stable-with-dissipation",
        "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "bc-correction" : "ElasticBounceLight",
        "time-integration" : "MidpointScheme",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "diffusive-term" : "MolteniDiffusiveTerm",
        "dp" : dp,
        "h-coef" :  3,
        "evaluation-function" : do_nothing
    })
