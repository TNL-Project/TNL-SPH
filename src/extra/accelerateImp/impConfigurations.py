def emptyTest( case_dir ):
    return 0, 0

configurations = []

for timeScheme in ["MidpointScheme", "MidpointSchemeWithEnergySecant"]:
    configurations.append({
        "case-tag" : f"damBreak2D_WCSPH-BI_hr-disipative_dp0.002_hdp2_{timeScheme}",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "diffusive-term" : "MolteniDiffusiveTerm",
        "time-integration" : f"{timeScheme}",
        "dp" : 0.002,
        "h-coef" : 2,
        "evaluation-function" : emptyTest
    })

for timeScheme in ["MidpointScheme", "MidpointSchemeWithEnergySecant"]:
    configurations.append({
        "case-tag" : f"damBreak2D_WCSPH-BI_hr-conservative_hdp4_{timeScheme}",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "viscous-term" : "None",
        "diffusive-term" : "None",
        "time-integration" : f"{timeScheme}",
        "dp" : 0.00075,
        "h-coef" : 4,
        "evaluation-function" : emptyTest
    })

configurations.append({
    "case-tag" : "damBreak2D_WCSPH-BI_default",
    "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
    "evaluation-function" : emptyTest
})

