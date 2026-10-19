"""wcsph_bi regression configurations (port of configurationsToTest.py:85-181)."""

configurations = [
    # dam break 2D
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_SVs-AV",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConsistent_numeric",
        "time-integration" : "SymplecticVerletScheme",
        "h-coef" :  2,
        "cfl" : 0.05,
        "bc-correction" : "ElasticBounce",
        "evaluation-function" : "damBreak2D_WCSPH_BI"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_SVs-MGVT",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConsistent_numeric",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "time-integration" : "SymplecticVerletScheme",
        "h-coef" :  2,
        "cfl" : 0.05,
        "bc-correction" : "ElasticBounce",
        "evaluation-function" : "damBreak2D_WCSPH_BI"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_hr-conservative",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "viscous-term" : "None",
        "diffusive-term" : "None",
        "dp" : 0.00075,
        "h-coef" : 4,
        "evaluation-function" : "damBreak2D_WCSPH_BI"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_hr-conservative-MGVT-DTNone",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "diffusive-term" : "None",
        "h-coef" : 2,
        "evaluation-function" : "damBreak2D_WCSPH_BI"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_hr-conservative-MGVT",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "diffusive-term" : "MolteniDiffusiveTerm",
        "h-coef" : 2,
        "evaluation-function" : "damBreak2D_WCSPH_BI"
    },
    ## dam break 3D
    #{
    #    "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_hr-conservative",
    #    "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
    #    "bc-type" : "BIConservative_numeric",
    #    "bc-correction" : "ElasticBounceLight",
    #    "time-integration" : "MidpointScheme",
    #    "viscous-term" : "None",
    #    "diffusive-term" : "None",
    #    "dp" : 0.01,
    #    "h-coef" :  3,
    #    "evaluation-function" : "damBreak3D_WCSPH_BI"
    #    # NOTE: Doesn't work! Requires midpoint?
    #},
    {
        "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_hr-MGVT",
        "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
        "bc-type" : "BIConsistent_numeric",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "dp" : 0.01,
        "h-coef" :  2,
        "cfl" : 0.05,
        "evaluation-function" : "damBreak3D_WCSPH_BI"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_conservative-MGVT",
        "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
        "bc-type" : "BIConservative_numeric",
        "viscous-term" : "PhysicalViscosity_MGVT",
        "h-coef" :  2,
        "cfl" : 0.05,
        "evaluation-function" : "damBreak3D_WCSPH_BI"
        # NOTE: With dp=0.01, this works even with default CFL, but for default dp=0.02,
        #       it requires lower CFL
    },
    {
        "case-tag" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI",
        "case" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI",
        "evaluation-function" : "poiseuilleFlowWithOpenBoundary2D_WCSPH_BI"
    },
    {
        "case-tag" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI",
        "case" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI",
        "evaluation-function" : "poiseuilleFlowWithPeriodicBoundary2D_WCSPH"
    }
]
