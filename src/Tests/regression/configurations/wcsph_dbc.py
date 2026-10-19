"""wcsph_dbc regression configurations (port of configurationsToTest.py:15-83)."""

configurations = [
    # dam break 2D
    {
        "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "bc-type" : "DBC",
        "h-coef" :  2**0.5,
        "evaluation-function" : "damBreak2D_WCSPH_DBC"
    },
    {
        "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC:MGVT",
        "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "bc-type" : "DBC",
        "h-coef" :  2**0.5,
        "viscous-term" : "PhysicalViscosity",
        "evaluation-function" : "damBreak2D_WCSPH_DBC"
    },
    {
        "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC:h-coef-2",
        "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "bc-type" : "DBC",
        "h-coef" :  2,
        "evaluation-function" : "damBreak2D_WCSPH_DBC"
    },
    {
        "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC:MDBC",
        "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "bc-type" : "MDBC",
        "h-coef" :  2**0.5,
        "evaluation-function" : "damBreak2D_WCSPH_MDBC"
    },
    #{
    # dam break 3D
    {
        "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
        "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
        "bc-type" : "DBC",
        "dp" : 0.01,
        "h-coef" :  2,
        "evaluation-function" : "damBreak3D_WCSPH_DBC"
    },
    {
        "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC:MGVT",
        "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
        "bc-type" : "DBC",
        "dp" : 0.01,
        "h-coef" :  2,
        "viscous-term" : "PhysicalViscosity",
        "evaluation-function" : "damBreak3D_WCSPH_DBC"
    },
    {
        "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC:MDBC",
        "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
        "bc-type" : "MDBC",
        "dp" : 0.01,
        "h-coef" :  2,
        "evaluation-function" : "damBreak3D_WCSPH_MDBC"
    },
    {
        "case-tag" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC",
        "case" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC",
        "evaluation-function" : "poiseuilleFlowWithOpenBoundary2D_WCSPH"
    },
    {
        "case-tag" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC",
        "case" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC",
        "evaluation-function" : "poiseuilleFlowWithPeriodicBoundary2D_WCSPH"
    }
]
