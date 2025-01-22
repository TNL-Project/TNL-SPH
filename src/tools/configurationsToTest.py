wcsph_rsph_configurations = [
        # dam break 2D
        {
            "case" : 'RSPH/damBreak2D_RSPH',
            "evaluation-function" : None
        }
]

wcsph_dbc_configurations = [
        # dam break 2D
        {
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "h-coef" :  2,
            "evaluation-function" : None
        },
        # dam break 3D
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "viscos-term" : "PhysicalViscosity",
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "diffusive_term" : "FourtakasDiffusiveTerm",
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC",
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC",
            "evaluation-function" : None
        }
]

wcsph_bi_configurations = [
        # TODO:
        # dam break 2D
        {
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "DBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : None
        },
        # dam break 3D
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "viscos-term" : "None",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "viscos-term" : "PhysicalViscosity_MGVT",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "viscos-term" : "ArtificialViscosity",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscos-term" : "None",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscos-term" : "PhysicalViscosity_MGVT",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscos-term" : "ArtificialViscosity",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI",
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI",
            "evaluation-function" : None
        }
]


test_configurations = [
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "viscous-term" : "PhysicalViscosity",
            "evaluation-function" : None
        }
]


#        {
#            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
#            "bc-type" : "DBC",
#            "dp" : 0.01,
#            "h-coef" :  2,
#            "diffusive_term" : "MolteniDiffusiveTerm",
#            "delta" : 0.1,
#            "viscos-term" : "ArtificialViscosity"
#            "alpha" : 0.02,
#        },
