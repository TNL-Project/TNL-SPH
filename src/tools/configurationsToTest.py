import sys
sys.path.append('evaluateExamplesMetrics')
import evaluateExamplesMetrics

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
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "viscous-term" : "ArtificialViscosity",
            "h-coef" :  2**0.5,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_DBC
        },
        {
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC:MGVT",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "h-coef" :  2**0.5,
            "viscous-term" : "PhysicalViscosity",
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_DBC
        },
        {
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC:h-coef-2",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "viscous-term" : "ArtificialViscosity",
            "h-coef" :  2,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_DBC
        },
        {
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC:MDBC",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_MDBC
        },
        #{
        # dam break 3D
        {
            "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : evaluateExamplesMetrics.damBreak3D_WCSPH_DBC
        },
        {
            "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC:MGVT",
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "viscous-term" : "PhysicalViscosity",
            "evaluation-function" : evaluateExamplesMetrics.damBreak3D_WCSPH_DBC
        },
        {
            "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC:MDBC",
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : evaluateExamplesMetrics.damBreak3D_WCSPH_MDBC
        },
        {
            "case-tag" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC",
            "case" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC",
            "evaluation-function" : evaluateExamplesMetrics.poiseuilleFlowWithOpenBoundary2D_WCSPH
        },
        {
            "case-tag" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC",
            "case" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC",
            "evaluation-function" : evaluateExamplesMetrics.poiseuilleFlowWithPeriodicBoundary2D_WCSPH
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
        #{
        #    "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
        #    "bc-type" : "BIConsistent_numeric",
        #    "viscous-term" : "None",
        #    "dp" : 0.01,
        #    "h-coef" :  2,
        #    "evaluation-function" : None
        #},
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "viscous-term" : "PhysicalViscosity_MGVT",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "viscous-term" : "ArtificialViscosity",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "None",
            "dp" : 0.01,
            "h-coef" :  2,
            "evaluation-function" : None
        },
        {
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "PhysicalViscosity_MGVT",
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
        # dam break 2D
        {
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_DBC
        }
        #{
        #    "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
        #    "bc-type" : "DBC",
        #    "dp" : 0.01,
        #    "h-coef" :  2,
        #    "viscous-term" : "PhysicalViscosity",
        #    "evaluation-function" : None
        #}
]


#        {
#            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
#            "bc-type" : "DBC",
#            "dp" : 0.01,
#            "h-coef" :  2,
#            "diffusive_term" : "MolteniDiffusiveTerm",
#            "delta" : 0.1,
#            "viscous-term" : "ArtificialViscosity"
#            "alpha" : 0.02,
#        },
