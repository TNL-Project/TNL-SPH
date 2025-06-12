import sys
sys.path.append('evaluateExamplesMetrics')
import evaluateExamplesMetrics

wcsph_rsph_configurations = [
        # dam break 2D
        {
            "case-tag" : "RSPH/damBreak2D_RSPH",
            "case" : "RSPH/damBreak2D_RSPH",
            "evaluation-function" : None
        }
]

wcsph_dbc_configurations = [
        # dam break 2D
        {
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
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
        # dam break 2D
        {
            "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_SVs-AV",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "time-integration" : "SymplecticVerletScheme",
            "h-coef" :  2,
            "cfl" : 0.05,
            "bc-correction" : "ElasticBounce",
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_BI
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
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_BI
        },
        {
            "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_hr-conservative",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "None",
            "diffusive-term" : "None",
            "dp" : 0.00075,
            "h-coef" : 4,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_BI
        },
        {
            "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_hr-conservative-MGVT-DTNone",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "PhysicalViscosity_MGVT",
            "diffusive-term" : "None",
            "h-coef" : 2,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_BI
        },
        {
            "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_hr-conservative-MGVT",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "PhysicalViscosity_MGVT",
            "diffusive-term" : "MolteniDiffusiveTerm",
            "h-coef" : 2,
            "evaluation-function" : evaluateExamplesMetrics.damBreak2D_WCSPH_BI
        },
        # dam break 3D
        {
            "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_hr-conservative",
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "bc-correction" : "ElasticBounceLight",
            "time-integration" : "MidpointScheme",
            "viscous-term" : "None",
            "diffusive-term" : "None",
            "dp" : 0.01,
            "h-coef" :  3,
            "evaluation-function" : evaluateExamplesMetrics.damBreak3D_WCSPH_BI
            # NOTE: Doesn't work! Requires midpoint?
        },
        {
            "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_hr-MGVT",
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConsistent_numeric",
            "viscous-term" : "PhysicalViscosity_MGVT",
            "dp" : 0.01,
            "h-coef" :  2,
            "cfl" : 0.05,
            "evaluation-function" : evaluateExamplesMetrics.damBreak3D_WCSPH_BI
        },
        {
            "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_conservative-MGVT",
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "PhysicalViscosity_MGVT",
            "h-coef" :  2,
            "cfl" : 0.05,
            "evaluation-function" : evaluateExamplesMetrics.damBreak3D_WCSPH_BI
            # NOTE: With dp=0.01, this works even with default CFL, but for default dp=0.02,
            #       it requires lower CFL
        },
        {
            "case-tag" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI",
            "case" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI",
            "evaluation-function" : evaluateExamplesMetrics.poiseuilleFlowWithOpenBoundary2D_WCSPH_BI
        },
        {
            "case-tag" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI",
            "case" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI",
            "evaluation-function" : evaluateExamplesMetrics.poiseuilleFlowWithPeriodicBoundary2D_WCSPH
        }
]

test_default_examples = [
        # RSPH:
        {
            "case-tag" : "RSPH/damBreak2D_RSPH_default",
            "case" : "RSPH/damBreak2D_RSPH",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        # WCSPH-DBC:
        {
            "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC_default",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        {
            "case-tag" : "WCSPH-DBC/damBreak3D_WCSPH-DBC_default",
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        {
            "case-tag" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC_default",
            "case" : "WCSPH-DBC/poiseuilleFlowWithOpenBoundary2D_WCSPH-DBC",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        {
            "case-tag" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC_default",
            "case" : "WCSPH-DBC/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-DBC",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        # WCSPH-BI:
        {
            "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_default",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        {
            "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_default",
            "case" : "WCSPH-BI/damBreak3D_WCSPH-BI",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        {
            "case-tag" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI_default",
            "case" : "WCSPH-BI/poiseuilleFlowWithOpenBoundary2D_WCSPH-BI",
            "evaluation-function" : evaluateExamplesMetrics.empty
        },
        {
            "case-tag" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI_default",
            "case" : "WCSPH-BI/poiseuilleFlowWithPeriodicBoundary2D_WCSPH-BI",
            "evaluation-function" : evaluateExamplesMetrics.empty
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
]
