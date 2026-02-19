import sys
sys.path.append('evaluateExamplesMetrics')

def do_nothing( case_dir ):
    return 0, 0

wcsph_rsph_configurations = [
        # dam break 2D
        {
            "case-tag" : "results_damBreak2D_RSPH",
            "case" : "RSPH/damBreak2D_RSPH",
            "evaluation-function" : do_nothing
        }
]

wcsph_dbc_configurations = [
        # dam break 2D
        {
            "case-tag" : "results_damBreak2D_WCSPH-DBC_DBC",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : do_nothing
        },
        {
            "case-tag" : "damBreak2D_WCSPH-DBC_MDBC",
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "h-coef" :  2**0.5,
            "evaluation-function" : do_nothing
        },
]

wcsph_bi_configurations = [
        # dam break 2D
        {
            "case-tag" : "results_damBreak2D_WCSPH-BI_Verlet",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "time-integration" : "VerletScheme",
            "h-coef" :  2,
            "speed-of-sound" : 34.3,
            "cfl" : 0.1,
            "bc-correction" : "ElasticBounce",
            "evaluation-function" : do_nothing
        },
        # dam break 2D
        {
            "case-tag" : "results_damBreak2D_WCSPH-BI_SymplecticVerlet",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "time-integration" : "SymplecticVerletScheme",
            "h-coef" :  2,
            "speed-of-sound" : 34.3,
            "cfl" : 0.1,
            "bc-correction" : "ElasticBounce",
            "evaluation-function" : do_nothing
        },
        # dam break 2D
        {
            "case-tag" : "results_damBreak2D_WCSPH-BI_RK4Scheme",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "time-integration" : "RK4Scheme",
            "h-coef" :  2,
            "speed-of-sound" : 34.3,
            "cfl" : 0.1,
            "bc-correction" : "ElasticBounce",
            "evaluation-function" : do_nothing
        },
]
