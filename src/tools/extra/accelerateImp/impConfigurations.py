def emptyTest( case_dir ):
    return 0, 0

configurations = [
        {
            "case-tag" : "damBreak2D_WCSPH-BI_hr-conservative_hdp4",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "bc-type" : "BIConservative_numeric",
            "viscous-term" : "None",
            "diffusive-term" : "None",
            "dp" : 0.00075,
            "h-coef" : 4,
            "evaluation-function" : emptyTest
        },
        {
            "case-tag" : "damBreak2D_WCSPH-BI_default",
            "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
            "evaluation-function" : emptyTest
        },
]
