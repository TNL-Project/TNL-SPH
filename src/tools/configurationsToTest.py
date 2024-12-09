case_configurations = [
        # RSPH model examples
        {
            "case" : 'RSPH/damBreak2D_RSPH'
        },
        # DBC model examples
        {
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "DBC",
            "h_coef" :  2**0.5
        },
        #{
        #    "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        #    "bc-type" : "DBC",
        #    "h_coef" :  2
        #},
        {
            "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
            "bc-type" : "MDBC",
            "h_coef" :  2**0.5
        },
        #{
        #    "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        #    "bc-type" : "MDBC",
        #    "h_coef" :  2
        #},
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "DBC"
        },
        {
            "case" : "WCSPH-DBC/damBreak3D_WCSPH-DBC",
            "bc-type" : "MDBC"
        }
]
