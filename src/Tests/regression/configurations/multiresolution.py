"""multiresolution regression configurations (new group)."""

configurations = [
    # multiresolution groups are structural smoke until H1/H2 freezes references (G21)
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution",
        "evaluation-function" : "empty",  # structural smoke until H1/H2 freezes references (G21)
        "target" : "damBreak2D_WCSPH-BI_MR_cuda",
        "init-script" : "init_multi_resolution.py"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution_local-timestepping",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI_multiresolution_local-timestepping",
        "evaluation-function" : "empty",  # structural smoke until H1/H2 freezes references (G21)
        "target" : "damBreak2D_WCSPH-BI_MR_LTS_cuda",
        "init-script" : "init_multi_resolution_lts.py"
    },
    {
        "case-tag" : "WCSPH-BI/damBreak3D_WCSPH-BI_multiresolution",
        "case" : "WCSPH-BI/damBreak3D_WCSPH-BI_multiresolution",
        "evaluation-function" : "empty",  # structural smoke until H1/H2 freezes references (G21)
        "target" : "damBreak3D_WCSPH-BI_MR_cuda",
        "init-script" : "init_multi_resolution.py"
    }
]
