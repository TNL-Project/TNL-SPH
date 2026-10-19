"""qa_smoke regression configurations (new group)."""

configurations = [
    # pinned local smoke set: the two default damBreak2D cases
    {
        "case-tag" : "WCSPH-DBC/damBreak2D_WCSPH-DBC_default",
        "case" : "WCSPH-DBC/damBreak2D_WCSPH-DBC",
        "evaluation-function" : "damBreak2D_WCSPH_DBC"  # the one seeded gated reference
    },
    {
        "case-tag" : "WCSPH-BI/damBreak2D_WCSPH-BI_default",
        "case" : "WCSPH-BI/damBreak2D_WCSPH-BI",
        "evaluation-function" : "empty"  # BI gated ref does not exist yet
    }
]
