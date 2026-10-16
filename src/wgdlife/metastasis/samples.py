"""Primary-site versus metastasis-site sample classification from a clinical sample table.
"""
import numpy as np
import pandas as pd

PRIMARY = "Primary"
METASTASIS = "Metastasis"
NOT_APPLICABLE = "Not Applicable"

def sample_type(cs: pd.DataFrame, col: str = "METASTATIC_SITE") -> pd.Series:
    """Per-sample tissue type from the clinical-sample table: "Primary" (METASTATIC_SITE == "Not Applicable"),
    "Metastasis" (any named site), NaN when METASTATIC_SITE is missing. Index and order follow `cs`."""
    site = cs[col]
    out = pd.Series(np.where(site.astype(str).str.strip().eq(NOT_APPLICABLE), PRIMARY, METASTASIS), index=cs.index, dtype=object)
    out[site.isna() | site.astype(str).str.strip().eq("")] = np.nan
    return out

