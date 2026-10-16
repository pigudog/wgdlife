"""The WGD criterion: length-weighted fraction of the autosomal genome with major copy number >= 2; WGD = fraction > 0.5.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

AUTOSOMES = frozenset(str(i) for i in range(1, 23))

def mcn_ge2_fraction_from_segments(seg: pd.DataFrame, *, chrom_col: str = "Chromosome", start_col: str = "Start",
                                   end_col: str = "End", major_col: str = "Major_Copy_Number",
                                   sample_col: str | None = None, autosomes=AUTOSOMES):
    """Length-weighted fraction of the autosomal genome with major copy number >= 2 (allele-specific segment
    profiles: TCGA ASCAT from the GDC, PCAWG consensus, OncoScan ASCAT); WGD = fraction > 0.5.

    seg : one segment per row with chromosome, start, end and major copy number (GDC ASCAT column names by default;
          "chr" prefixes are tolerated). With ``sample_col`` the frame may hold many samples and a Series is returned.
    """
    s = seg[[c for c in (sample_col, chrom_col, start_col, end_col, major_col) if c]].copy()
    s["_c"] = s[chrom_col].astype(str).str.replace("chr", "", regex=False)
    s = s[s["_c"].isin(autosomes)]
    s["_len"] = (pd.to_numeric(s[end_col]) - pd.to_numeric(s[start_col])).clip(lower=0).astype(float)
    s["_hit"] = s["_len"].where(pd.to_numeric(s[major_col], errors="coerce") >= 2, 0.0)
    if sample_col is None:
        tot = s["_len"].sum()
        return float(s["_hit"].sum() / tot) if tot > 0 else float("nan")
    g = s.groupby(sample_col)[["_hit", "_len"]].sum()
    return (g["_hit"] / g["_len"].replace(0, np.nan)).rename("mcn_ge2_frac")

def call_wgd_from_fraction(frac, thresh: float = 0.5):
    """WGD status from the autosomal major-CN >= 2 fraction (> 0.5 = WGD)."""
    return (pd.Series(frac) > thresh).astype(int) if not np.isscalar(frac) else int(frac > thresh)
