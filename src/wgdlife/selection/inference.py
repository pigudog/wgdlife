"""Per-locus inference: the empirical P values of the WGD-anchored statistic and the stratum index arrays.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

def empirical_p(S_obs, S_sim):
    """Empirical one-sided P values (enrichment, depletion) of the observed statistic against the null draws."""
    n = len(S_sim)
    pe = (1 + int((S_sim >= S_obs).sum())) / (n + 1)
    pd_ = (1 + int((S_sim <= S_obs).sum())) / (n + 1)
    return pe, pd_

def strata(cancer, n):
    """Index arrays of the strata (one per distinct label; a single block when cancer is None)."""
    if cancer is None:
        return [np.arange(n)]
    cancer = np.asarray(cancer)
    return [np.where(cancer == c)[0] for c in pd.unique(cancer)]
