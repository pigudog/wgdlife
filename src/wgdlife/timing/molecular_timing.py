"""Timing of panel mutations relative to WGD from the mutant copy number on allele-specific copy-number states.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

def panel_classify(df, *, mult_col="expected_alt_copies", tcn_col="tcn", lcn_col="lcn",
                   clon_col="clonality", pre_thresh=2.0):
    """Pre/post-WGD call per mutation on allele-specific-CN-annotated panel data (after Bielski et al. 2018), using
    the pipeline's own mutant-copy column (``expected_alt_copies``, an integer). The caller restricts to WGD+ samples
    first (timing is relative to a truncal doubling).

    Classes:
      pre_wgd balanced doubled state (major == minor >= 2) or 2+0, mult >= pre_thresh
      post_wgd balanced doubled state or 2+0, mult < pre_thresh
      pre_wgd_unbal unbalanced state with major >= 2, mult >= pre_thresh (one-sided call; reported only)
      ambiguous unbalanced state with major >= 2, mult < pre_thresh (a pre-WGD mutation that lost one copy
                     and a post-WGD mutation reach the same endpoint)
      subclonal SUBCLONAL calls (after the most recent common ancestor, hence after a truncal WGD)
      uninformative everything else (major < 2, missing values, INDETERMINATE)"""
    d = df.copy()
    tcn = pd.to_numeric(d[tcn_col], errors="coerce")
    lcn = pd.to_numeric(d[lcn_col], errors="coerce")
    mult = pd.to_numeric(d[mult_col], errors="coerce")
    major = tcn - lcn
    d["major"] = major
    clonal = d[clon_col].astype(str).str.upper().eq("CLONAL")
    subcl = d[clon_col].astype(str).str.upper().eq("SUBCLONAL")
    inf = clonal & major.ge(2) & mult.notna() & lcn.notna()
    # states where pre and post are parsimony-separable: balanced doubled states (2+2, 3+3, ...) and the
    # cnLOH-doubled 2+0 state (mult >= 2 = mutation and LOH before WGD; mult <= 1 = post-WGD mutation)
    bal = inf & ((major == lcn) | ((major == 2) & (lcn == 0)))
    unb = inf & ~bal
    pre = mult >= pre_thresh
    d["timing"] = np.select(
        [bal & pre, bal & ~pre, unb & pre, unb & ~pre, subcl],
        ["pre_wgd", "post_wgd", "pre_wgd_unbal", "ambiguous", "subclonal"],
        default="uninformative")
    # timed CN state (e.g. 2+2, 3+3, 2+0): pre/post base rates differ by state, so downstream tests stratify by it
    st = major.astype("Int64").astype(str) + "+" + lcn.astype("Int64").astype(str)
    d["timing_state"] = np.where(bal, st, None)
    return d
