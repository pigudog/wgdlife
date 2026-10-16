"""Arm-level allele-state timing rule for chromosome-arm losses: per (sample, arm) the modal (major, minor)
configuration covering >= 80% of the informative arm length; (2,0) = before WGD, (2,1) = after WGD, otherwise
ambiguous. Intended for samples with a single WGD (filter upstream).
"""

import numpy as np
import pandas as pd

from wgdlife.cohort import cna_events6 as c6

def prep_segments(seg, samples, build="hg19"):
    """Arm-annotated segment table (cna_events6.prep_segments) in the geometry of ``build`` (hg19 | hg38). The
    switch is process-global (cna_events6.set_genome); every timing platform of the article is GRCh37, hence the
    default."""
    c6.set_genome(build)
    return c6.prep_segments(seg, samples)

ARMS = ["17p", "9p", "4q", "1q", "8q", "17q", "16q", "8p", "18q", "13q", "3p"]

def wilson(k, n, z=1.959964):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)

def arm_configs(s):
    """per (sample, arm): modal (major, minor) covering >= 80% of the lcn-informative arm length."""
    s = s.dropna(subset=["lcn"]).copy()
    s["major"] = (s.tcn - s.lcn).astype(int); s["minor"] = s.lcn.astype(int)
    g = s.groupby(["SAMPLE_ID", "arm", "major", "minor"])["len"].sum().reset_index()
    g = g.join(g.groupby(["SAMPLE_ID", "arm"])["len"].sum().rename("tot"), on=["SAMPLE_ID", "arm"])
    g = g.sort_values(["SAMPLE_ID", "arm", "len"], ascending=[True, True, False]).drop_duplicates(["SAMPLE_ID", "arm"])
    g = g[g.len / g.tot >= 0.8].copy()
    g["tcn"] = g.major + g.minor
    g["timing"] = np.where((g.major == 2) & (g.minor == 0), "pre", np.where((g.major == 2) & (g.minor == 1), "post", "amb"))
    g["loss"] = g.tcn < 4; g["loh"] = g.minor == 0; g["gain"] = g.tcn > 4
    return g

def arm_timing(g, platform):
    """frac_pre = n(2+0) / (n(2+0) + n(2+1)) per arm, with Wilson 95% CI."""
    rows = []
    for a in ARMS:
        x = g[g.arm == a]; k = int((x.timing == "pre").sum()); m = int((x.timing == "post").sum())
        lo, hi = wilson(k, k + m)
        rows.append(dict(platform=platform, arm=a, n_pre=k, n_post=m, frac_pre=k / (k + m) if k + m else np.nan, lo=lo, hi=hi))
    return pd.DataFrame(rows)
