"""Frequency-stratified empirical FDR for the WGD-anchored selection statistic: the shrunk normalized effect size
(nES) and its FDR (nFDR2), with the null effect sizes shrunk by the same per-locus mean as the observed ones.
"""
from __future__ import annotations
import numpy as np

SIN45 = np.sin(np.pi / 4)
FREQ_CATS = (0.0, 0.02, 0.05, 0.10, 1.0)

def estimate_fdr2(obs, exp, nsim, maxFDR=0.10):
    """Top-down empirical FDR within one frequency bin.
    obs : |nES| of observed features (1D).
    exp : |rES| pooled across the bin's features x simulations (1D null distribution).
    nsim: number of (retained) simulations. Returns per-feature FDR (features not reached in the
    top-down scan, i.e. past the first crossing of maxFDR, stay 1.0)."""
    obs = np.asarray(obs, float); exp = np.asarray(exp, float)
    no = len(obs); all_fdr = np.ones(no)
    if no == 0:
        return all_fdr
    order = np.argsort(-obs)                       # descending scan order
    vals = obs[order]
    obs_sorted_asc = np.sort(obs)
    exp_sorted_asc = np.sort(exp)
    ne = len(exp_sorted_asc)
    n_obs_ge = no - np.searchsorted(obs_sorted_asc, vals, side="left")   # #obs >= val
    n_exp_ge = ne - np.searchsorted(exp_sorted_asc, vals, side="left")   # #null >= val
    fdr_seq = np.minimum((n_exp_ge / nsim) / np.maximum(n_obs_ge, 1), 1.0)
    cross = np.where(fdr_seq >= maxFDR)[0]
    stop = cross[0] if len(cross) else no - 1      # assign up to & incl. first crossing, rest stay 1
    idx = order[:stop + 1]
    all_fdr[idx] = fdr_seq[:stop + 1]
    return all_fdr

def nes_and_fdr(wES, null_ES, freq, *, maxFDR=0.10, freq_cats=FREQ_CATS):
    """Compute nES + frequency-stratified empirical FDR for one null (one cohort x class).
      wES : (G,) observed effect size (O - E), scaled by sin(pi/4).
      null_ES : (G, nsim) per-simulation null ES (S_k - E), scaled by sin(pi/4).
      freq : (G,) locus frequency (fraction of samples altered), the stratifier.
    Returns dict(nES, nFDR2, mean_absES). Significance := nFDR2 <= maxFDR."""
    wES = np.asarray(wES, float) * SIN45; null_ES = np.asarray(null_ES, float) * SIN45; freq = np.asarray(freq, float)
    mean_absES = np.abs(null_ES).mean(1)                              # mean_k |rES_k| per locus
    nES = np.sign(wES) * np.clip(np.abs(wES) - mean_absES, 0.0, None) # shrinkage
    nsim = null_ES.shape[1]
    nFDR2 = np.ones(len(wES))
    for i in range(1, len(freq_cats)):
        sel = (freq > freq_cats[i - 1]) & (freq <= freq_cats[i])
        if sel.sum() > 1:
            obs = np.abs(nES[sel])
            # null ES shrunk by the same per-locus mean |rES| as the observed
            exp = np.clip(np.abs(null_ES[sel]) - mean_absES[sel][:, None], 0.0, None).ravel()
            nFDR2[sel] = estimate_fdr2(obs, exp, nsim, maxFDR=maxFDR)
    return {"nES": nES, "nFDR2": nFDR2, "mean_absES": mean_absES}

def call_type(nES):
    """CO (co-occurrence with WGD) if nES > 0, ME (mutual exclusivity with WGD) if nES < 0."""
    return np.where(nES > 0, "CO", np.where(nES < 0, "ME", "none"))
