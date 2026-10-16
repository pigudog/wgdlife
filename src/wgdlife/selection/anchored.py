"""WGD-anchored selection: the per-locus null and the significance call from the retained null residuals.
"""
from __future__ import annotations
import numpy as np, zlib
from .null import opportunity, fixed_k_null
from .inference import empirical_p
from . import fdr

def stable_seed(*parts) -> int:
    """16-bit null-seed key for (cohort, class, ...), identical in every process.

    Python's hash() on str / tuples is salted per process (PYTHONHASHSEED); crc32 gives the same seed in every
    process. Callers use ``stable_seed(...) * 1000 + locus_index``.
    """
    return zlib.crc32("|".join(map(str, parts)).encode()) & 0xffff

def anchored_locus(x, wgd, cidx, burden, seed, nsim):
    """WGD-anchored null for ONE locus.
      x : (n,) 0/1 observed events for the locus.
      wgd : (n,) 0/1 anchor. cidx: list of index arrays (cancer x panel strata).
      burden : (n,) same-class per-sample event burden (opportunity source).
    Returns (O, E, sd, p_enrich, p_deplete, resid) with resid = (S_k - E) as float32 (len nsim)."""
    x = np.asarray(x, float)
    Bmg = np.clip(np.asarray(burden, float) - x, 0.0, None)   # leave-one-out background burden
    q = opportunity(Bmg + x)
    ww = np.asarray(wgd, float)
    Sobs, S = fixed_k_null(x, ww, cidx, q, seed, nsim)
    S = np.asarray(S, float)
    Em = max(float(S.mean()), 1e-9)
    pe, pd_ = empirical_p(Sobs, S)
    return float(Sobs), Em, float(S.std()), float(pe), float(pd_), (S - Em).astype(np.float32)

def significance(O, E, freq, null_ES, *, maxFDR=0.10, freq_cats=fdr.FREQ_CATS):
    """Significance from the retained sufficient statistics (no re-simulation).
      O, E : (G,) observed / null-mean counts. freq: (G,) locus frequency.
      null_ES: (G, nsim) retained per-sim residuals (S_k - E).
    Returns dict(wES, nES, nFDR2, mean_absES, sig, type). sig := nFDR2 <= maxFDR."""
    O = np.asarray(O, float); E = np.asarray(E, float)
    wES = O - E
    r = fdr.nes_and_fdr(wES, null_ES, freq, maxFDR=maxFDR, freq_cats=freq_cats)
    return {"wES": wES, "nES": r["nES"], "nFDR2": r["nFDR2"], "mean_absES": r["mean_absES"],
            "sig": r["nFDR2"] <= maxFDR, "type": fdr.call_type(r["nES"])}
