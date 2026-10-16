"""The null: fixed-K, opportunity-proportional, WGD-blind redistribution of events.
"""
from __future__ import annotations

import numpy as np

def opportunity(burden):
    """Opportunity q_t ∝ the sample's full same-class burden B_t (including the locus itself; callers pass
    B_{-g} + x). Within a block, w_t = n * B_t / sum(B). Same functional form for every event class; only the
    input burden differs (mutation count, focal or broad burden)."""
    return np.clip(np.asarray(burden, float), 0.0, None)

def fixed_k_null(x, ww, cidx, q, seed, nsim, chunk=1000):
    """Exact fixed-K conditional null for one locus (WGD is a fixed anchor, not simulated).

    x : (n,) 0/1 observed events.
    ww : (n,) 0/1 anchor (WGD) per sample; observed and simulated statistics share this vector.
    cidx : list of per-cancer index arrays (K_c fixed within each block).
    q : (n,) opportunity ∝ burden; drives redistribution within each cancer.
    Row-sum-corrected redistribution: within a block E_t = K_c * q_t/sum(q), draw U~Uniform(0,1) and keep the
    K_c cells with the largest (E_t - U_t) (gene frequency exact, burden approximately preserved).
    Returns (S_obs, S_sim[array length nsim]).
    """
    Sobs = float((ww * x).sum())
    rng = np.random.RandomState(seed)
    parts = []
    for i in range(0, nsim, chunk):
        b = min(chunk, nsim - i)
        Sc = np.zeros(b)
        for idx in cidx:
            K = int(x[idx].sum())
            if K == 0:
                continue
            qc = q[idx]; wwc = ww[idx]; nc = len(idx)
            if K >= nc:
                Sc += wwc.sum(); continue
            E = K * (qc / qc.sum())                    # E_t = K_c * mu_t, sum_t E = K_c
            keys = E[None, :] - rng.rand(b, nc)         # row-sum-corrected redistribution: E - Uniform(0,1)
            top = np.argpartition(-keys, K - 1, axis=1)[:, :K]
            Sc += wwc[top].sum(1)
        parts.append(Sc)
    return Sobs, np.concatenate(parts)


