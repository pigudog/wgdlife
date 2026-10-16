"""Selection engine on synthetic data: planted co-occurrence and exclusivity are found and null genes are not; seeds
are process-independent; the null is burden-neutral; nFDR2 compares shrunk observed and null effect sizes.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from wgdlife.selection.anchored import anchored_locus, significance
from wgdlife.selection.fdr import SIN45, estimate_fdr2, nes_and_fdr
from wgdlife.selection.inference import strata

def _synthetic(n=600, seed=7):
    rng = np.random.default_rng(seed)
    wgd = (rng.random(n) < 0.35).astype(float)
    cancer = np.array(["typeA"] * (n // 2) + ["typeB"] * (n - n // 2))
    X = np.zeros((n, 4))
    X[:, 0] = rng.random(n) < np.where(wgd == 1, 0.40, 0.10)     # enriched in WGD+
    X[:, 1] = rng.random(n) < np.where(wgd == 1, 0.05, 0.30)     # depleted in WGD+
    X[:, 2] = rng.random(n) < 0.15                                # null
    X[:, 3] = rng.random(n) < 0.20                                # null
    return X, wgd, cancer

def _run(seed_data, seed_engine, nsim=2000):
    X, wgd, cancer = _synthetic(seed=seed_data)
    cidx = strata(cancer, len(wgd)); bur = X.sum(1)
    rows = [anchored_locus(X[:, g], wgd, cidx, bur, seed_engine + g, nsim) for g in range(X.shape[1])]
    O = np.array([r[0] for r in rows]); E = np.array([r[1] for r in rows])
    p2 = np.array([min(1.0, 2 * min(r[3], r[4])) for r in rows])
    sg = significance(O, E, X.mean(0), np.vstack([r[5] for r in rows]))
    return O / E, p2, sg

def test_planted_enrichment_and_depletion():
    ER, p2, sg = _run(7, 0)
    assert ER[0] > 1.2 and p2[0] < 0.05 and sg["type"][0] == "CO"
    assert ER[1] < 0.8 and p2[1] < 0.05 and sg["type"][1] == "ME"
    assert np.isfinite(ER).all() and np.isfinite(p2).all() and sg["nFDR2"].shape == (4,)

def test_null_genes_not_called():
    ER, p2, sg = _run(11, 1)
    assert p2[2] > 0.05 and p2[3] > 0.05

def test_stable_seed_is_process_independent():
    """Null seeds must not depend on PYTHONHASHSEED (hash() is salted per process)."""
    code = ("import numpy as np; from wgdlife.selection.anchored import stable_seed, anchored_locus;"
            "r=np.random.RandomState(0); n=400; x=(r.rand(n)<0.2).astype(float); w=(r.rand(n)<0.4).astype(float);"
            "c=[np.arange(n)]; b=r.gamma(2,1,n);"
            "s=stable_seed('cohort','broad_gain')*1000+7;"
            "o=anchored_locus(x,w,c,b,s,200);"
            "print(s, float(o[0]), float(o[1]))")
    outs = []
    for hs in ("1", "12345"):
        env = dict(os.environ, PYTHONHASHSEED=hs, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
        outs.append(subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True,
                                   check=True).stdout.strip())
    assert outs[0] == outs[1], outs

def test_null_is_burden_neutral():
    """With no locus-specific WGD effect and only higher burden in WGD+, the null must not call CO."""
    rng = np.random.default_rng(1); n, G = 2000, 30
    wgd = (rng.random(n) < 0.25).astype(float)
    B = np.minimum(rng.poisson(np.where(wgd == 1, 12.0, 3.0)), G)
    p = rng.gamma(2, 1, G); p /= p.sum()
    X = np.zeros((n, G))
    for i in range(n):
        if B[i] > 0:
            X[i, rng.choice(G, size=B[i], replace=False, p=p)] = 1
    bur = X.sum(1)
    O = E = 0.0
    for g in range(G):
        o, e, *_ = anchored_locus(X[:, g], wgd, [np.arange(n)], bur, 1000 + g, 200)
        O += o; E += e
    assert abs(O / E - 1) < 0.03, O / E

def test_nfdr2_shrinks_null_like_observed():
    """nFDR2 must compare shrunk |nES| with the equally shrunk null."""
    rng = np.random.default_rng(3); G, K = 60, 300
    null = rng.normal(0, 1, (G, K)); obs = rng.normal(0, 1, G); obs[:8] += 6
    freq = np.full(G, 0.03)
    r = nes_and_fdr(obs, null, freq)
    m = np.abs(null * SIN45).mean(1)
    exp = np.clip(np.abs(null * SIN45) - m[:, None], 0, None).ravel()
    ref = estimate_fdr2(np.abs(r["nES"]), exp, K)
    assert np.allclose(r["nFDR2"], ref)
