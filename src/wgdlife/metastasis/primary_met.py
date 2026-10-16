"""Primary-versus-metastasis comparison of WGD-associated events: cancer-type- and size-matched resampling of
primary and metastatic samples, an admixed null, and the category rule.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from wgdlife.selection.anchored import stable_seed
from wgdlife.selection.fdr import SIN45

RUNS, FRAC, MINCLASS, SEED = 10, 0.8, 31, 7
SIG_NFDR2 = 0.10
CAT_NSIG_HI, CAT_NSIG_LO, NULL_PCTL = 6, 2, 99

def per_cancer_k(Pby, Mby, frac=FRAC, minclass=MINCLASS):
    """{matching unit -> per-run draw size}: a unit (cancer type) is kept when it holds >= minclass primary (P) and
    >= minclass metastatic (M) samples; k = round(frac * the smaller arm). Unit order follows Pby."""
    out = {}
    for u in Pby:
        n = min(len(Pby[u]), len(Mby.get(u, ())))
        if n >= minclass:
            out[u] = round(frac * n)
    return out

def draw_idsets(Pby, Mby, n_ca, posn, seed=SEED, runs=RUNS):
    """{(group, run): row-index array into the pooled sample order} for the groups P and M (k samples per unit) and
    the admixed controls A and B (k // 2 from P and the remainder from M per unit). One default_rng(seed) is consumed
    group-major (P, M, A, B), run-major, unit-major, so every draw is reproducible."""
    rng = np.random.default_rng(seed)
    ids = {}
    for g in ("P", "M", "A", "B"):
        for r in range(runs):
            out = []
            for u, k in n_ca.items():
                if g == "P":
                    out += list(rng.choice(Pby[u], k, replace=False))
                elif g == "M":
                    out += list(rng.choice(Mby[u], k, replace=False))
                else:
                    kp = k // 2
                    out += list(rng.choice(Pby[u], kp, replace=False))
                    out += list(rng.choice(Mby[u], k - kp, replace=False))
            ids[(g, r)] = np.array([posn[s] for s in out])
    return ids

def locus_seed(gtag, run, modality, g):
    """Deterministic per-locus null seed: stable_seed(group, run, modality) * 1000 + locus index."""
    return stable_seed(gtag, run, modality) * 1000 + int(g)

def is_sig(df):
    return df["nFDR2"] <= SIG_NFDR2

def _testable(runs):
    """Per-row testability: nFDR2 is estimated within frequency bins and a bin holding a single locus is
    skipped (nFDR2 stays 1.0 = 'not reached', not 'tested and null'; see selection.fdr.nes_and_fdr). Runs
    tables carry a `testable` column (= n_in_bin > 1); tables without it treat every row as testable."""
    if "testable" in runs.columns:
        return runs["testable"].astype(bool)
    return pd.Series(True, index=runs.index)

def per_group_stats(runs, group, sig_only=False, es=None, cnt=None, tcnt=None, nes=None):
    """Per-locus statistics of one resample group (the shared body of `aggregate`).
    Returns (modality, gene, <es>, <cnt>, <nes>, <tcnt>): <es> = mean wES over the rows kept (ALL runs the locus
    appears in, or only the significant ones when sig_only), <cnt> = # significant
    AND testable runs, <tcnt> = (runs of this group) − (runs where the locus sat alone in its frequency bin)."""
    es = es or f"es_{group}"
    cnt = cnt or f"{group}_count"
    tcnt = tcnt or f"{group}_test"
    nes = nes or f"n{es}"
    tst = _testable(runs)
    r = runs.assign(_test=tst.astype(int), _untest=(~tst).astype(int), _sig=(is_sig(runs) & tst).astype(int))
    x = r[r.group == group]
    nrun = int(x.run.nunique())
    if sig_only:
        x = x[x._sig == 1]
    s = x.groupby(["modality", "gene"]).agg(
        **{es: ("wES", "mean"), cnt: ("_sig", "sum"), "_u": ("_untest", "sum"), nes: ("nES", "mean")}).reset_index()
    s[tcnt] = nrun - s.pop("_u")
    return s

def aggregate(runs):
    """runs (group, run, modality, gene, wES, nES, p, nFDR2[, n_in_bin, testable]) -> per-group mean
    wES (over ALL runs the locus appears in), # significant runs (p_count/m_count: testable runs only) and
    # testable runs (p_test/m_test = RUNS minus the runs where the locus sat alone in its frequency bin; a
    locus absent from a run — no event in that draw — counts as tested-and-null, so
    tables without `testable` give p_test = m_test = RUNS). Returns (obs, nul) frames for `categorize`.
    Null: per admixed run keep only the loci significant in that run; es_A (es_B) = mean wES over the A (B)
    runs where significant; loci significant in >= 1 A run AND >= 1 B run (inner merge)."""
    def agg(g, es, cnt, tcnt, sig_only=False):
        return per_group_stats(runs, g, sig_only, es, cnt, tcnt)

    P, M = agg("P", "es_P", "p_count", "p_test"), agg("M", "es_M", "m_count", "m_test")
    obs = P.merge(M, on=["modality", "gene"], how="outer").fillna({"p_count": 0, "m_count": 0, "p_test": 0, "m_test": 0})
    obs = obs[obs[["es_P", "es_M"]].notna().all(1)].copy()
    obs["dist_obs"] = (obs.es_M - obs.es_P) * SIN45
    obs = obs[(obs.p_count > 0) | (obs.m_count > 0)].copy()
    A, B = agg("A", "es_A", "a_sig", "a_test", True), agg("B", "es_B", "b_sig", "b_test", True)
    nul = A.merge(B, on=["modality", "gene"], how="inner")
    nul = nul[(nul.a_sig > 0) | (nul.b_sig > 0)].copy()
    nul["dist"] = (nul.es_B - nul.es_A) * SIN45
    return obs, nul

def null_threshold(nul, pctl=NULL_PCTL):
    """|dist| threshold from the admix null: the order statistic sort(abs(distance))[round(pctl / 100 * n)]
    (round = half-to-even)."""
    if not len(nul):
        return np.inf
    d = np.sort(nul.dist.abs().to_numpy())
    return float(d[max(int(round(pctl / 100 * len(d))), 1) - 1])

def categorize(obs, nul, hi=CAT_NSIG_HI, lo=CAT_NSIG_LO, pctl=NULL_PCTL):
    """Category rule, nested inside the distance gate (|dist| > thr, else "none"). Returns (obs with `cat`, threshold).
    Enriched by MAGNITUDE (Menr if |es_M| > |es_P|, Penr if |es_P| > |es_M|), overridden in order by
    Pspec (p >= hi & m <= lo), Mspec (m >= hi & p <= lo), switch (both >= hi, opposite sign).
    Testable-aware: a side-specific call needs the OTHER side testable in >= hi runs (p_test / m_test from
    `aggregate`; a locus alone in its frequency bin has nFDR2 = 1 by construction, which is 'not reached', not
    'null'); otherwise the locus stays enriched. obs without p_test/m_test = every run testable."""
    thr = null_threshold(nul, pctl)
    if "p_test" not in obs.columns:
        obs = obs.assign(p_test=np.inf, m_test=np.inf)

    def cat(r):
        if not abs(r.dist_obs) > thr:
            return "none"
        c = "Menr" if abs(r.es_M) > abs(r.es_P) else "Penr"
        if r.p_count >= hi and r.m_count <= lo and r.m_test >= hi:
            c = "Pspec"
        if r.m_count >= hi and r.p_count <= lo and r.p_test >= hi:
            c = "Mspec"
        if r.p_count >= hi and r.m_count >= hi and np.sign(r.es_P * r.es_M) == -1:
            c = "switch"
        return c

    obs = obs.copy()
    obs["cat"] = obs.apply(cat, axis=1) if len(obs) else []
    return obs, thr
