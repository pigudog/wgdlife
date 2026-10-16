"""wgdlife.timing.consensus — consensus timing of an event relative to WGD across platforms and cancer types.

Three-level hierarchy:
  patient         soft P(event before WGD) with resolvable mass w = P_pre + P_post ("unknown" is dropped as weight,
                  never hard-thresholded); multi-sample patients averaged upstream.
  cancer type     within platform: effective counts S = sum(w P), N = sum(w) -> logit((S + .5)/(N + 1)),
                  var = 1/(S + .5) + 1/(N - S + .5); empirical-Bayes normal-normal shrinkage across types within
                  platform x event; patient heterogeneity = expected minority share.
  consensus       DerSimonian-Laird random effects over TYPE cells (never over raw pooled patients); I2 over types;
                  cross-platform I2 over the per-platform random-effects estimates.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

MIN_PAT_TYPE = 5


def expit(x):
    return 1 / (1 + np.exp(-np.asarray(x, float)))


def dl_random_effects(y, v):
    """DerSimonian-Laird random effects. Returns (estimate, se, tau2, I2); tau2 and I2 are nan for k < 2."""
    y, v = np.asarray(y, float), np.asarray(v, float)
    w = 1 / v
    if len(y) == 1:
        return y[0], np.sqrt(v[0]), np.nan, np.nan
    fe = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - fe) ** 2); dfq = len(y) - 1
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - dfq) / c) if c > 0 else 0.0
    ws = 1 / (v + tau2)
    est = np.sum(ws * y) / np.sum(ws)
    I2 = max(0.0, (Q - dfq) / Q) if Q > 0 else 0.0
    return est, np.sqrt(1 / np.sum(ws)), tau2, I2


def type_cells(pat: pd.DataFrame, min_pat: int = MIN_PAT_TYPE) -> pd.DataFrame:
    """Cancer-type level: one cell per (platform, cancer, event) from patient rows [platform, patient, cancer, event, P, w]."""
    rows = []
    for (pf, ca, ev), g in pat.groupby(["platform", "cancer", "event"]):
        if len(g) < min_pat:
            continue
        S = float((g.w * g.P).sum()); N = float(g.w.sum())
        p = (S + .5) / (N + 1)
        minority = (float((g.w * (1 - g.P)).sum()) if p > .5 else S) / N
        rows.append(dict(platform=pf, cancer=ca, event=ev, n_pat=len(g), N_eff=N, P_type=p,
                         y=np.log(p / (1 - p)), v=1 / (S + .5) + 1 / (N - S + .5), minority_share=minority))
    T = pd.DataFrame(rows)
    if not len(T):
        return T
    T["P_type_eb"] = np.nan
    for (pf, ev), g in T.groupby(["platform", "event"]):          # EB shrinkage within platform x event
        mu, se, tau2, _ = dl_random_effects(g.y, g.v)
        tau2 = 0.0 if not np.isfinite(tau2) else tau2
        yeb = (g.y / g.v + mu / tau2) / (1 / g.v + 1 / tau2) if tau2 > 0 else pd.Series(mu, index=g.index)
        T.loc[g.index, "P_type_eb"] = expit(yeb)
    return T


def consensus(pat: pd.DataFrame, min_pat: int = MIN_PAT_TYPE) -> pd.DataFrame:
    """Consensus level: one row per event — P_consensus with 95% CI, I2 over types and over platforms, per-platform estimates."""
    T = type_cells(pat, min_pat)
    rows = []
    for ev, g in T.groupby("event"):
        est, se, tau2, I2 = dl_random_effects(g.y, g.v)
        r = dict(event=ev, k_types=len(g), k_platforms=g.platform.nunique(), n_pat=int(g.n_pat.sum()),
                 P_consensus=float(expit(est)), ci_lo=float(expit(est - 1.96 * se)), ci_hi=float(expit(est + 1.96 * se)),
                 I2_types=I2)
        pe, pv = [], []
        for pf, gp in g.groupby("platform"):
            e2, s2, _, i2 = dl_random_effects(gp.y, gp.v)
            r[f"P_{pf}"] = float(expit(e2)); r[f"k_types_{pf}"] = len(gp); r[f"I2_types_{pf}"] = i2
            pe.append(e2); pv.append(s2 ** 2)
        r["I2_platforms"] = dl_random_effects(pe, pv)[3] if len(pe) > 1 else np.nan
        side_c = expit(est) > .5
        r["n_types_pre"] = int((g.P_type > .5).sum()); r["n_types_post"] = int((g.P_type <= .5).sum())
        r["frac_types_minority_side"] = float(((g.P_type > .5) != side_c).mean())
        W = pat[pat.event == ev]
        mass_min = (W.w * (1 - W.P)).sum() if side_c else (W.w * W.P).sum()
        r["patient_minority_share"] = float(mass_min / W.w.sum()) if W.w.sum() > 0 else np.nan
        rows.append(r)
    return pd.DataFrame(rows)
