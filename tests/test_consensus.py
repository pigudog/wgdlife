"""Consensus meta-analysis contracts (timing.consensus)."""
import numpy as np
import pandas as pd

from wgdlife.timing.consensus import consensus, dl_random_effects, expit, type_cells


def test_dl_reduces_to_fixed_effect_when_homogeneous():
    y, v = np.array([0.5, 0.5, 0.5]), np.array([0.1, 0.1, 0.1])
    est, se, tau2, i2 = dl_random_effects(y, v)
    assert abs(est - 0.5) < 1e-12 and abs(se - np.sqrt(0.1 / 3)) < 1e-12 and tau2 == 0.0 and i2 == 0.0


def test_dl_single_cell_passthrough():
    est, se, tau2, i2 = dl_random_effects([0.3], [0.04])
    assert est == 0.3 and se == 0.2 and np.isnan(tau2) and np.isnan(i2)


def _patients(n_pre, n_post, platform, cancer, event, seed):
    rng = np.random.default_rng(seed)
    P = np.r_[rng.uniform(0.8, 1.0, n_pre), rng.uniform(0.0, 0.2, n_post)]
    return pd.DataFrame(dict(platform=platform, patient=[f"{platform}{i}" for i in range(len(P))], cancer=cancer,
                             event=event, P=P, w=1.0))


def test_consensus_recovers_the_majority_side_and_counts_strata():
    pat = pd.concat([_patients(40, 5, "MSK", "Breast", "TP53", 1), _patients(30, 4, "MSK", "Lung", "TP53", 2),
                     _patients(20, 3, "TCGA", "Breast", "TP53", 3), _patients(2, 30, "MSK", "Breast", "loss_3q", 4)])
    T = type_cells(pat)
    assert set(T.platform) == {"MSK", "TCGA"} and T.loc[T.event == "TP53", "n_pat"].sum() == 102
    C = consensus(pat).set_index("event")
    assert C.loc["TP53", "P_consensus"] > 0.7 and C.loc["TP53", "k_types"] == 3 and C.loc["TP53", "k_platforms"] == 2   # +.5 pseudo-counts and EB shrinkage pull it below the raw 0.88
    assert C.loc["TP53", "ci_lo"] > 0.5 and C.loc["loss_3q", "P_consensus"] < 0.2
    assert abs(expit(0.0) - 0.5) < 1e-12
