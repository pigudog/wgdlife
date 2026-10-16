"""primary_met contracts: constants, draw order, locus seeds, category rule, aggregation (data-independent).
"""

import numpy as np
import pandas as pd

from wgdlife.metastasis import primary_met as PM


def test_pm_constants():
    assert (PM.RUNS, PM.FRAC, PM.MINCLASS, PM.SEED) == (10, 0.8, 31, 7)
    assert (PM.CAT_NSIG_HI, PM.CAT_NSIG_LO, PM.NULL_PCTL) == (6, 2, 99)
    assert PM.SIG_NFDR2 == 0.10


def test_draw_idsets_reproduces_draw_order():
    Pby = {"a": [f"p{i}" for i in range(40)], "b": [f"q{i}" for i in range(30)]}
    Mby = {"a": [f"m{i}" for i in range(50)], "b": [f"n{i}" for i in range(10)]}
    n_ca = PM.per_cancer_k(Pby, Mby)
    assert n_ca == {"a": 32}                       # b skipped: min(30, 10) < 31
    pool = [s for ca in n_ca for s in Pby[ca] + Mby[ca]]
    posn = {s: i for i, s in enumerate(pool)}
    ids = PM.draw_idsets(Pby, Mby, n_ca, posn)
    rng = np.random.default_rng(7)
    for g in ("P", "M", "A", "B"):
        for r in range(10):
            ref = []
            for ca, k in n_ca.items():
                if g in ("P", "M"):
                    ref += [posn[s] for s in rng.choice({"P": Pby, "M": Mby}[g][ca], k, replace=False)]
                else:                                  # admixed control: half P, half M, remainder to M
                    kp = k // 2
                    ref += [posn[s] for s in rng.choice(Pby[ca], kp, replace=False)]
                    ref += [posn[s] for s in rng.choice(Mby[ca], k - kp, replace=False)]
            assert np.array_equal(ids[(g, r)], np.array(ref))


def test_locus_seed_deterministic():
    assert PM.locus_seed("P", 3, "snv", 12) == PM.locus_seed("P", 3, "snv", 12)
    assert PM.locus_seed("P", 3, "snv", 12) != PM.locus_seed("M", 3, "snv", 12)


def test_categorize_rule():
    obs = pd.DataFrame({"modality": ["x"] * 7, "gene": list("abcdefg"),
                        "es_P": [10, 1, 10, -5, 1, 50, 10], "es_M": [-10, 10, 1, -5, 1, 1, 12],
                        "p_count": [8, 1, 9, 0, 3, 3, 1], "m_count": [8, 9, 0, 0, 3, 3, 9]})
    obs["dist_obs"] = (obs.es_M - obs.es_P) * PM.SIN45
    nul = pd.DataFrame({"dist": np.linspace(-5, 5, 101)})
    out, thr = PM.categorize(obs, nul)
    assert np.isclose(thr, 5.0)                    # order statistic sort(|dist|)[round(0.99 * 101) - 1]
    assert list(out["cat"]) == ["switch", "Mspec", "Pspec", "none", "none", "Penr", "none"]


def test_aggregate_counts_significant_testable_runs():
    rows = []
    for g in "PMAB":
        for r in range(10):
            rows.append((g, r, "snv", "TP53", 5.0, 1.0, 0.001, 0.05 if r < 4 else 0.5, 3, r != 0))
    runs = pd.DataFrame(rows, columns=["group", "run", "modality", "gene", "wES", "nES", "p", "nFDR2", "n_in_bin", "testable"])
    obs, nul = PM.aggregate(runs)
    assert int(obs.p_count.iloc[0]) == 3 and int(obs.p_test.iloc[0]) == 9     # run 0 is significant but not testable
    assert int(obs.m_count.iloc[0]) == 3 and len(nul) == 1
    obs2, _ = PM.aggregate(runs.drop(columns="testable"))
    assert int(obs2.p_count.iloc[0]) == 4 and int(obs2.p_test.iloc[0]) == 10
