"""Example: WGD-associated events in primary-site versus metastasis-site samples of one MSK-IMPACT cancer type.
"""
import argparse
import multiprocessing as mp
import os
import time
from collections import defaultdict
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import pandas as pd

from wgdlife.cohort import caller_wgd as W
from wgdlife.cohort import cna_events6 as c6
from wgdlife.cohort.cohorts import snv_bundle
from wgdlife.cohort.constants import CANCER_TYPES_12
from wgdlife.metastasis import primary_met as PM
from wgdlife.metastasis.samples import sample_type
from wgdlife.selection.anchored import anchored_locus, significance
from wgdlife.selection.fdr import FREQ_CATS
from wgdlife.selection.inference import strata

HERE = Path(__file__).resolve().parent
CHRLEN, CEN = c6.GENOME_CONSTANTS["hg19"]   # MSK-IMPACT FACETS coordinates are GRCh37
NE_ONCOTREE = {"SCBC", "UNEC", "GINETES", "STSC", "HNNE", "ALUCA", "LUNE", "LUCA", "LNET", "SCCO", "PANET", "PRNE", "PRSCC", "SCLC", "NEC", "NET"}
HYPER_FC, MINCLASS, MAXFDR = 10.0, 31, 0.10
CLASSES = ["broad_gain", "broad_loh", "focal_amp", "focal_del", "focal_loh"]
MODS = CLASSES + ["snv"]
RUN_COLS = ["group", "run", "modality", "gene", "freq", "O", "E", "wES", "nES", "p", "nFDR2", "n_in_bin", "testable"]
G = {}

def _worker(task):
    gtag, r, m, g = task
    ids, cidx, cidx_snv = G["IDSETS"][(gtag, r)]
    x = G["X"][m][ids, g]
    O, E, sd, pe, pdl, resid = anchored_locus(x, G["WGD"][ids], cidx_snv if m == "snv" else cidx, G["BUR"][m][ids],
                                              PM.locus_seed(gtag, r, m, g), G["NSIM"])
    return gtag, r, m, g, float(x.mean()), O, E, pe, pdl, resid

def bin_sizes(freq):
    b = np.digitize(np.asarray(freq, float), FREQ_CATS[1:-1], right=True)
    return pd.Series(b).map(pd.Series(b).value_counts()).to_numpy()

def locus_label(gene_bed):
    """gene -> 10-Mb arm bin label from the bed's own coordinates (e.g. '17p:0-10Mb')."""
    lab = {}
    for g, ch, mid in zip(gene_bed.gene, gene_bed.chrom.astype(str), gene_bed.mid):
        if ch in CEN and np.isfinite(mid):
            arm = "p" if mid < CEN[ch] else "q"; b = int(mid // 10e6) * 10
            lab[g] = f"{ch}{arm}:{b}-{b + 10}Mb"
    return lab

def is_terminal(unit):
    import re
    m = re.match(r"^(\d+)([pq]):(\d+)-(\d+)Mb$", str(unit))
    if not m:
        return False
    ch, arm, b0, b1 = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
    return (arm == "p" and b0 == 0) or (arm == "q" and b1 * 1e6 >= CHRLEN[ch])

def collapse_focal_loh(runs, gene_bed):
    """focal LOH genes -> 10-Mb arm loci (representative = largest O per group / run / locus); telomere-proximal loci dropped."""
    L = locus_label(gene_bed)
    fl = runs[runs.modality == "focal_loh"].copy(); fl["locus"] = fl.gene.map(L)
    fl = fl[fl.locus.notna() & ~fl.locus.map(is_terminal)]
    fl = fl.sort_values("O", ascending=False).groupby(["group", "run", "locus"], sort=False).head(1)
    fl["gene"] = fl.locus + " (" + fl.gene + ")"
    return pd.concat([runs[runs.modality != "focal_loh"], fl.drop(columns="locus")], ignore_index=True)

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cancer-type", default="Breast Cancer")
    ap.add_argument("--runs", type=int, default=PM.RUNS); ap.add_argument("--nsim", type=int, default=1000)
    ap.add_argument("--njob", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--data", default=os.environ.get("WGDLIFE_DATA", "data")); ap.add_argument("--out", default="results")
    a = ap.parse_args()
    c6.set_genome("hg19")   # MSK-IMPACT FACETS segments are GRCh37
    if a.cancer_type not in CANCER_TYPES_12:
        raise SystemExit(f"--cancer-type must be one of {list(CANCER_TYPES_12)}")
    d = Path(a.data) / "msk_impact_50k" / "derived"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tag = a.cancer_type.replace(" ", "_")
    t = time.time()
    S = pd.read_csv(d / "samples.tsv", sep="\t", low_memory=False)
    base = S[S.CANCER_TYPE.eq(a.cancer_type) & S.MSI_TYPE.eq("Stable") & S.wgd.astype(str).str.lower().isin(["true", "false"])
             & S.DISEASE_STATUS.isin(["Primary", "Metastasis"]) & ~S.ONCOTREE_CODE.isin(NE_ONCOTREE)].copy()
    base["arm"] = sample_type(base)                      # primary-site / metastasis-site by METASTATIC_SITE
    base = base[base.arm.notna()]
    base["lineage"] = base.CANCER_TYPE_DETAILED.fillna(base.CANCER_TYPE + " (unspecified)")
    base["strata"] = base.lineage.astype(str) + "|" + base.GENE_PANEL.astype(str)
    units = sorted(base.lineage.unique())
    Pby = {u: base[(base.lineage == u) & (base.arm == "Primary")].SAMPLE_ID.tolist() for u in units}
    Mby = {u: base[(base.lineage == u) & (base.arm == "Metastasis")].SAMPLE_ID.tolist() for u in units}
    n_ca = PM.per_cancer_k(Pby, Mby, minclass=MINCLASS)
    if not n_ca:
        raise SystemExit("no detailed type with > 30 primary-site and > 30 metastasis-site samples")
    pool = [s for u in n_ca for g in ("P", "M") for s in {"P": Pby, "M": Mby}[g][u]]
    posn = {s: i for i, s in enumerate(pool)}
    b = base.set_index("SAMPLE_ID").reindex(pool)
    WGD = b.wgd.astype(str).str.lower().eq("true").astype(int).to_numpy(float); STR = b.strata.to_numpy()
    print(f"{a.cancer_type}: matching units " + ", ".join(f"{u} ({len(Pby[u])} P / {len(Mby[u])} M -> k {k})" for u, k in n_ca.items())
          + f"; per run N = {sum(n_ca.values())} per arm; pool {len(pool)} (WGD+ {WGD.mean():.3f})", flush=True)
    seg = pd.read_csv(d / "segments.tsv.gz", sep="\t", dtype={"chrom": str, "cn_state": str})
    seg = seg[seg.SAMPLE_ID.isin(set(pool))]
    gene_bed = pd.read_csv(d / "gene_bed_hg19.tsv", sep="\t", dtype={"chrom": str})
    bundles, _, _, _ = W.build_seven_wgd(seg, pool, gene_bed, b.wgd.astype(str).str.lower().eq("true").astype(int).to_dict(),
                                         pd.to_numeric(b.ploidy, errors="coerce").to_dict())
    target = sorted(set(open(HERE.parent / "preprocess" / "target_genes_505.txt").read().split()))
    mm = pd.read_csv(d / "mutations.tsv.gz", sep="\t", usecols=["sample", "gene"])
    mm = mm[mm["sample"].isin(set(pool))]
    Xs, bs = snv_bundle(mm, pool, target)
    bs = bs.to_numpy(float)
    med = pd.Series(bs).groupby(STR).transform("median").to_numpy(float)
    hyper = bs > HYPER_FC * np.maximum(med, 1.0)
    STR_SNV = np.where(hyper, STR.astype(str) + "|HYPER", STR)
    X, BUR, UNITS = {}, {}, {}
    for m in CLASSES:
        X[m] = np.ascontiguousarray(np.asarray(bundles[m]["X"], float)); BUR[m] = np.asarray(bundles[m]["burden"], float); UNITS[m] = list(bundles[m]["units"])
    X["snv"] = np.ascontiguousarray(Xs.to_numpy(float)); BUR["snv"] = bs; UNITS["snv"] = list(Xs.columns)
    ids = PM.draw_idsets(Pby, Mby, n_ca, posn, runs=a.runs)
    IDSETS = {k: (v, strata(STR[v], len(v)), strata(STR_SNV[v], len(v))) for k, v in ids.items()}
    G.update(X=X, BUR=BUR, WGD=WGD, IDSETS=IDSETS, NSIM=a.nsim)
    tasks = [(gt, r, m, int(g)) for (gt, r), (v, _c, _cs) in IDSETS.items() for m in MODS for g in np.where(X[m][v].sum(0) > 0)[0]]
    print(f"{len(tasks)} (arm x run x modality x locus) tasks, null size {a.nsim}, {a.njob} processes ({time.time() - t:.0f}s)", flush=True)
    with mp.Pool(a.njob) as pl:
        res = pl.map(_worker, tasks, chunksize=16)
    byk = defaultdict(list)
    for rr in res:
        byk[rr[:3]].append(rr)
    rows = []
    for (gt, r, m), rr in byk.items():
        rr.sort(key=lambda z: z[3])
        freq = np.array([z[4] for z in rr]); O = np.array([z[5] for z in rr]); E = np.array([z[6] for z in rr])
        p2 = np.minimum(1.0, 2 * np.minimum([z[7] for z in rr], [z[8] for z in rr]))
        sg = significance(O, E, freq, np.vstack([z[9] for z in rr]), maxFDR=MAXFDR)
        nbin = bin_sizes(freq)
        for i, z in enumerate(rr):
            rows.append((gt, r, m, UNITS[m][z[3]], freq[i], O[i], E[i], float(O[i] - E[i]), float(sg["nES"][i]), float(p2[i]),
                         float(sg["nFDR2"][i]), int(nbin[i]), bool(nbin[i] > 1)))
    runs = collapse_focal_loh(pd.DataFrame(rows, columns=RUN_COLS), gene_bed)
    runs.to_csv(out / f"msk_primary_met_runs_{tag}.tsv", sep="\t", index=False)
    hi = max(1, round(PM.CAT_NSIG_HI * a.runs / PM.RUNS)); lo = round(PM.CAT_NSIG_LO * a.runs / PM.RUNS)
    obs, nul = PM.aggregate(runs)
    obs, thr = PM.categorize(obs, nul, hi=hi, lo=lo)
    obs.to_csv(out / f"msk_primary_met_{tag}.tsv", sep="\t", index=False)
    print(f"\n{len(obs)} loci; null threshold (99th percentile of |es_A - es_B|) {thr:.2f}; categories: {obs['cat'].value_counts().to_dict()}")
    print(obs[obs["cat"].ne("none")].sort_values("dist_obs", key=np.abs, ascending=False)[["modality", "gene", "es_P", "es_M", "p_count", "m_count", "dist_obs", "cat"]]
          .head(20).to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    pub = HERE / "expected" / "msk_primary_met_published.tsv"
    if pub.exists():
        from scipy.stats import spearmanr
        P = pd.read_csv(pub, sep="\t")
        j = obs.merge(P, on=["modality", "gene"], suffixes=("", "_pub"))
        rp = spearmanr(j.es_P, j.es_P_pub).statistic; rm = spearmanr(j.es_M, j.es_M_pub).statistic
        rd = spearmanr(j.es_M - j.es_P, j.es_M_pub - j.es_P_pub).statistic
        print(f"\ncheck against the published cohort-wide table: {len(j)} shared loci; Spearman rho es_P {rp:.2f}, es_M {rm:.2f}, "
              f"es_M - es_P {rd:.2f} (one cancer type against all 12; agreement in rank is expected, identity is not)")
    print(f"\ndone ({time.time() - t:.0f}s) -> {out}")

if __name__ == "__main__":
    main()
