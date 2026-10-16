"""Example: WGD-anchored selection in one MSK-IMPACT cancer type, from the public release.
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
from wgdlife.selection.anchored import anchored_locus, significance, stable_seed
from wgdlife.selection.inference import strata

HERE = Path(__file__).resolve().parent
NE_ONCOTREE = {"SCBC", "UNEC", "GINETES", "STSC", "HNNE", "ALUCA", "LUNE", "LUCA", "LNET", "SCCO", "PANET", "PRNE", "PRSCC", "SCLC", "NEC", "NET"}
DETAILED_FLOOR, MINN, HYPER_FC = 50, 80, 10.0
NSIM_CNA, NSIM_SNV, MAXFDR, MIN_FREQ = 10000, 2000, 0.10, 0.005
CLASSES = ["broad_gain", "broad_loh", "focal_amp", "focal_del", "focal_loh"]
G = {}

def discovery_primaries(samples: pd.DataFrame) -> pd.DataFrame:
    """Primary tumours of the 12 carcinoma types with a FACETS-QC-passed WGD call; neuroendocrine histologies dropped;
    lineage stratum = detailed type when it holds >= 50 such tumours, else '<type> (other)', x gene panel, with MSI-high
    tumours as their own stratum."""
    m = samples[samples.CANCER_TYPE.isin(CANCER_TYPES_12) & samples.facets_qc.astype(str).str.lower().eq("true")
                & samples.wgd.notna() & samples.DISEASE_STATUS.eq("Primary") & ~samples.ONCOTREE_CODE.isin(NE_ONCOTREE)].copy()
    m["w"] = m.wgd.astype(str).str.lower().eq("true").astype(int)
    dn = m.groupby(["CANCER_TYPE", "CANCER_TYPE_DETAILED"]).SAMPLE_ID.transform("size")
    m["lineage"] = np.where(dn >= DETAILED_FLOOR, m.CANCER_TYPE_DETAILED, m.CANCER_TYPE + " (other)")
    m["stratum"] = m.lineage.astype(str) + "|" + m.GENE_PANEL.astype(str)
    m.loc[m.MSI_TYPE.eq("Instable"), "stratum"] += "|MSI"
    return m

def _cna_worker(task):
    ck, cls, g = task
    ids = G["ids"][ck]
    x = G["X"][cls][ids][:, g]
    O, E, sd, pe, pdl, resid = anchored_locus(x, G["wgd"][ids], G["cidx"][ck], G["bur"][cls][ids],
                                              stable_seed(ck, cls) * 1000 + g, G["nsim"])
    return ck, cls, g, float(x.mean()), O, E, resid, float(min(1.0, 2 * min(pe, pdl)))

def _snv_worker(g):
    x = G["Xs"][:, g]
    O, E, sd, pe, pdl, resid = anchored_locus(x, G["wgd_s"], G["cidx_s"], G["bur_s"], 777000 + g, G["nsim_s"])
    return g, float(x.mean()), O, E, resid, float(min(1.0, 2 * min(pe, pdl)))

def _table(ck, cls, rows):
    rows.sort(key=lambda r: r[0])
    units = [r[1] for r in rows]
    freq = np.array([r[2] for r in rows]); O = np.array([r[3] for r in rows]); E = np.array([r[4] for r in rows])
    resid = np.vstack([r[5] for r in rows]); p2 = np.array([r[6] for r in rows])
    sg = significance(O, E, freq, resid, maxFDR=MAXFDR)
    ER = O / np.clip(E, 1e-9, None)
    return pd.DataFrame({"cohort": ck, "class": cls, "unit": units, "frequency": freq, "O": O, "E": E, "p_two": p2, "ER": ER,
                         "wES": sg["wES"], "nES": sg["nES"], "nFDR2": sg["nFDR2"], "sig_type": sg["type"],
                         "sig10": sg["nFDR2"] <= MAXFDR, "log2ER": np.log2(np.clip(ER, 1e-9, None))})

def compare(res: pd.DataFrame, expected: Path) -> None:
    if not expected.exists():
        print(f"(no expected table at {expected})")
        return
    exp = pd.read_csv(expected, sep="\t")
    m = res.merge(exp, on=["cohort", "class", "unit"], suffixes=("", "_pub"))
    miss = len(exp) - len(m)
    print(f"check against the published values ({expected.name}): {len(m)} loci matched, {miss} published loci not reproduced")
    for c in ("O", "E", "p_two", "nES", "nFDR2"):
        d = (m[c] - m[f"{c}_pub"]).abs()
        print(f"  max |delta {c}| = {d.max():.3g}  (identical: {int((d < 1e-9).sum())}/{len(m)})")
    same = (m.sig10.astype(bool) == m.sig10_pub.astype(bool)).mean()
    print(f"  significance calls agree for {100 * same:.1f}% of loci")

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cancer-type", default="Breast Cancer", help="one of the 12 MSK cancer types (default Breast Cancer); 'none' = pan-cancer only")
    ap.add_argument("--pancancer", action="store_true", help="also run the pan-cancer cohort (all 12 types)")
    ap.add_argument("--njob", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--nsim-cna", type=int, default=NSIM_CNA); ap.add_argument("--nsim-snv", type=int, default=NSIM_SNV)
    ap.add_argument("--data", default=os.environ.get("WGDLIFE_DATA", "data"))
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    c6.set_genome("hg19")   # MSK-IMPACT FACETS segments are GRCh37
    d = Path(a.data) / "msk_impact_50k" / "derived"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    S = pd.read_csv(d / "samples.tsv", sep="\t", low_memory=False)
    m = discovery_primaries(S)
    seg = pd.read_csv(d / "segments.tsv.gz", sep="\t", dtype={"chrom": str, "cn_state": str})
    m = m[m.SAMPLE_ID.isin(set(seg.SAMPLE_ID))]
    samples = sorted(m.SAMPLE_ID)
    m = m.set_index("SAMPLE_ID").reindex(samples)
    WGD = m.w.to_numpy(float); STR = m.stratum.to_numpy(); CT = m.CANCER_TYPE.to_numpy()
    print(f"discovery primaries: {len(samples):,} (WGD+ {WGD.mean():.3f}; {m.stratum.nunique()} strata) ({time.time() - t:.0f}s)", flush=True)
    gene_bed = pd.read_csv(d / "gene_bed_hg19.tsv", sep="\t", dtype={"chrom": str})
    bundles, _, _, _ = W.build_seven_wgd(seg[seg.SAMPLE_ID.isin(set(samples))], samples, gene_bed, m.w.to_dict(),
                                         pd.to_numeric(m.ploidy, errors="coerce").to_dict())
    cohorts = {}
    if a.pancancer or a.cancer_type == "none":
        cohorts["pancancer"] = np.arange(len(samples))
    if a.cancer_type != "none":
        if a.cancer_type not in CANCER_TYPES_12:
            raise SystemExit(f"--cancer-type must be one of {list(CANCER_TYPES_12)}")
        ids = np.where(CT == a.cancer_type)[0]
        if not (len(ids) >= MINN and WGD[ids].sum() >= 10 and (1 - WGD[ids]).sum() >= 10):
            raise SystemExit(f"{a.cancer_type}: fewer than {MINN} tumours or fewer than 10 per WGD class")
        cohorts[a.cancer_type] = ids
    G.update(ids=cohorts, cidx={k: strata(STR[v], len(v)) for k, v in cohorts.items()},
             X={c: np.ascontiguousarray(bundles[c]["X"]) for c in CLASSES},
             bur={c: np.asarray(bundles[c]["burden"], float) for c in CLASSES},
             wgd=WGD, nsim=a.nsim_cna)
    units = {c: list(bundles[c]["units"]) for c in CLASSES}
    tasks = [(ck, cls, int(g)) for ck, ids in cohorts.items() for cls in CLASSES for g in np.where(G["X"][cls][ids].sum(0) > 0)[0]]
    print(f"copy-number loci: {len(tasks)} (cohort x class x locus) tasks, null size {a.nsim_cna}, {a.njob} processes ({time.time() - t:.0f}s)", flush=True)
    with mp.Pool(a.njob) as pool:
        res = pool.map(_cna_worker, tasks, chunksize=8)
    rows = defaultdict(list)
    for ck, cls, g, freq, O, E, resid, p2 in res:
        rows[(ck, cls)].append((g, units[cls][g], freq, O, E, resid, p2))
    tables = [_table(ck, cls, rr) for (ck, cls), rr in rows.items()]
    print(f"copy-number classes done ({time.time() - t:.0f}s); mutations...", flush=True)
    # non-silent mutations in the panel genes: opportunity = each tumour's count over the tested genes;
    # ultramutators (count > 10 x their stratum median) form their own stratum
    target = sorted(set(open(HERE.parent / "preprocess" / "target_genes_505.txt").read().split()))
    mm = pd.read_csv(d / "mutations.tsv.gz", sep="\t", usecols=["sample", "gene"])
    for ck, ids in cohorts.items():
        sc = [samples[i] for i in ids]
        mc = mm[mm["sample"].isin(set(sc))]
        Xdf, bur = snv_bundle(mc, sc, target)
        bur = bur.to_numpy(float)
        strata_snv = STR[ids].astype(object)
        med = pd.Series(bur).groupby(strata_snv).transform("median").to_numpy(float)
        hyper = bur > HYPER_FC * np.maximum(med, 1.0)
        strata_snv = np.where(hyper, strata_snv.astype(str) + "|HYPER", strata_snv)
        X = Xdf.to_numpy(float); keep = [j for j in range(X.shape[1]) if X[:, j].mean() >= MIN_FREQ]
        X = np.ascontiguousarray(X[:, keep]); su = [Xdf.columns[j] for j in keep]
        G.update(Xs=X, wgd_s=WGD[ids], cidx_s=strata(strata_snv, len(sc)), bur_s=bur, nsim_s=a.nsim_snv)
        with mp.Pool(a.njob) as pool:
            rs = pool.map(_snv_worker, range(len(su)), chunksize=8)
        tables.append(_table(ck, "snv", [(g, su[g], freq, O, E, resid, p2) for g, freq, O, E, resid, p2 in rs]))
        print(f"  {ck}: {len(su)} genes tested, {int(hyper.sum())} ultramutators ({time.time() - t:.0f}s)", flush=True)
    res = pd.concat(tables, ignore_index=True)
    for ck in cohorts:
        f = out / f"msk_selection_{ck.replace(' ', '_')}.tsv"
        r = res[res.cohort.eq(ck)]
        r.to_csv(f, sep="\t", index=False)
        sig = r[r.sig10]
        print(f"\n{ck}: {len(r)} loci tested, {len(sig)} at FDR <= {MAXFDR} ({int((sig.sig_type == 'CO').sum())} co-occurring with WGD, "
              f"{int((sig.sig_type == 'ME').sum())} mutually exclusive) -> {f}")
        print(sig.sort_values("nFDR2")[["class", "unit", "frequency", "O", "E", "ER", "p_two", "nFDR2", "sig_type"]].head(15)
              .to_string(index=False, float_format=lambda v: f"{v:.3g}"))
        compare(r, HERE / "expected" / f"msk_selection_{ck.replace(' ', '_')}.tsv")
    print(f"\ndone ({time.time() - t:.0f}s)")

if __name__ == "__main__":
    main()
