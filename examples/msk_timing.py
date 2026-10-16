"""Example: timing of mutations and arm losses relative to WGD in MSK-IMPACT, from the public release.
"""
import argparse
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from wgdlife.cohort.constants import CANCER_TYPES_12
from wgdlife.timing.arm_state import arm_configs, prep_segments, wilson
from wgdlife.timing.molecular_timing import panel_classify

HERE = Path(__file__).resolve().parent
ONCOGENIC = {"Oncogenic", "Likely Oncogenic"}
EPOCH = {"pre_wgd": "pre", "post_wgd": "post", "subclonal": "subclonal"}

def master_table(d: Path) -> pd.DataFrame:
    """One row per timed event (sample, event, class, epoch) in the WGD+ carcinomas."""
    t = time.time()
    S = pd.read_csv(d / "samples.tsv", sep="\t", low_memory=False)
    cs = S[S.facets_qc.astype(str).str.lower().eq("true") & S.CANCER_TYPE.isin(CANCER_TYPES_12)]
    wpos = cs[cs.wgd.astype(str).str.lower().eq("true")]
    sset = set(wpos.SAMPLE_ID)
    print(f"WGD+ carcinomas with a QC-passed call: {len(sset):,} samples, {wpos.PATIENT_ID.nunique():,} patients ({time.time() - t:.0f}s)", flush=True)
    target = set(open(HERE.parent / "preprocess" / "target_genes_505.txt").read().split())
    mm = pd.read_csv(d / "mutations.tsv.gz", sep="\t", low_memory=False)
    mm = mm[mm["sample"].isin(sset) & mm.gene.isin(target)]
    sc = panel_classify(mm)
    sc = sc[sc.timing.isin(EPOCH)]
    snv = pd.DataFrame({"sample": sc["sample"].values, "event": sc.gene.values, "class": "SNV", "epoch": sc.timing.map(EPOCH).values,
                        "is_oncogenic": mm.loc[sc.index, "ONCOGENIC"].isin(ONCOGENIC).values})
    print(f"timed mutations: {len(snv):,} ({time.time() - t:.0f}s)", flush=True)
    seg = pd.read_csv(d / "segments.tsv.gz", sep="\t", dtype={"chrom": str, "cn_state": str})
    s = prep_segments(seg[seg.SAMPLE_ID.isin(sset)], sorted(sset))
    bl = arm_configs(s)
    bl = bl[bl.timing.isin(["pre", "post"])]
    loss = pd.DataFrame({"sample": bl.SAMPLE_ID.values, "event": bl.arm.values, "class": "broad_loss",
                         "epoch": bl.timing.values, "is_oncogenic": True})
    print(f"timed arm losses: {len(loss):,} ({time.time() - t:.0f}s)", flush=True)
    M = pd.concat([snv, loss], ignore_index=True)
    info = cs.set_index("SAMPLE_ID")
    M["patient"] = M["sample"].map(info.PATIENT_ID)
    M["cancer"] = M["sample"].map(info.CANCER_TYPE)
    site = M["sample"].map(info.METASTATIC_SITE).fillna("Not Applicable").astype(str)   # specimen site, not the patient's disease status
    M["status"] = np.where(site.isin(["", "Not Applicable"]), "Primary", "Metastasis")
    return M

def cell(d: pd.DataFrame):
    """patients timed pre or post, P(pre-WGD) and its Wilson 95% interval."""
    k = int((d.epoch == "pre").sum()); n = int(d.epoch.isin(["pre", "post"]).sum())
    if n == 0:
        return n, np.nan, np.nan, np.nan
    p = k / n; lo, hi = wilson(k, n)
    return n, p, lo, hi

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("WGDLIFE_DATA", "data"))
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    d = Path(a.data) / "msk_impact_50k" / "derived"
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    M = master_table(d)
    M.to_csv(out / "msk_timing_events.tsv", sep="\t", index=False)
    # one record per patient x event; mutations restricted to oncogenic calls (the baseline row keeps all of them)
    M = M.drop_duplicates(["patient", "class", "event"])
    M_all = M
    M = M[M["class"].ne("SNV") | M.is_oncogenic]
    exp = pd.read_csv(HERE / "expected" / "msk_timing.tsv", sep="\t")
    rows = []
    for r in exp.itertuples(index=False):
        if r.event.startswith("All non-synonymous"):
            sub = M_all[M_all["class"].eq("SNV")]
        elif r.event.startswith("All timed arm losses"):
            sub = M[M["class"].eq("broad_loss")]
        elif r.event_class == "Arm loss":
            sub = M[M["class"].eq("broad_loss") & M.event.eq(r.event.replace(" loss", ""))]
        else:
            sub = M[M["class"].eq("SNV") & M.event.eq(r.event)]
        n, p, lo, hi = cell(sub); _, pp, _, _ = cell(sub[sub.status.eq("Primary")]); _, pm, _, _ = cell(sub[sub.status.eq("Metastasis")])
        rows.append(dict(event=r.event, event_class=r.event_class, n=n, P_pre=p, ci_lo=lo, ci_hi=hi, P_pre_primary=pp, P_pre_metastasis=pm,
                         n_pub=r.n, P_pre_pub=r.P_pre, P_pre_primary_pub=r.P_pre_primary, P_pre_metastasis_pub=r.P_pre_metastasis))
    T = pd.DataFrame(rows)
    T.to_csv(out / "msk_timing.tsv", sep="\t", index=False)
    print("\nevent timing relative to WGD (this run | published):")
    print(T.assign(P=T.P_pre.map(lambda v: f"{v:.2f}"), CI=T.apply(lambda r: f"({r.ci_lo:.2f}-{r.ci_hi:.2f})", axis=1))
          [["event", "event_class", "n", "n_pub", "P", "CI", "P_pre_pub", "P_pre_primary", "P_pre_metastasis"]]
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    same_n = int((T.n == T.n_pub).sum()); same_p = int((np.round(T.P_pre, 2) == T.P_pre_pub).sum())
    print(f"\ncheck: patient counts identical for {same_n}/{len(T)} rows; P(pre-WGD) identical at two decimals for {same_p}/{len(T)} rows -> {out / 'msk_timing.tsv'}")

if __name__ == "__main__":
    main()
