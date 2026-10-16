"""Arm-annotated allele-specific segments (centromere-split arm assignment, sample and arm baselines), the
proximity-locus map of panel genes, and the per-class event bundles consumed by the selection engine.
"""
from __future__ import annotations
import numpy as np, pandas as pd

# hg38 chromosome length + centromere midpoint (bp), autosomes
_CHRLEN = {"1":248956422,"2":242193529,"3":198295559,"4":190214555,"5":181538259,"6":170805979,
 "7":159345973,"8":145138636,"9":138394717,"10":133797422,"11":135086622,"12":133275309,
 "13":114364328,"14":107043718,"15":101991189,"16":90338345,"17":83257441,"18":80373285,
 "19":58617616,"20":64444167,"21":46709983,"22":50818468}
_CEN = {"1":123e6,"2":93e6,"3":90e6,"4":50e6,"5":48e6,"6":60e6,"7":60e6,"8":45e6,"9":43e6,"10":39e6,
 "11":53e6,"12":35e6,"13":17e6,"14":17e6,"15":19e6,"16":37e6,"17":24e6,"18":18e6,"19":26e6,
 "20":28e6,"21":12e6,"22":15e6}
ARMLEN = {c+"p":_CEN[c] for c in _CEN}; ARMLEN.update({c+"q":_CHRLEN[c]-_CEN[c] for c in _CEN})

GENOME_CONSTANTS = {
    "hg38": (dict(_CHRLEN), dict(_CEN)),
    "hg19": ({"1":249250621,"2":243199373,"3":198022430,"4":191154276,"5":180915260,"6":171115067,
              "7":159138663,"8":146364022,"9":141213431,"10":135534747,"11":135006516,"12":133851895,"13":115169878,
              "14":107349540,"15":102531392,"16":90354753,"17":81195210,"18":78077248,"19":59128983,"20":63025520,
              "21":48129895,"22":51304566},
             {"1":125.0e6,"2":93.3e6,"3":91.0e6,"4":50.4e6,"5":48.4e6,"6":61.0e6,"7":59.9e6,"8":45.6e6,
              "9":49.0e6,"10":40.2e6,"11":53.7e6,"12":35.8e6,"13":17.9e6,"14":17.6e6,"15":19.0e6,"16":36.6e6,
              "17":24.0e6,"18":17.2e6,"19":26.5e6,"20":27.5e6,"21":13.2e6,"22":14.7e6}),
}
CURRENT_GENOME = "hg38"

def set_genome(build: str) -> str:
    """Switch the arm geometry (chromosome lengths, centromere midpoints, ARMLEN) to ``build`` (hg38 | hg19).

    The switch is process-global (the module dictionaries are rewritten in place): every later call of prep_segments,
    gene_loci or the callers in this process uses it. Call it once per cohort and never mix builds in one process
    without switching in between."""
    global CURRENT_GENOME
    L, C = GENOME_CONSTANTS[build]
    _CHRLEN.clear(); _CHRLEN.update(L); _CEN.clear(); _CEN.update(C)
    ARMLEN.clear(); ARMLEN.update({c + "p": _CEN[c] for c in _CEN}); ARMLEN.update({c + "q": _CHRLEN[c] - _CEN[c] for c in _CEN})
    CURRENT_GENOME = build
    return build

BROAD_FRAC = 0.80      # >=80% of covered arm territory -> broad
FOCAL_ARM_MAX = 0.50   # a segment is focal if it spans < 0.5 of its arm

ACROCENTRIC_P = {"13p", "14p", "15p", "21p", "22p"}

def prep_segments(seg, samples):
    """Normalise an allele-specific segment table: autosomes only, centromere-spanning segments cut at the centromere
    into a p and a q piece, acrocentric p arms (13p 14p 15p 21p 22p) dropped, arm assignment, integer CN, per-sample
    and per-arm length-weighted baselines. Returns the annotated frame."""
    s = seg.copy()
    if "c" not in s.columns:
        s["c"] = s["chrom"].astype(str).str.replace("chr", "", regex=False)
    s = s[s.SAMPLE_ID.isin(set(samples)) & s.c.isin(_CEN)].copy()
    for col in [c for c in ("seg_start", "seg_end", "seg_length", "tcn", "lcn") if c in s.columns]:
        s[col] = pd.to_numeric(s[col], errors="coerce")
    s = s.dropna(subset=["tcn", "seg_start", "seg_end"])
    s["tcn"] = s["tcn"].astype(int)
    cen = s.c.map(_CEN)
    span = (s.seg_start < cen) & (s.seg_end > cen)
    left = s[span].copy(); left["seg_end"] = cen[span]
    right = s[span].copy(); right["seg_start"] = cen[span]
    s = pd.concat([s[~span], left, right], ignore_index=True)
    s["seg_length"] = s.seg_end - s.seg_start
    s["mid"] = (s.seg_start + s.seg_end) / 2
    s["arm"] = s.c + np.where(s["mid"] < s.c.map(_CEN), "p", "q")
    s["armlen"] = s.arm.map(ARMLEN)
    s = s[~s.arm.isin(ACROCENTRIC_P)].copy()
    s["len"] = (s.seg_end - s.seg_start).clip(lower=0)
    # SAMPLE baseline = autosomal length-weighted modal tcn
    sb = (s.groupby(["SAMPLE_ID", "tcn"])["len"].sum().reset_index()
          .sort_values(["SAMPLE_ID", "len", "tcn"], ascending=[True, False, True])
          .drop_duplicates("SAMPLE_ID").set_index("SAMPLE_ID")["tcn"])
    s["Bs"] = s.SAMPLE_ID.map(sb)
    # ARM baseline = per-arm length-weighted modal tcn
    ab = (s.groupby(["SAMPLE_ID", "arm", "tcn"])["len"].sum().reset_index()
          .sort_values(["SAMPLE_ID", "arm", "len", "tcn"], ascending=[True, True, False, True])
          .drop_duplicates(["SAMPLE_ID", "arm"]).set_index(["SAMPLE_ID", "arm"])["tcn"])
    s["Ba"] = s.set_index(["SAMPLE_ID", "arm"]).index.map(ab)
    return s

def _bundle_broad(broad, samples, arm_col):
    """{X, burden, units, meta} for a broad class. unit = arm; burden = # broad arms of class."""
    idx = list(samples)
    hit = broad[broad[arm_col]]
    X = (hit.assign(v=1).pivot_table(index="SAMPLE_ID", columns="arm", values="v",
                                     aggfunc="max", fill_value=0).reindex(idx, fill_value=0))
    burden = hit.groupby("SAMPLE_ID").size().reindex(idx, fill_value=0)
    units = list(X.columns)
    meta = pd.DataFrame({"unit": units, "scale": "broad",
                         "freq": [float((X[u] > 0).mean()) for u in units]})
    return {"X": X.to_numpy(float), "burden": burden.to_numpy(float),
            "units": units, "samples": idx, "meta": meta}

def gene_loci(gene_bed, W=3e6, drivers=None):
    """Deterministic, event-independent locus map from panel-gene proximity: on each arm, genes sorted by position
    are grouped greedily so that a locus spans at most ``W`` bp from its first gene (bounded span, not single
    linkage); a gene-dense cluster collapses to one locus. Returns {gene -> locus label}; the label is a ``drivers``
    member if the locus holds one, else its alphabetically first gene."""
    gb = gene_bed.copy()
    gb["c"] = gb["chrom"].astype(str).str.replace("chr", "", regex=False)
    gb = gb[gb.c.isin(_CEN)].copy()
    gb["mid"] = pd.to_numeric(gb["mid"], errors="coerce")
    gb = gb.dropna(subset=["mid"])
    gb["arm"] = gb.c + np.where(gb["mid"] < gb.c.map(_CEN), "p", "q")
    drivers = set(drivers) if drivers else set()
    g2loc = {}
    def flush(mem):
        if not mem: return
        drv = [g for g in mem if g in drivers]
        lab = sorted(drv)[0] if drv else sorted(mem)[0]
        for g in mem: g2loc[g] = lab
    for arm, grp in gb.groupby("arm"):
        grp = grp.sort_values("mid"); members = []; anchor = None
        for _, r in grp.iterrows():
            if anchor is not None and r["mid"] - anchor > W:   # bounded span from the locus anchor
                flush(members); members = []; anchor = None
            if anchor is None: anchor = r["mid"]
            members.append(r["gene"])
        flush(members)
    return g2loc

def _bundle_focal(focal, gene_bed, samples, flag, g2loc, min_freq=0.005):
    """{X, burden, units, meta} for a focal class. X unit = proximity locus of panel genes (``g2loc`` from
    gene_loci); burden = # distinct loci altered (what X counts); tested units = loci altered in >= min_freq of
    the samples."""
    idx = list(samples)
    ev = focal[focal[flag]].copy()
    gb = gene_bed.rename(columns={"chrom": "c", "mid": "gmid"}).copy()
    gb["c"] = gb.c.astype(str).str.replace("chr", "", regex=False)
    gb["gmid"] = pd.to_numeric(gb["gmid"], errors="coerce")
    hits = ev.merge(gb[["gene", "c", "gmid"]], on="c")
    hits = hits[(hits.gmid >= hits.seg_start) & (hits.gmid <= hits.seg_end)].copy()
    if len(hits) == 0:
        return {"X": np.zeros((len(idx), 0)), "burden": np.zeros(len(idx)),
                "units": [], "samples": idx, "meta": pd.DataFrame(columns=["unit","scale","freq","members"])}
    hits["locus"] = hits.gene.map(lambda g: g2loc.get(g, g))
    full = (hits.assign(v=1).pivot_table(index="SAMPLE_ID", columns="locus", values="v",
                                         aggfunc="max", fill_value=0).reindex(idx, fill_value=0))
    burden = full.sum(1)
    X = full.loc[:, full.mean() >= min_freq]
    units = list(X.columns)
    memb = hits.groupby("locus")["gene"].apply(lambda s: ",".join(sorted(set(s)))).to_dict()
    meta = pd.DataFrame({"unit": units, "scale": "focal",
                         "freq": [float((X[u] > 0).mean()) for u in units],
                         "members": [memb.get(u, u)[:120] for u in units]})
    return {"X": X.to_numpy(float), "burden": burden.to_numpy(float),
            "units": units, "samples": idx, "meta": meta}

# canonical focal-amplification drivers (for locus labelling only; peak calls are CN-based)
AMP_DRIVERS = {"MYC","MYCN","MYCL","CCND1","CCND2","CCND3","CCNE1","CDK4","CDK6","MDM2","MDM4",
               "ERBB2","EGFR","FGFR1","FGFR2","FGFR3","FGFR4","KRAS","NRAS","MET","KIT","PDGFRA",
               "AR","BCL2","BCL6","MCL1","TERT","SOX2","GATA6","NKX2-1","YAP1","AKT1","AKT2","RICTOR"}
