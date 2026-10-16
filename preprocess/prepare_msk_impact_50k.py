"""Turn the three MSK-IMPACT 50K release files into the compact tables the examples read.
"""
import argparse
import os
import time
from pathlib import Path

import pandas as pd

from wgdlife.cohort.constants import SILENT_CLASSES, current_symbol

CLIN_COLS = ["SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE", "DISEASE_STATUS",
             "METASTATIC_SITE", "GENE_PANEL", "MSI_TYPE", "TMB", "facets_qc", "purity", "wgd", "ploidy"]
MAF_COLS = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "tcn", "lcn", "expected_alt_copies",
            "clonality", "ONCOGENIC"]
HERE = Path(__file__).resolve().parent

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("WGDLIFE_DATA", "data"))
    a = ap.parse_args()
    src = Path(a.data) / "msk_impact_50k"
    out = src / "derived"
    out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    # clinical
    cs = pd.read_csv(src / "data_clinical_sample.txt", sep="\t", comment="#", low_memory=False, usecols=CLIN_COLS)
    cs.to_csv(out / "samples.tsv", sep="\t", index=False)
    print(f"samples.tsv: {len(cs):,} samples ({time.time() - t:.0f}s)", flush=True)
    # segments + gene coordinates from the gene-level FACETS table (streamed; 3.5 GB)
    cols = ["SAMPLE_ID", "gene", "chrom", "gene_start", "gene_end", "seg_start", "seg_end", "seg_length", "tcn.em", "lcn.em", "cn_state"]
    parts, genes = [], []
    for i, ch in enumerate(pd.read_csv(src / "data_facets_gene_level.txt", sep="\t", usecols=cols, dtype={"chrom": str, "cn_state": str},
                                       low_memory=False, chunksize=2_000_000)):
        genes.append(ch[["gene", "chrom", "gene_start", "gene_end"]].drop_duplicates())
        parts.append(ch.drop(columns=["gene", "gene_start", "gene_end"]).drop_duplicates(["SAMPLE_ID", "chrom", "seg_start", "seg_end"]))
        if i % 10 == 0:
            print(f"  gene-level chunk {i} ({time.time() - t:.0f}s)", flush=True)
    seg = pd.concat(parts, ignore_index=True).drop_duplicates(["SAMPLE_ID", "chrom", "seg_start", "seg_end"])
    seg["chrom"] = seg.chrom.astype(str).str.replace("chr", "", regex=False)
    for c in ["seg_start", "seg_end", "seg_length", "tcn.em", "lcn.em"]:
        seg[c] = pd.to_numeric(seg[c], errors="coerce")
    seg = seg.rename(columns={"tcn.em": "tcn", "lcn.em": "lcn"})[["SAMPLE_ID", "chrom", "seg_start", "seg_end", "seg_length", "tcn", "lcn", "cn_state"]]
    seg.to_csv(out / "segments.tsv.gz", sep="\t", index=False)
    print(f"segments.tsv.gz: {len(seg):,} segments in {seg.SAMPLE_ID.nunique():,} samples ({time.time() - t:.0f}s)", flush=True)
    g = pd.concat(genes).drop_duplicates()
    g["gene"] = current_symbol(g.gene)
    g["chrom"] = g.chrom.astype(str).str.replace("chr", "", regex=False)
    g["mid"] = (g.gene_start + g.gene_end) / 2
    target = set(open(HERE / "target_genes_505.txt").read().split())
    bed = g[g.gene.isin(target)].sort_values(["chrom", "mid"])[["gene", "chrom", "mid"]]
    bed.to_csv(out / "gene_bed_hg19.tsv", sep="\t", index=False)
    print(f"gene_bed_hg19.tsv: {len(bed)} panel genes with release coordinates ({len(target - set(bed.gene))} target genes absent from the release)", flush=True)
    # mutations (non-silent)
    parts = []
    for ch in pd.read_csv(src / "data_mutations_ccf.maf.txt", sep="\t", comment="#", low_memory=False, usecols=MAF_COLS, chunksize=1_000_000):
        parts.append(ch[~ch.Variant_Classification.isin(SILENT_CLASSES)])
    mm = pd.concat(parts, ignore_index=True).rename(columns={"Tumor_Sample_Barcode": "sample", "Hugo_Symbol": "gene"})
    mm["gene"] = current_symbol(mm.gene)   # release symbols -> current HGNC, as for the gene bed
    mm = mm[["sample", "gene"] + [c for c in MAF_COLS if c not in ("Tumor_Sample_Barcode", "Hugo_Symbol")]]   # documented column order
    mm.to_csv(out / "mutations.tsv.gz", sep="\t", index=False)
    print(f"mutations.tsv.gz: {len(mm):,} non-silent calls in {mm['sample'].nunique():,} samples ({time.time() - t:.0f}s)")

if __name__ == "__main__":
    main()
