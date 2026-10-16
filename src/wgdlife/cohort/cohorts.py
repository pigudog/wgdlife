"""The SNV event matrix and per-sample mutation burden consumed by the selection engine.
"""
from __future__ import annotations

import pandas as pd

def snv_bundle(mut_long, samples, genes):
    """SNV: X = 0/1 over the tested ``genes``; burden = each sample's count of non-synonymous calls within the tested
    genes (the opportunity). ``mut_long``: [sample, gene], already non-synonymous-filtered."""
    idx = list(samples)
    sub = mut_long[mut_long["gene"].isin(set(genes))]
    burden = sub.groupby("sample").size().reindex(idx, fill_value=0)
    X = (sub.assign(v=1).pivot_table(index="sample", columns="gene", values="v",
                                     aggfunc="max", fill_value=0).reindex(idx, fill_value=0)
         if not sub.empty else pd.DataFrame(0.0, index=idx, columns=[]))
    return X, burden
