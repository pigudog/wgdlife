"""The WGD criterion: length-weighted major-CN >= 2 fraction of the autosomal genome.
"""
from pathlib import Path

import pandas as pd

from wgdlife.cohort.wgd_label import call_wgd_from_fraction, mcn_ge2_fraction_from_segments

EX = Path(__file__).resolve().parents[1] / "examples" / "tcga_wgd_example"

def test_segment_fraction_is_length_weighted_and_autosomal():
    seg = pd.DataFrame({"Chromosome": ["chr1", "chr1", "chrX"], "Start": [0, 100, 0], "End": [100, 400, 1000],
                        "Major_Copy_Number": [2, 1, 2]})
    f = mcn_ge2_fraction_from_segments(seg)
    assert abs(f - 0.25) < 1e-12                                          # 100 of 400 autosomal bp; chrX ignored
    assert call_wgd_from_fraction(f) == 0 and call_wgd_from_fraction(0.51) == 1

def test_public_example_reproduces_published_calls():
    exp = pd.read_csv(EX / "expected.tsv", sep="\t")
    for _, r in exp.iterrows():
        seg = pd.read_csv(EX / "segments" / f"{r['sample']}.ascat3.seg.txt", sep="\t")
        f = mcn_ge2_fraction_from_segments(seg)
        assert abs(f - r.mcn_ge2_frac) < 1e-6 and call_wgd_from_fraction(f) == r.wgd, r["sample"]
