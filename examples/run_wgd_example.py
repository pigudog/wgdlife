"""Minimal public example: WGD criterion + arm-loss timing on 24 TCGA ASCAT3 profiles (see README.md)."""
from pathlib import Path

import pandas as pd

from wgdlife.cohort.wgd_label import call_wgd_from_fraction, mcn_ge2_fraction_from_segments
from wgdlife.timing import arm_state

HERE = Path(__file__).resolve().parent / "tcga_wgd_example"


def main():
    expected = pd.read_csv(HERE / "expected.tsv", sep="\t")
    rows, segs = [], []
    for _, r in expected.iterrows():
        seg = pd.read_csv(HERE / "segments" / f"{r['sample']}.ascat3.seg.txt", sep="\t")
        frac = mcn_ge2_fraction_from_segments(seg)
        rows.append(dict(sample=r["sample"], cancer=r["cancer"], mcn_ge2_frac=frac, wgd=call_wgd_from_fraction(frac),
                         expected_frac=r["mcn_ge2_frac"], expected_wgd=r["wgd"]))
        s = seg.rename(columns={"Chromosome": "chrom", "Start": "seg_start", "End": "seg_end", "Copy_Number": "tcn",
                                "Minor_Copy_Number": "lcn"})
        s["chrom"] = s.chrom.astype(str).str.replace("chr", "", regex=False)
        s["seg_length"] = (s.seg_end - s.seg_start).clip(lower=0)
        s["SAMPLE_ID"] = r["sample"]
        segs.append(s[["SAMPLE_ID", "chrom", "seg_start", "seg_end", "seg_length", "tcn", "lcn"]])
    out = pd.DataFrame(rows)
    out["abs_diff"] = (out.mcn_ge2_frac - out.expected_frac).abs()
    print(out.to_string(index=False, float_format=lambda x: f"{x:.6f}"))
    assert (out.wgd == out.expected_wgd).all(), "WGD calls differ from the published calls"
    assert out.abs_diff.max() < 1e-6, "fractions differ from the published values"
    print(f"\nWGD calls reproduced for {len(out)} tumours (max |diff| in fraction = {out.abs_diff.max():.1e})")

    # arm-loss timing by allele state in the WGD+ tumours: 2+0 = before WGD, 2+1 = after WGD
    seg = pd.concat(segs, ignore_index=True)
    wgd_pos = out.loc[out.wgd == 1, "sample"]
    s = arm_state.prep_segments(seg, wgd_pos, build="hg38")   # GDC ASCAT3 segments are GRCh38
    cfg = arm_state.arm_configs(s)
    timing = arm_state.arm_timing(cfg, "TCGA-example")
    print("\nArm-loss timing in the 12 WGD+ tumours (n_pre = 2+0, n_post = 2+1):")
    print(timing.to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()
