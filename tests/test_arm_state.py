"""Arm-level allele-state timing rule: 2+0 -> before WGD, 2+1 -> after WGD, >= 80% of the arm."""
import pandas as pd

from wgdlife.timing import arm_state


def _seg(sample, chrom, start, end, tcn, lcn):
    return dict(SAMPLE_ID=sample, chrom=chrom, seg_start=start, seg_end=end, seg_length=end - start, tcn=tcn, lcn=lcn)


def test_rule_on_whole_arm_states():
    seg = pd.DataFrame([_seg("a", "17", 1, 22_000_000, 2, 0),        # 17p 2+0 -> pre
                        _seg("a", "17", 26_000_000, 83_000_000, 3, 1),  # 17q 2+1 -> post
                        _seg("b", "17", 1, 22_000_000, 4, 2),          # 17p 2+2 -> amb
                        _seg("b", "17", 26_000_000, 83_000_000, 3, 1)])
    s = arm_state.prep_segments(seg, ["a", "b"])
    cfg = arm_state.arm_configs(s).set_index(["SAMPLE_ID", "arm"])
    assert cfg.loc[("a", "17p"), "timing"] == "pre" and cfg.loc[("a", "17q"), "timing"] == "post"
    assert cfg.loc[("b", "17p"), "timing"] == "amb"
    t = arm_state.arm_timing(cfg.reset_index(), "toy").set_index("arm")
    assert (t.loc["17p", "n_pre"], t.loc["17p", "n_post"]) == (1, 0) and (t.loc["17q", "n_pre"], t.loc["17q", "n_post"]) == (0, 2)
