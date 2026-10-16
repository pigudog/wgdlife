"""Arm assignment in cna_events6.prep_segments: centromere-spanning segments are split, acrocentric p arms dropped.
"""
import pandas as pd
import pytest

from wgdlife.cohort import cna_events6 as c6

@pytest.fixture(autouse=True)
def _hg38_geometry():
    """These synthetic segments are written in hg38 coordinates; pin the geometry so the test does not depend on
    whichever platform last called set_genome() in the same process."""
    c6.set_genome("hg38")
    yield

def _seg():
    return pd.DataFrame({"SAMPLE_ID": ["s1"] * 3, "chrom": ["9", "13", "1"],
                         "seg_start": [1, 1, 1], "seg_end": [138_000_000, 114_000_000, 248_000_000],
                         "seg_length": [138_000_000, 114_000_000, 248_000_000], "tcn": [1, 2, 2], "lcn": [0, 1, 1]})

def test_whole_chromosome_lands_on_both_arms_and_acrocentric_p_is_dropped():
    s = c6.prep_segments(_seg(), ["s1"])
    assert {"9p", "9q", "1p", "1q", "13q"} <= set(s.arm) and "13p" not in set(s.arm)
    assert (s.seg_end - s.seg_start).sum() < _seg().seg_length.sum()   # acrocentric p piece removed
