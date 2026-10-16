"""WGD-aware caller contract: thresholds and the expected ploidy baseline.
"""
import numpy as np

from wgdlife.cohort import caller_wgd as W

def test_thresholds():
    assert W.AMP_FOLD_LOG2 == 1.0 and W.DEL_FOLD_LOG2 == -2.0 and W.AMP_FLOOR == 2

def test_expected_ploidy_baseline():
    wgd = np.array([0, 1, 1, 1]); ploidy = np.array([2.0, 3.2, 5.9, 6.2])
    assert list(W.expected_ploidy_baseline(wgd, ploidy)) == [2.0, 4.0, 4.0, 8.0]
