"""Mutation timing on allele-specific copy-number states (timing.molecular_timing.panel_classify).
"""
import pandas as pd

from wgdlife.timing.molecular_timing import panel_classify

def test_panel_classify_states():
    d = pd.DataFrame({"tcn": [4, 4, 3, 3, 2, 4, 1], "lcn": [2, 2, 1, 1, 0, 2, 0],
                      "expected_alt_copies": [2, 1, 2, 1, 2, 2, 1],
                      "clonality": ["CLONAL"] * 5 + ["SUBCLONAL", "CLONAL"]})
    t = panel_classify(d)
    assert list(t.timing) == ["pre_wgd", "post_wgd", "pre_wgd_unbal", "ambiguous", "pre_wgd", "subclonal", "uninformative"]
    assert list(t.timing_state[:2]) == ["2+2", "2+2"] and t.timing_state[4] == "2+0" and t.timing_state[2] is None
