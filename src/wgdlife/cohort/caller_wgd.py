"""WGD-aware CNA event caller: focal amplification and deletion by the doubling-invariant fold change relative to
the arm, with the expected-ploidy baseline (2 / 4 / 8 by WGD multiplicity); broad calls and focal LOH from the same
arm-annotated segments.
"""
from __future__ import annotations
import numpy as np
from wgdlife.cohort.cna_events6 import (
    BROAD_FRAC, FOCAL_ARM_MAX, prep_segments, gene_loci, _bundle_broad, _bundle_focal, AMP_DRIVERS,
)

AMP_FOLD_LOG2 = 1.0     # focal_amp: log2(TCN/Ba) >= +1  (TCN >= 2*Ba)
DEL_FOLD_LOG2 = -2.0    # focal_del relative arm: log2(TCN/Ba) <= -2  (TCN <= Ba/4); plus tcn==0 HOMDEL
AMP_FLOOR = 2           # focal_amp also needs tcn >= expected ploidy baseline + 2

def expected_ploidy_baseline(wgd, ploidy):
    """Expected ploidy baseline by WGD multiplicity: 2 without WGD, 4 after one doubling, 8 at ploidy >= 6."""
    return np.where(wgd == 0, 2.0, np.where(np.asarray(ploidy, float) < 6, 4.0, 8.0))

def wgd_broad_calls(s, Bcol="Bexp"):
    """Per (SAMPLE_ID, arm): broad_gain/broad_loss (tcn ≷ expected baseline Bcol) / broad_loh
    (minor==0). ≥80% covered arm territory. Returns one row per (sample, arm)."""
    s = s.copy()
    s["up"] = np.where(s.tcn > s[Bcol], s.len, 0.0)
    s["down"] = np.where(s.tcn < s[Bcol], s.len, 0.0)
    s["loh"] = np.where(s.lcn.fillna(-1) == 0, s.len, 0.0)
    g = s.groupby(["SAMPLE_ID", "arm"]).agg(covlen=("len", "sum"), up=("up", "sum"),
                                            dn=("down", "sum"), lh=("loh", "sum")).reset_index()
    g = g[g.covlen > 0]
    g["broad_gain"] = g.up / g.covlen >= BROAD_FRAC
    g["broad_loss"] = g.dn / g.covlen >= BROAD_FRAC
    g["broad_loh"] = g.lh / g.covlen >= BROAD_FRAC
    return g

def wgd_focal_events(s, broad):
    """Focal AMP/DEL by the doubling-invariant fold log2(TCN/Ba); focal_loh (minor==0) unchanged.
    Ba==0 → skip fold (arm homozygously deleted, handled by HOMDEL); tcn==0 → HOMDEL branch.
    focal_amp additionally requires tcn >= Bexp + AMP_FLOOR, where Bexp is the sample's expected ploidy baseline
    (2 / 4 / 8): a segment merely spared on a lost arm (Ba = 1, tcn = 2-4) is not an amplification; the floor makes
    the event a high-copy state relative to the genome's own ploidy (diploid >= 4, doubled >= 6)."""
    bmap = broad.set_index(["SAMPLE_ID", "arm"])["broad_loh"].to_dict()
    f = s[s.len < FOCAL_ARM_MAX * s.armlen].copy()
    Ba = f.Ba.to_numpy(float); tcn = f.tcn.to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        l2 = np.where(Ba > 0, np.log2(np.where(tcn > 0, tcn, np.nan) / Ba), np.nan)
    f["focal_amp"] = (l2 >= AMP_FOLD_LOG2) & (tcn >= f.Bexp.to_numpy(float) + float(AMP_FLOOR))
    f["focal_del"] = (tcn == 0) | (l2 <= DEL_FOLD_LOG2)
    arm_is_broad_loh = np.array([bmap.get((sid, a), False) for sid, a in zip(f.SAMPLE_ID, f.arm)])
    f["focal_loh"] = (f.lcn.fillna(-1) == 0) & (f.tcn >= 1) & (~arm_is_broad_loh)
    return f[f.focal_amp | f.focal_del | f.focal_loh].copy()

def build_seven_wgd(seg, samples, gene_bed, wgd_map, ploidy_map):
    """Event bundles ({X, burden, units, samples, meta} per class) from allele-specific segments: broad_gain,
    broad_loss and broad_loh at the sample baseline, focal_amp, focal_del and focal_loh at the arm baseline, for the
    samples in ``samples``; wgd_map / ploidy_map give each sample's WGD status and ploidy (the expected baseline).
    Returns (bundles, annotated segments, broad calls, focal events)."""
    s = prep_segments(seg, samples)
    s["w"] = s.SAMPLE_ID.map(wgd_map); s["pl"] = s.SAMPLE_ID.map(ploidy_map)
    s["Bexp"] = expected_ploidy_baseline(s.w.to_numpy(), s.pl.to_numpy())
    broad = wgd_broad_calls(s)
    focal = wgd_focal_events(s, broad)
    g2loc = gene_loci(gene_bed, drivers=AMP_DRIVERS)
    bundles = {
        "broad_gain": _bundle_broad(broad, samples, "broad_gain"),
        "broad_loss": _bundle_broad(broad, samples, "broad_loss"),
        "broad_loh":  _bundle_broad(broad, samples, "broad_loh"),
        "focal_amp":  _bundle_focal(focal, gene_bed, samples, "focal_amp", g2loc=g2loc),
        "focal_del":  _bundle_focal(focal, gene_bed, samples, "focal_del", g2loc=g2loc),
        "focal_loh":  _bundle_focal(focal, gene_bed, samples, "focal_loh", g2loc=g2loc),
    }
    return bundles, s, broad, focal
