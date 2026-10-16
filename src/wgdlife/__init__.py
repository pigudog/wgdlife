"""wgdlife — the life history of whole-genome doubling across human cancers.

The computational methods of Zhang et al. (2026): the WGD criterion and CNA event decomposition (`cohort`), the
WGD-anchored selection framework (`selection`), event timing relative to WGD with its consensus
meta-analysis (`timing`; the convergence test of the event spectrum is the R engine in `R/`), and the
matched primary–metastasis comparison (`metastasis`). Engines only: the study's data-loading pipelines and the
figure code live elsewhere (see README).
"""
__version__ = "1.0.0"
from . import cohort, selection, timing, metastasis

__all__ = ["cohort", "selection", "timing", "metastasis"]
