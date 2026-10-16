# wgdlife

Computational methods of *The life history of whole-genome doubling across human cancers* (Zhang et al., 2026):
WGD classification, WGD-anchored selection analysis, event timing with its consensus meta-analysis, and the
matched primary–metastasis comparison. Python 3.10; the convergence test (`R/`) is R.

## Environment

```bash
conda create -n wgdlife python=3.10.20 -y
conda activate wgdlife
pip install -r requirements.lock.txt
pip install -e .
pip install pytest==9.1.1            # only for `make test`
Rscript R/requirements.R        # only for R/convergence_test.R (R 4.4.3, data.table 1.18.6.1)
```

## Configuration

Everything that reads data takes its directory from `WGDLIFE_DATA` (default: `data/` relative to the current working
directory, so run the commands below from the repository root):
`preprocess/` writes the MSK-IMPACT tables to `data/msk_impact_50k/`, and the R convergence test reads
`fig3_convergence_events_long.tsv` and `fig3_convergence_event_universe.tsv` (built by the companion repository
`wgdlife_analysis`; see `R/README.md`) from `data/`.

## Genome build

Arm geometry defaults to GRCh38 (TCGA ASCAT3). For GRCh37 data (the MSK-IMPACT release, ICGC/PCAWG, TRACERx) call
`wgdlife.cohort.cna_events6.set_genome("hg19")` before `prep_segments` or `build_seven_wgd`; the MSK-IMPACT examples do this, and
`wgdlife.timing.arm_state.prep_segments(seg, samples, build=...)` switches the geometry itself (hg19 by default).

## Minimal example

```bash
make smoke      # python examples/run_wgd_example.py
make test       # pytest
```

`examples/tcga_wgd_example/` holds 24 open-access TCGA allele-specific copy-number profiles (GDC ASCAT3). The
example computes the length-weighted fraction of the autosomal genome with major copy number ≥ 2, calls WGD above 0.5,
checks the calls against the published values and times the recurrent arm losses by allele state.

```python
import pandas as pd
from wgdlife.cohort import mcn_ge2_fraction_from_segments, call_wgd_from_fraction
seg = pd.read_csv("examples/tcga_wgd_example/segments/TCGA-04-1331-01A.ascat3.seg.txt", sep="\t")
call_wgd_from_fraction(mcn_ge2_fraction_from_segments(seg))      # 1
```

Modules: `wgdlife.cohort` (WGD criterion, copy-number event classes, event matrices), `wgdlife.selection`
(the WGD-anchored selection test), `wgdlife.timing` (mutation and arm-loss timing, consensus),
`wgdlife.metastasis` (matched primary–metastasis comparison), `R/` (convergence test).

## Reproducing the MSK-IMPACT analyses

The selection, timing and primary–metastasis analyses of the article run on the MSK-IMPACT 50K release, which is
public under CC BY-NC-ND 4.0 (verbatim non-commercial sharing only, no derived works), so the derived tables are not
shipped here. `preprocess/` downloads it from Zenodo and writes the compact
tables the examples read (see `preprocess/README.md`); each example then applies the package to those tables and
checks its output against the published values shipped in `examples/expected/`.

```bash
make fetch        # python preprocess/fetch_msk_impact_50k.py    4.2 GB from Zenodo, MD5-verified
make prepare      # python preprocess/prepare_msk_impact_50k.py  about 5 minutes
python examples/msk_timing.py                                   # Supplementary Table 4, MSK-IMPACT columns
python examples/msk_selection.py --cancer-type "Breast Cancer"  # Fig. 2 statistic for one cancer type
python examples/msk_primary_met.py --cancer-type "Breast Cancer" --runs 3 --nsim 200   # Fig. 4 procedure, quick look
```

| Example | What it reproduces | Check against `examples/expected/` | Runtime (8 cores) |
|---|---|---|---|
| `msk_timing.py` | timing of oncogenic mutations and arm losses relative to WGD in every WGD+ carcinoma (10,914 samples; Supplementary Table 4, MSK-IMPACT columns) | patient counts identical and P(pre-WGD) identical to two decimals (28/28 rows) | ~20 s |
| `msk_selection.py` | the WGD-anchored selection test in the primary tumours (patients with localized disease, `DISEASE_STATUS`) of one cancer type: six event classes, lineage × panel strata, the same null and seeds; `--pancancer` adds the 12-type run | O, E, P and FDR identical for every locus of the cancer type | 1–3 min per type on 8 cores; ~1 h pan-cancer |
| `msk_primary_met.py` | the matched primary-site / metastasis-site (`METASTATIC_SITE`) resampling and category rule applied to one cancer type with the article's settings | rank correlation with the cohort-wide published effect sizes (not identity: the published run matches 37 detailed types across 12 cancers at once); breast cancer at `--runs 3 --nsim 200` gives Spearman ρ ≈ 0.9 for the primary and ≈ 0.8 for the metastasis effect sizes | ~1 min quick look; ~10 min at the published settings |

## Data

| Dataset | Source | Access | Used by |
|---|---|---|---|
| TCGA allele-specific copy number, 24 tumours (GDC ASCAT3) | Genomic Data Commons | open; redistributed in `examples/tcga_wgd_example/` | `examples/run_wgd_example.py`, tests |
| MSK-IMPACT 50K (clinical sample table, FACETS gene-level copy number, mutations) | Zenodo, https://doi.org/10.5281/zenodo.18445440 | CC BY-NC-ND 4.0; downloaded by `preprocess/` | the three MSK-IMPACT examples |
| Timed-event tables of the convergence test (`fig3_convergence_events_long.tsv`, `fig3_convergence_event_universe.tsv`) | built from the MSK-IMPACT 50K release and TCGA by the companion repository `wgdlife_analysis` | not shipped: the MSK-IMPACT rows are derived data of a CC BY-NC-ND release | `R/convergence_test_run.R` |

### Data availability

The sequencing data generated in this study are being deposited in the Genome Sequence Archive for Human (GSA-Human)
under submission ID subHRA032831; access is controlled and they are not used by this package. Previously published
datasets used here were obtained from the repositories above: the MSK-IMPACT 50K cohort through Zenodo
(https://doi.org/10.5281/zenodo.18445440) and TCGA allele-specific copy-number profiles through the Genomic Data Commons.
Accession numbers, source links and access conditions for every dataset analysed in the article are listed in its
Supplementary Table 13.
Figure code is in the companion repository `wgdlife_analysis` (https://github.com/pigudog/wgdlife_analysis).


## Licence

MIT.
