# R/ — the convergence test (Fig. 3g; Supplementary Tables 5 and 6)

`convergence_test.R` implements the estimators (cancer-type-explained variance R² of rarefied event spectra =
Euclidean PERMANOVA R², normalized mutual information, leave-one-out nearest-centroid classification), the
within-sample bootstrap and the within-sample permutation of the pre/post label. `convergence_test_run.R` is the
driver that produced the published runs (main analysis, confounder controls, TCGA replication); it reads the two
timed-event tables `fig3_convergence_events_long.tsv` and `fig3_convergence_event_universe.tsv`, which the companion
repository `wgdlife_analysis` (https://github.com/pigudog/wgdlife_analysis) builds from the MSK-IMPACT 50K release and
TCGA (not shipped here: the MSK-IMPACT rows are derived data of a CC BY-NC-ND release), and writes `fig3_convergence_main.tsv`.

Input format (tab-separated). `fig3_convergence_events_long.tsv`: one row per timed event with the columns
`sample`, `cancer_type`, `event`, `kind` (`arm_loss` | `snv` | `cnloh_state`, the last only in the all-events control run), `tier` (`pre` | `post`), `platform`
(`msk` | `tcga`), `currency` (`strict` | `all`). `fig3_convergence_event_universe.tsv`: one row per event of the
universe with `event`, `pre`, `post` (counts), `platform`, `currency` (plus the bookkeeping columns `type_ok` and
`in_matched_universe`). Both are read at `source()` time, so the
engine needs the tables to load.

Outputs, written next to the inputs: `fig3_convergence_main.tsv` (every statistic per run), `fig3_convergence_bytype.tsv`
(per cancer type), `fig3_convergence_byevent.tsv` (per event) and `fig3_convergence_primary_draws.tsv` (bootstrap and null draws
of the main run). The full driver takes 8–14 min on 8 cores.
`convergence_test.R` runs the draws with `parallel::mclapply` on `NCORE` cores (8 by default; edit the constant; no
forking on Windows) under the L'Ecuyer RNG stream, so results are reproducible for a given core count.

Run with R 4.4.3 and data.table 1.18.6.1 from the repository root, with the two tables in `data/` (or `WGDLIFE_DATA`
set to their directory). `Rscript R/requirements.R` installs data.table 1.18.6.1 only when that version is absent;
building it from source needs a C compiler (Rtools on Windows, the Xcode command-line tools on macOS):

```bash
Rscript R/convergence_test_run.R
```
