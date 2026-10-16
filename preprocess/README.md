# preprocess/ — public input data for the MSK-IMPACT examples

The three MSK-IMPACT examples run on the MSK-IMPACT 50K release (Bandlamudi et al., *Cancer Cell* 2026), deposited at
Zenodo under https://doi.org/10.5281/zenodo.18445440 (CC BY-NC-ND 4.0). That licence permits use with attribution but
not redistribution, so no part of the release is included here: `fetch_msk_impact_50k.py` downloads the three files
below into `$WGDLIFE_DATA/msk_impact_50k/` (default `data/msk_impact_50k/`) and verifies each against the MD5 listed
by the record; `prepare_msk_impact_50k.py` then writes the compact tables the examples read.

| Release file | Size | MD5 | Provides |
|---|---|---|---|
| `data_clinical_sample.txt` | 9.9 MB | `627040924f7f88ee571b88749b830cbd` | cancer type and detailed type, OncoTree code, disease status and metastatic site, gene panel, MSI status, FACETS QC, WGD call, ploidy, purity, TMB |
| `data_facets_gene_level.txt` | 3.53 GB | `775f174f3442a7cba7e151548bed1d2d` | allele-specific FACETS segments (total and minor copy number) and the hg19 coordinates of the panel genes |
| `data_mutations_ccf.maf.txt` | 620 MB | `1249ca252de6938148679edf230177cc` | mutation calls with mutant-copy multiplicity, clonality and OncoKB annotation |

| Derived table (`derived/`) | Size | Content |
|---|---|---|
| `samples.tsv` | 8 MB | one row per sample (54,331) with the annotation columns above |
| `segments.tsv.gz` | 42 MB | 2.8 million de-duplicated segments: `SAMPLE_ID, chrom, seg_start, seg_end, seg_length, tcn, lcn, cn_state` |
| `gene_bed_hg19.tsv` | 10 KB | 491 panel genes with their release coordinates (`gene, chrom, mid`); release symbols mapped to current HGNC symbols |
| `mutations.tsv.gz` | 3 MB | 471,285 non-silent calls: `sample, gene, Variant_Classification, tcn, lcn, expected_alt_copies, clonality, ONCOGENIC` |

`target_genes_505.txt` is the union of the genes targeted by the MSK-IMPACT panel versions; it defines the gene
universe of the mutation analyses.

```bash
python preprocess/fetch_msk_impact_50k.py       # 4.2 GB download; skips files already present with the right MD5
python preprocess/prepare_msk_impact_50k.py     # about 5 minutes; streams the 3.5 GB gene-level table
```
