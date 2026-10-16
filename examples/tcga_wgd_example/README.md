# Public example: WGD calling and arm-loss timing on 24 TCGA tumours

`segments/` holds 24 open-access GDC ASCAT3 allele-specific copy-number profiles (TCGA-OV and TCGA-BRCA; file names
= TCGA sample barcode), 188 KB in total. `expected.tsv` gives, for each, the WGD status and the major-CN ≥ 2
fraction published for the discovery cohort of the article. `run_wgd_example.py` recomputes both from the segment
files with `wgdlife.cohort.wgd_label` and then applies the arm-level allele-state timing rule
(`wgdlife.timing.arm_state`) to the recurrent arm losses of the WGD+ tumours.
