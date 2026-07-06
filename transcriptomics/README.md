# Transcriptomics pipeline

R pipeline that builds the mechanosensitive gene panel and fits its
fibrosis-stage expression trajectories from three human liver RNA-seq cohorts
(GSE130970, GSE135251, GSE162694). It runs from the uploaded gene list through
the per-gene model fits used to validate the virtual-cell model.

## Scripts (run in order)

| Script | Purpose |
|---|---|
| `0_update_gene_list.R` | Add TEAD1–4 / SRF, fix symbols, annotate Ensembl IDs → `data/0_genes_nucleo_all_updated.tsv` |
| `1_Get_collapsse_data.R` | Download the three GEO datasets, collapse to common genes, build the count and metadata masters |
| `1_1_dds.R` | Per-dataset DESeq2 normalization over F0–F4, sex effect removed |
| `1_2_Stage_mean_genes.R` | Mean panel expression per fibrosis stage × dataset |
| `1_3_Fit_filter.R` | Fit linear / power-law / sigmoid per gene (best by AIC), keep direction-concordant genes |

## Sample selection

Only fibrosis-graded NAFLD/NASH biopsies enter the fit. Histologically normal
controls (`Control_Normal`) are excluded, and F0 here means NAFLD with steatosis
but no fibrosis — not a healthy baseline. Samples with missing sex are dropped
before DESeq2, so **not every F0 biopsy is used** (e.g. GSE130970 keeps 18 of 21).

`data/samples_used_metadata.tsv` lists the **391 samples** that survive filtering
and feed the DESeq2 fits:

| Dataset | F0 | F1 | F2 | F3 | F4 | Total |
|---|---|---|---|---|---|---|
| GSE130970 | 18 | 28 | 9 | 14 | 2 | 71 |
| GSE135251 | 38 | 48 | 54 | 54 | 14 | 208 |
| GSE162694 | 35 | 30 | 27 | 8 | 12 | 112 |
| **All** | **91** | **106** | **90** | **76** | **28** | **391** |

Columns: `sample_id`, `dataset`, `stage`, `sex`, `age`, `steatosis_grade`,
`nas_score`, `true_control`, `disease`, `histology_group`.

## Notes

- Scripts still `setwd()` to the original project directory and use relative
  `1_GEO_data/` and `2_DGE/` output paths; adjust those to run here.
- Cohort background is documented in [`../data/RANseq_datasets_info.md`](../data/RANseq_datasets_info.md).
