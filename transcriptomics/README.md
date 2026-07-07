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
| `1_3_Fit_filter.R` | Fit linear / power-law / sigmoid per gene (best by AIC), keep direction-concordant genes → `results/model_fits/conserved_genes.tsv` |
| `2_validate_predictions.py` | **Validation.** Cross-check the observed `conserved_genes.tsv` against the model's predictions read live from `hepatic_panel.py` / `gene_module.py` (direction + shape), with binomial tests vs chance → `results/validation_results.{tsv,json}` |

Each script resolves its paths relative to its own location, so run them from
anywhere (`Rscript transcriptomics/0_update_gene_list.R`) without editing paths.

## Validation (the falsifiable test)

Scripts 0–1_3 produce what the data *observe*; `2_validate_predictions.py` turns
that into a test of the model. It reads the model's **pre-registered** predictions
straight from the code — for each gene, an up/down **direction** (from its
mechanotransduction role) and a linear/power/sigmoid **shape** (from a fixed
category→shape rule in `src/hepatic_panel.py`, not hand-assigned per gene) — and
compares them to the observed fibrosis-stage trajectories:

```bash
python transcriptomics/2_validate_predictions.py --extended
# defaults: --conserved results/model_fits/conserved_genes.tsv,
#           --out results/validation_results, --phenotype hepatocyte
```

Because the predictions are read live from the model, the test can never silently
drift from the model it validates, and both assignments are made independently of
the observed data — so this is a genuine test, not a fit.

**Result (extended hepatic panel, 42 genes across the three cohorts):**

| Prediction | Correct | Rate | Chance | p-value |
|---|---|---|---|---|
| **Direction** (up/down) | 41/42 | 0.976 | 0.5 | 9.8 × 10⁻¹² |
| **Shape** (linear/power/sigmoid) | 14/42 | 0.333 | 0.333 | 0.56 |

The model predicts the **direction** of the stiffness response almost perfectly
and far above chance. The exact **shape** is not resolved: with only five
fibrosis stages (F0–F4) the AIC fit rarely separates saturating power-law from
sigmoid (27/42 observed fits are power-law), so shape agreement sits at chance —
reported honestly rather than over-claimed. The single direction miss is PTK2
(FAK), plausibly reflecting negative feedback in advanced fibrosis.

## Layout

```
transcriptomics/
├── data/                 gene panel and sample metadata
├── geo/
│   ├── raw_downloads/     GEO downloads (not tracked)
│   └── R_objects/         per-dataset metadata + count objects
└── results/
    ├── normalization/     per-dataset DESeq2 output (meta_filtered_*.tsv)
    ├── stage_means/       mean_expr_per_stage.tsv
    ├── model_fits/        model_fits.tsv, conserved_genes.tsv
    └── validation_results.{tsv,json}   model-vs-observed comparison
```

Small summary tables are tracked. Heavy, regenerable artifacts (raw downloads,
count matrices, `dds_*.rds`, `norm_counts_*.tsv`, `.RData`) are git-ignored and
rebuilt by rerunning the scripts.

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

Cohort background is documented in [`../data/RANseq_datasets_info.md`](../data/RANseq_datasets_info.md).
