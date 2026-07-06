# RNA-seq datasets — background

Background notes on the three public human liver RNA-seq cohorts used for
fibrosis-stage transcriptomic validation of the mechanosensitive gene panel.

This file documents **what each cohort is** (source, platform, etiology). The
processing pipeline, sample selection, and per-stage sample counts that actually
enter the fits are in [`../transcriptomics/README.md`](../transcriptomics/README.md).

## Cohorts

All three are bulk RNA-seq of snap-frozen human liver biopsies from NAFLD/NASH
patients, staged across the fibrosis spectrum (F0–F4). Together they let the
model's per-gene predictions be tested against independent, differently-sourced
transcriptional trajectories.

### GSE130970
- **Source:** [GSE130970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130970)
- **Cohort:** NAFLD severity, bulk liver RNA-seq
- **Notes:** biopsies graded for steatosis, activity, and fibrosis; spans F0–F4.

### GSE135251
- **Source:** [GSE135251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251)
- **Cohort:** large staged NAFLD/NASH cohort (Govaere et al.)
- **Notes:** one of the largest staged NAFLD liver transcriptome resources;
  fibrosis stage annotated per sample.

### GSE162694
- **Source:** [GSE162694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162694)
- **Cohort:** NAFLD fibrosis, staged biopsies
- **Notes:** independent staged cohort used to cross-check stage trajectories.

## How they are used

The datasets are downloaded, collapsed to common genes, DESeq2-normalized per
dataset (with the sex effect removed), and summarized as mean panel expression
per fibrosis stage. Each gene in the panel has a **pre-registered** response
shape (linear / power-law / sigmoid), assigned from its mechanotransduction role
*before* fitting; the pipeline keeps only genes whose measured fibrosis-stage
trajectory is direction-concordant with that prediction. This makes the
comparison a falsifiable test of the virtual-cell gene predictions rather than a
post-hoc fit.

> **Note on F0.** Here F0 means NAFLD with steatosis but *no fibrosis* — not a
> healthy baseline. Histologically normal controls are excluded from the fits.
> See the [transcriptomics pipeline](../transcriptomics/README.md) for the exact
> sample selection and the 391-sample per-stage breakdown.
