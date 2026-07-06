## ==========================================================================
##  DESeq2 — Per-dataset normalization
##
##  Samples kept: F0 + F1 + F2 + F3 + F4  (F0 = NAFLD, no fibrosis, baseline)
##  Objective: normalize counts across the full fibrosis spectrum (F0–F4)
##  and export normalized count matrices for downstream per-stage averaging.
##
##  Outputs (per dataset):
##    dds_GSE130970.rds         
##    dds_GSE135251.rds
##    dds_GSE162694.rds
##    norm_counts_GSE130970.tsv  
##    norm_counts_GSE135251.tsv
##    norm_counts_GSE162694.tsv
## ==========================================================================
setwd("~/Cell-nucleus-mechanical-and-transcriptomic-modeling/")
library(DESeq2)
library(dplyr)
library(limma)

DIR_OBJ <- "./1_GEO_data/R_objects"
DIR_OUT <- "./2_DGE/2_2_per_dataset_dds"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
##  1. LOAD
## ---------------------------------------------------------------------------

counts <- read.table(
  file.path(DIR_OBJ, "1_1_counts_master.tsv"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)

metadata <- read.table(
  file.path(DIR_OBJ, "1_1_metadata_master.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

stopifnot(all(colnames(counts) == metadata$sample_id))

table(metadata$histology_group, metadata$stage)

## ---------------------------------------------------------------------------
##  2. SUBSET: keep F0 through F4 (drop Control_Normal)
##  F0 = NAFLD patients with steatosis but no fibrosis.
## ---------------------------------------------------------------------------

valid_groups <- c("F0", "F1", "F2", "F3", "F4")
keep         <- metadata$histology_group %in% valid_groups

metadata_filt <- metadata[keep, ]
counts_filt   <- counts[, keep]

message("Samples after subsetting (F0–F4): ", ncol(counts_filt))
print(table(metadata_filt$dataset, metadata_filt$histology_group))

## ---------------------------------------------------------------------------
##  3. HELPER: build and run DESeq2, then extract normalized counts
##  Design ~ sex + stage; sex is a covariate, stage numeric 0-4.
##  Output = log2(counts / size factors + 1) with sex effect removed.
## ---------------------------------------------------------------------------

run_deseq2 <- function(meta_sub, counts_all, dataset_name, dir_out) {
  
  message("\n=== ", dataset_name, " ===")
  
  ## -- subset to this dataset ------------------------------------------------
  idx      <- meta_sub$dataset == dataset_name
  meta_ds  <- meta_sub[idx, ]
  mat_ds   <- as.matrix(counts_all[, idx])
  storage.mode(mat_ds) <- "integer"
  
  ## -- drop samples with missing sex ----------------------------------------
  ok       <- !is.na(meta_ds$sex)
  if (any(!ok)) message("  Dropping ", sum(!ok), " sample(s) with missing sex")
  meta_ds  <- meta_ds[ok, ]
  mat_ds   <- mat_ds[, ok]
  
  ## -- low-count filter: keep genes with >= 10 counts in >= 3 samples -------
  mat_ds   <- mat_ds[rowSums(mat_ds >= 10) >= 3, ]
  
  message("  Samples : ", ncol(mat_ds))
  message("  Genes   : ", nrow(mat_ds))
  print(table(meta_ds$histology_group))
  
  ## -- colData ---------------------------------------------------------------
  #  F0 samples carry stage == 0 and sit at the intercept.
  coldata <- data.frame(
    row.names = meta_ds$sample_id,
    sex       = factor(meta_ds$sex,   levels = c("Male", "Female")),
    stage     = factor(meta_ds$stage, levels = c(0, 1, 2, 3, 4))
  )
  
  ## -- DESeq2 ----------------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = mat_ds,
    colData   = coldata,
    design    = ~ sex + stage
  )
  dds <- DESeq(dds)
  
  ## -- Save dds object -------------------------------------------------------
  rds_path <- file.path(dir_out, paste0("dds_", dataset_name, ".rds"))
  saveRDS(dds, rds_path)
  message("  Saved → ", rds_path)
  
  ## -- Extract normalized counts (counts / size factors) --------------------
  norm_mat <- counts(dds, normalized = TRUE)   # genes × samples, double

  ## -- Log2 transform + remove sex effect (limma::removeBatchEffect) --------
  #  Regress out sex so per-stage means track fibrosis, not sex balance.
  #  Log2(x+1) first, as removeBatchEffect assumes additive normal error.
  log_mat           <- log2(norm_mat + 1)
  log_mat_corrected <- removeBatchEffect(log_mat, batch = coldata$sex)

  tsv_path <- file.path(dir_out, paste0("norm_counts_", dataset_name, ".tsv"))
  write.table(
    log_mat_corrected,
    file  = tsv_path,
    sep   = "\t",
    quote = FALSE,
    col.names = NA
  )
  message("  Saved → ", tsv_path)

  ## -- Return list for optional downstream use in this session --------------
  invisible(list(dds = dds, norm_counts = log_mat_corrected, meta = meta_ds))
}

## ---------------------------------------------------------------------------
##  4. RUN PER DATASET
## ---------------------------------------------------------------------------

res_130 <- run_deseq2(metadata_filt, counts_filt, "GSE130970", DIR_OUT)
res_135 <- run_deseq2(metadata_filt, counts_filt, "GSE135251", DIR_OUT)
res_162 <- run_deseq2(metadata_filt, counts_filt, "GSE162694", DIR_OUT)

## ---------------------------------------------------------------------------
##  5. QUICK SANITY CHECK
##     Verify that column names in each normalized matrix match the metadata
##     stored inside the corresponding dds object.
## ---------------------------------------------------------------------------

for (res in list(res_130, res_135, res_162)) {
  dds_cols  <- colnames(res$norm_counts)
  meta_rows <- rownames(colData(res$dds))
  stopifnot(all(dds_cols == meta_rows))
}
message("\nAll column–metadata alignments OK.")

## ---------------------------------------------------------------------------
##  6. SAVE SESSION IMAGE (optional)
## ---------------------------------------------------------------------------

save.image(file.path(DIR_OUT, "2_2_per_dataset_dds.RData"))
message("\nSession image saved → ", file.path(DIR_OUT, "2_2_per_dataset_dds.RData"))

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
##  7. SAVE FILTERED METADATA PER DATASET
##     Reflects exact samples present in each normalized counts matrix
## ---------------------------------------------------------------------------
load(file.path(DIR_OUT, "2_2_per_dataset_dds.RData"))


for (res in list(res_130, res_135, res_162)) {
  ds_name  <- unique(res$meta$dataset)
  out_path <- file.path(DIR_OUT, paste0("meta_filtered_", ds_name, ".tsv"))
  write.table(
    res$meta,
    file      = out_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  message("  Saved → ", out_path)
  
  # Cross-check: counts columns == meta rows
  stopifnot(all(colnames(res$norm_counts) == res$meta$sample_id))
  message("  OK: column–metadata alignment confirmed for ", ds_name)
}
