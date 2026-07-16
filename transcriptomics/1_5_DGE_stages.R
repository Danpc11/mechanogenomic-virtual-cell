## ==========================================================================
##  Script 1_5_DGE_stages — Differential expression across fibrosis stages
##
##  All genes are tested, not only the nuclear mechanics list.
##  Samples: F0 + F1 + F2 + F3 + F4. F0 = NAFLD with steatosis and no
##  fibrosis, used as the baseline. Control_Normal samples are excluded.
##
##  Two levels of evidence are combined:
##
##  1. Joint model over the three datasets: ~ dataset + sex + group.
##     Pooling 391 samples gives the power needed for the early stages,
##     where each cohort on its own is underpowered. From this single fit
##     come the LRT (any change across stages) and the 10 pairwise
##     contrasts: F1/F2/F3/F4 vs F0, F2/F3/F4 vs F1, F3/F4 vs F2, F4 vs F3.
##
##  2. Per-dataset models: ~ sex + group, one per cohort. These are not
##     used for testing, only to record the sign of the fold change in
##     each cohort. The joint model removes the main dataset effect but
##     not cohort-specific stage effects, so a call is only trusted when
##     every cohort moves in the same direction as the joint estimate.
##
##  Outputs (results/dge_stages):
##    dge_contrasts_joint.tsv.gz      joint statistics, every gene and contrast
##    dge_contrasts_per_dataset.tsv.gz  per-cohort fold changes
##    dge_lrt_joint.tsv               joint LRT per gene
##    dge_calls.tsv                   calls passing all filters
##    dge_gene_summary.tsv            per-gene summary and class
##    dge_genes_per_contrast.tsv      calls per contrast
##    dge_intersections.tsv           set sizes
## ==========================================================================

# paths relative to this script
.file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
BASE  <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(tibble)
library(vroom)

DIR_OBJ <- file.path(BASE, "geo", "R_objects")
DIR_OUT <- file.path(BASE, "results", "dge_stages")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

PADJ_CUT <- 0.05
LFC_CUT  <- log2(1.5)

STAGES   <- c("F0", "F1", "F2", "F3", "F4")
DATASETS <- c("GSE130970", "GSE135251", "GSE162694")

## ---------------------------------------------------------------------------
##  1. Load and subset
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

keep          <- metadata$histology_group %in% STAGES
metadata_filt <- metadata[keep, ]
counts_filt   <- counts[, keep]

## drop samples without sex
ok            <- !is.na(metadata_filt$sex)
if (any(!ok)) message("Dropping ", sum(!ok), " sample(s) with missing sex")
metadata_filt <- metadata_filt[ok, ]
counts_filt   <- counts_filt[, ok]

message("Samples kept (F0-F4): ", ncol(counts_filt))
print(table(metadata_filt$dataset, metadata_filt$histology_group))

## ---------------------------------------------------------------------------
##  2. Contrast list: every pair of stages, later stage as numerator
## ---------------------------------------------------------------------------

pairs_mat <- combn(STAGES, 2)
contrasts_tbl <- data.frame(
  numerator   = pairs_mat[2, ],
  denominator = pairs_mat[1, ],
  stringsAsFactors = FALSE
)
contrasts_tbl$contrast <- paste0(contrasts_tbl$numerator, "_vs_",
                                 contrasts_tbl$denominator)
contrasts_tbl$type <- ifelse(contrasts_tbl$denominator == "F0",
                             "baseline", "between_strata")

## pull every pairwise contrast out of a fitted object
extract_contrasts <- function(dds, label) {
  lapply(seq_len(nrow(contrasts_tbl)), function(i) {
    r <- results(dds, contrast = c("group", contrasts_tbl$numerator[i],
                                            contrasts_tbl$denominator[i]))
    as.data.frame(r) %>%
      rownames_to_column("Ensembl_ID") %>%
      transmute(
        source   = label,
        contrast = contrasts_tbl$contrast[i],
        type     = contrasts_tbl$type[i],
        Ensembl_ID,
        baseMean,
        log2FC   = log2FoldChange,
        lfcSE,
        pvalue,
        padj
      )
  }) %>% bind_rows()
}

## ---------------------------------------------------------------------------
##  3. Joint model over the three datasets
## ---------------------------------------------------------------------------

message("\n=== Joint model ===")

mat_all <- as.matrix(counts_filt)
storage.mode(mat_all) <- "integer"
mat_all <- mat_all[rowSums(mat_all >= 10) >= 3, ]

cd_all <- data.frame(
  row.names = metadata_filt$sample_id,
  dataset = factor(metadata_filt$dataset, levels = DATASETS),
  sex     = factor(metadata_filt$sex, levels = c("Male", "Female")),
  group   = factor(metadata_filt$histology_group, levels = STAGES)
)

message("  Samples : ", ncol(mat_all))
message("  Genes   : ", nrow(mat_all))

dds_joint <- DESeqDataSetFromMatrix(mat_all, cd_all, ~ dataset + sex + group)

## Wald fit for the pairwise contrasts
dds_joint_wald <- DESeq(dds_joint)
joint_pairs    <- extract_contrasts(dds_joint_wald, "joint")

## LRT against the model without stage
dds_joint_lrt <- DESeq(dds_joint, test = "LRT", reduced = ~ dataset + sex)
joint_lrt <- as.data.frame(results(dds_joint_lrt)) %>%
  rownames_to_column("Ensembl_ID") %>%
  transmute(Ensembl_ID, baseMean, stat_lrt = stat,
            pvalue_lrt = pvalue, padj_lrt = padj)

## ---------------------------------------------------------------------------
##  4. Per-dataset models, used only for the direction of the fold change
## ---------------------------------------------------------------------------

run_dataset <- function(dataset_name) {

  message("\n=== ", dataset_name, " ===")

  idx     <- metadata_filt$dataset == dataset_name
  meta_ds <- metadata_filt[idx, ]
  mat_ds  <- as.matrix(counts_filt[, idx])
  storage.mode(mat_ds) <- "integer"
  mat_ds  <- mat_ds[rowSums(mat_ds >= 10) >= 3, ]

  cd_ds <- data.frame(
    row.names = meta_ds$sample_id,
    sex   = factor(meta_ds$sex, levels = c("Male", "Female")),
    group = factor(meta_ds$histology_group, levels = STAGES)
  )

  message("  Samples : ", ncol(mat_ds))
  message("  Genes   : ", nrow(mat_ds))

  dds <- DESeq(DESeqDataSetFromMatrix(mat_ds, cd_ds, ~ sex + group))
  extract_contrasts(dds, dataset_name)
}

per_dataset <- bind_rows(lapply(DATASETS, run_dataset))

## ---------------------------------------------------------------------------
##  5. Direction agreement between each cohort and the joint estimate
## ---------------------------------------------------------------------------

dir_check <- joint_pairs %>%
  select(contrast, Ensembl_ID, joint_lfc = log2FC) %>%
  inner_join(
    per_dataset %>% select(contrast, Ensembl_ID, dataset = source, log2FC),
    by = c("contrast", "Ensembl_ID")
  ) %>%
  filter(!is.na(log2FC), !is.na(joint_lfc)) %>%
  group_by(contrast, Ensembl_ID) %>%
  summarise(
    n_ds_tested   = n(),
    n_ds_same_dir = sum(sign(log2FC) == sign(joint_lfc)),
    .groups = "drop"
  ) %>%
  mutate(all_same_dir = n_ds_tested > 0 & n_ds_same_dir == n_ds_tested)

joint_pairs <- joint_pairs %>%
  left_join(dir_check, by = c("contrast", "Ensembl_ID")) %>%
  mutate(
    n_ds_tested   = replace_na(n_ds_tested, 0L),
    n_ds_same_dir = replace_na(n_ds_same_dir, 0L),
    all_same_dir  = replace_na(all_same_dir, FALSE),
    sig_joint = !is.na(padj) & padj < PADJ_CUT &
                !is.na(log2FC) & abs(log2FC) >= LFC_CUT,
    call = sig_joint & all_same_dir,
    dir  = ifelse(!call, NA_character_, ifelse(log2FC > 0, "up", "down"))
  )

## ---------------------------------------------------------------------------
##  6. Gene symbols
## ---------------------------------------------------------------------------

sym_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(joint_pairs$Ensembl_ID),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>%
  dplyr::rename(Ensembl_ID = ENSEMBL, Gene_symbol = SYMBOL) %>%
  distinct(Ensembl_ID, .keep_all = TRUE)

joint_pairs <- left_join(joint_pairs, sym_map, by = "Ensembl_ID")
joint_lrt   <- left_join(joint_lrt,   sym_map, by = "Ensembl_ID")

calls <- joint_pairs %>%
  filter(call) %>%
  select(Ensembl_ID, Gene_symbol, contrast, type, dir, baseMean,
         log2FC, lfcSE, padj, n_ds_tested, n_ds_same_dir) %>%
  arrange(contrast, padj)

message("\nCalls (joint significant and same direction in every cohort): ",
        nrow(calls))

## ---------------------------------------------------------------------------
##  7. Per-gene summary
##
##  baseline      : gene differs from F0 in at least one stage
##  between_strata: gene differs between two fibrotic stages
##  class         : intersection of both sets
## ---------------------------------------------------------------------------

base_set <- calls %>%
  filter(type == "baseline") %>%
  group_by(Ensembl_ID) %>%
  summarise(
    n_baseline_sig   = n(),
    baseline_dirs    = paste(sort(unique(dir)), collapse = "/"),
    first_stage_diff = min(as.integer(sub("^F(\\d)_vs_F0$", "\\1", contrast))),
    baseline_hits    = paste(contrast, collapse = ","),
    .groups = "drop"
  )

strata_set <- calls %>%
  filter(type == "between_strata") %>%
  group_by(Ensembl_ID) %>%
  summarise(
    n_strata_sig = n(),
    strata_dirs  = paste(sort(unique(dir)), collapse = "/"),
    strata_hits  = paste(contrast, collapse = ","),
    .groups = "drop"
  )

gene_summary <- joint_lrt %>%
  mutate(sig_lrt = !is.na(padj_lrt) & padj_lrt < PADJ_CUT) %>%
  select(Ensembl_ID, Gene_symbol, baseMean, padj_lrt, sig_lrt) %>%
  left_join(base_set,   by = "Ensembl_ID") %>%
  left_join(strata_set, by = "Ensembl_ID") %>%
  mutate(
    n_baseline_sig = replace_na(n_baseline_sig, 0L),
    n_strata_sig   = replace_na(n_strata_sig, 0L),
    class = case_when(
      n_baseline_sig > 0 & n_strata_sig > 0 ~ "baseline_and_strata",
      n_baseline_sig > 0                    ~ "baseline_only",
      n_strata_sig   > 0                    ~ "strata_only",
      TRUE                                  ~ "not_significant"
    )
  ) %>%
  arrange(desc(class == "baseline_and_strata"), padj_lrt)

## ---------------------------------------------------------------------------
##  8. Set sizes
## ---------------------------------------------------------------------------

intersections <- gene_summary %>%
  summarise(
    genes_tested        = n(),
    lrt_sig             = sum(sig_lrt),
    baseline_any        = sum(n_baseline_sig > 0),
    strata_any          = sum(n_strata_sig > 0),
    baseline_and_strata = sum(class == "baseline_and_strata"),
    baseline_only       = sum(class == "baseline_only"),
    strata_only         = sum(class == "strata_only")
  ) %>%
  pivot_longer(everything(), names_to = "set", values_to = "n_genes")

print(as.data.frame(intersections))

per_contrast <- calls %>%
  count(contrast, type, dir) %>%
  pivot_wider(names_from = dir, values_from = n, values_fill = 0)

print(as.data.frame(per_contrast))

cat("\nfirst stage differing from F0:\n")
print(table(gene_summary$first_stage_diff, useNA = "no"))

## ---------------------------------------------------------------------------
##  9. Write outputs
## ---------------------------------------------------------------------------

vroom_write(joint_pairs,  file.path(DIR_OUT, "dge_contrasts_joint.tsv.gz"))
vroom_write(per_dataset,  file.path(DIR_OUT, "dge_contrasts_per_dataset.tsv.gz"))
vroom_write(joint_lrt,    file.path(DIR_OUT, "dge_lrt_joint.tsv"))
vroom_write(calls,        file.path(DIR_OUT, "dge_calls.tsv"))
vroom_write(gene_summary, file.path(DIR_OUT, "dge_gene_summary.tsv"))
vroom_write(per_contrast, file.path(DIR_OUT, "dge_genes_per_contrast.tsv"))
vroom_write(intersections, file.path(DIR_OUT, "dge_intersections.tsv"))

message("\nDone. Files written to ", DIR_OUT)
