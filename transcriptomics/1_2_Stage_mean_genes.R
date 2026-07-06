## ==========================================================================
##  Script 2_3 — Mean normalized expression per stage × dataset
## ==========================================================================

setwd("~/Cell-nucleus-mechanical-and-transcriptomic-modeling/")
library(vroom)
library(dplyr)
library(tidyr)
library(tibble)

DIR_DDS <- "./2_DGE/2_2_per_dataset_dds"
DIR_OUT <- "./2_DGE/2_3_mean_per_stage"
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
##  1. Gene panel — collapse duplicates by Ensembl_ID
## ---------------------------------------------------------------------------

genes_interest <- vroom("./0_genes_nucleo_all_updated.tsv", show_col_types = FALSE)
genes_interest <- as.data.frame(genes_interest)   # forzar data.frame plano

# Check for duplicates
dups <- dplyr::count(genes_interest, Ensembl_ID) %>% dplyr::filter(n > 1)
if (nrow(dups) > 0) {
  message("Ensembl_IDs duplicados en genes_interest:")
  print(dups)
} else {
  message("Sin duplicados en Ensembl_ID.")
}

# Collapse anyway, safe even without duplicates
genes_ref <- genes_interest %>%
  dplyr::group_by(Ensembl_ID) %>%
  dplyr::summarise(
    Gene_symbol               = dplyr::first(Gene_symbol),
    Role_in_nuclear_mechanics = paste(unique(Role_in_nuclear_mechanics), collapse = " | "),
    .groups = "drop"
  ) %>%
  as.data.frame()

target_genes <- genes_ref$Ensembl_ID
message("Genes en panel: ", length(target_genes))

## ---------------------------------------------------------------------------
##  2. Helper
## ---------------------------------------------------------------------------

mean_by_stage <- function(dataset_name) {
  
  message("\n--- Procesando ", dataset_name, " ---")
  
  ## -- load counts -----------------------------------------------------------
  norm_raw <- vroom(
    file.path(DIR_DDS, paste0("norm_counts_", dataset_name, ".tsv")),
    show_col_types = FALSE
  )
  norm_raw        <- as.data.frame(norm_raw)
  colnames(norm_raw)[1] <- "Ensembl_ID"
  rownames(norm_raw)    <- norm_raw$Ensembl_ID
  norm_raw$Ensembl_ID   <- NULL
  
  ## -- load meta -------------------------------------------------------------
  meta <- vroom(
    file.path(DIR_DDS, paste0("meta_filtered_", dataset_name, ".tsv")),
    show_col_types = FALSE
  )
  meta <- as.data.frame(meta)
  
  stopifnot(all(colnames(norm_raw) == meta$sample_id))
  
  ## -- filter genes ----------------------------------------------------------
  present <- intersect(target_genes, rownames(norm_raw))
  missing <- setdiff(target_genes, rownames(norm_raw))
  message("  Genes encontrados : ", length(present), " / ", length(target_genes))
  if (length(missing) > 0) {
    syms <- genes_ref$Gene_symbol[genes_ref$Ensembl_ID %in% missing]
    message("  Genes ausentes    : ", paste(syms, collapse = ", "))
  }
  
  norm_sub <- norm_raw[present, , drop = FALSE]
  
  ## -- long + means ----------------------------------------------------------
  #  norm_counts_*.tsv holds log2(counts/sizefactor + 1) with sex effect
  #  removed. Average directly on the log2 scale.
  long <- norm_sub %>%
    tibble::rownames_to_column("Ensembl_ID") %>%
    tidyr::pivot_longer(-Ensembl_ID, names_to = "sample_id", values_to = "log2_expr") %>%
    dplyr::left_join(
      dplyr::select(meta, sample_id, histology_group, stage),
      by = "sample_id"
    ) %>%
    dplyr::group_by(Ensembl_ID, histology_group, stage) %>%
    dplyr::summarise(
      mean_expr = mean(log2_expr),
      sd_expr   = sd(log2_expr),
      n_samples = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(dataset = dataset_name)
  
  return(long)
}

## ---------------------------------------------------------------------------
##  3. Run and combine
## ---------------------------------------------------------------------------

results <- dplyr::bind_rows(
  mean_by_stage("GSE130970"),
  mean_by_stage("GSE135251"),
  mean_by_stage("GSE162694")
) %>%
  dplyr::left_join(genes_ref, by = "Ensembl_ID") %>%
  dplyr::arrange(dataset, Gene_symbol, stage) %>%
  dplyr::select(
    dataset, Ensembl_ID, Gene_symbol,
    histology_group, stage, mean_expr, sd_expr, n_samples
  )

## ---------------------------------------------------------------------------
##  4. Preview and save
## ---------------------------------------------------------------------------

message("\n--- Preview LMNA ---")
print(dplyr::filter(results, Gene_symbol == "LMNA"), n = 15)

message("\n--- Filas por dataset × grupo ---")
print(table(results$dataset, results$histology_group))

write.table(
  results,
  file      = file.path(DIR_OUT, "mean_expr_per_stage.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("\nGuardado → ", file.path(DIR_OUT, "mean_expr_per_stage.tsv"))
message("Dimensiones: ", nrow(results), " filas × ", ncol(results), " columnas")