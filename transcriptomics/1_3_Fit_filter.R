## ==========================================================================
##  Script 1_3_Fit_filter — Model fitting + Direction filtering
##  Normalization: (FX - F0) / F0  |  Fit range: F1-F4 (x = 1,2,3,4)
##  F0 is only the normalization baseline, not a fit point.
## ==========================================================================

# paths relative to this script
.file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
BASE  <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
library(dplyr)
library(tidyr)
library(vroom)

DIR_IN   <- file.path(BASE, "results", "stage_means")
DIR_FITS <- file.path(BASE, "results", "model_fits")
dir.create(DIR_FITS, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
##  1. Load data
## ---------------------------------------------------------------------------

df <- vroom(file.path(DIR_IN, "mean_expr_per_stage.tsv"), show_col_types = FALSE) %>%
  as.data.frame()

## ---------------------------------------------------------------------------
##  2. Normalize F1-F4 to F0: norm_expr = (FX - F0) / F0
##  Positive = above baseline, negative = below. Genes without a valid
##  F0 are dropped.
## ---------------------------------------------------------------------------

f0_ref <- df %>%
  dplyr::filter(stage == 0, !is.na(mean_expr), mean_expr > 0) %>%
  dplyr::select(dataset, Gene_symbol, Ensembl_ID, f0 = mean_expr)

df_norm <- df %>%
  dplyr::filter(stage %in% c(1, 2, 3, 4)) %>%
  dplyr::inner_join(f0_ref, by = c("dataset", "Gene_symbol", "Ensembl_ID")) %>%
  dplyr::mutate(norm_expr = (mean_expr - f0) / f0) %>%
  dplyr::select(-f0)

## ---------------------------------------------------------------------------
##  3. Fit function — pick the best model by AIC
##     Input: x = 1,2,3,4   |   y = (FX - F0) / F0
## ---------------------------------------------------------------------------

fit_models <- function(sub) {

  fit_data <- sub %>%
    dplyr::filter(!is.na(norm_expr))

  xs <- fit_data$stage
  ys <- fit_data$norm_expr
  n  <- nrow(fit_data)

  if (n < 4) return(NULL)

  results <- list()

  ## Linear: y = a + b*x
  tryCatch({
    m   <- lm(ys ~ xs)
    rss <- sum(resid(m)^2)
    results$linear <- list(
      model  = "linear",
      aic    = n * log(rss / n) + 2 * 3,
      r2     = summary(m)$r.squared,
      coef_a = unname(coef(m)[1]),
      coef_b = unname(coef(m)[2]),
      coef_c = NA_real_
    )
  }, error = function(e) NULL)

  ## Power-law: y = a * x^b  (x >= 1, no offset)
  tryCatch({
    a0 <- if (abs(mean(ys)) < 1e-10) 0.01 else mean(ys) / mean(xs)
    m  <- nls(ys ~ a * xs^b,
              start   = list(a = a0, b = 1),
              control = nls.control(maxiter = 500, warnOnly = TRUE))
    rss <- sum(resid(m)^2)
    results$powerlaw <- list(
      model  = "power_law",
      aic    = n * log(rss / n) + 2 * 3,
      r2     = 1 - rss / sum((ys - mean(ys))^2),
      coef_a = unname(coef(m)["a"]),
      coef_b = unname(coef(m)["b"]),
      coef_c = NA_real_
    )
  }, error = function(e) NULL)

  ## Sigmoid: y = a / (1 + exp(-b*(x-c)))
  tryCatch({
    lin_slope <- coef(lm(ys ~ xs))[2]
    b0 <- if (!is.na(lin_slope) && lin_slope != 0) sign(lin_slope) else 1
    a0 <- sign(mean(ys)) * max(abs(ys))
    m  <- nls(ys ~ a / (1 + exp(-b * (xs - c))),
              start   = list(a = a0, b = b0, c = median(xs)),
              control = nls.control(maxiter = 500, warnOnly = TRUE))
    rss <- sum(resid(m)^2)
    results$sigmoid <- list(
      model  = "sigmoid",
      aic    = n * log(rss / n) + 2 * 4,
      r2     = 1 - rss / sum((ys - mean(ys))^2),
      coef_a = unname(coef(m)["a"]),
      coef_b = unname(coef(m)["b"]),
      coef_c = unname(coef(m)["c"])
    )
  }, error = function(e) NULL)

  if (length(results) == 0) return(NULL)

  aics <- sapply(results, `[[`, "aic")
  best <- results[[which.min(aics)]]
  best$n_points <- n

  as.data.frame(best)
}

## ---------------------------------------------------------------------------
##  4. Apply per gene x dataset
## ---------------------------------------------------------------------------

fits <- df_norm %>%
  dplyr::group_by(Gene_symbol, Ensembl_ID, dataset) %>%
  dplyr::group_modify(~ fit_models(.x)) %>%
  dplyr::ungroup()

write.table(fits, file.path(DIR_FITS, "model_fits.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("model_fits.tsv guardado (", nrow(fits), " filas)")

## ---------------------------------------------------------------------------
##  5. Net direction f(4) - f(1): difference between fit-range ends
##  Positive = model rises F1 to F4 vs F0 baseline. Negative = falls.
## ---------------------------------------------------------------------------

fits <- fits %>%
  dplyr::mutate(
    delta_F4_F1 = dplyr::case_when(
      model == "linear"    ~ 3 * coef_b,
      model == "power_law" ~ coef_a * (4^coef_b) - coef_a * (1^coef_b),
      model == "sigmoid"   ~ coef_a / (1 + exp(-coef_b * (4 - coef_c))) -
                             coef_a / (1 + exp(-coef_b * (1 - coef_c))),
      TRUE ~ NA_real_
    ),
    direction = ifelse(delta_F4_F1 > 0, "up", "down")
  )

## Genes where all 3 datasets agree on direction
genes_concordant <- fits %>%
  dplyr::group_by(Gene_symbol, Ensembl_ID) %>%
  dplyr::summarise(
    n_datasets   = dplyr::n(),
    n_up         = sum(direction == "up"),
    n_down       = sum(direction == "down"),
    dominant_dir = dplyr::case_when(
      n_up   == 3 ~ "up",
      n_down == 3 ~ "down",
      TRUE        ~ "discordant"
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(dominant_dir != "discordant")

## Final table: one row per gene, columns per dataset
conserved <- fits %>%
  dplyr::inner_join(
    genes_concordant %>% dplyr::select(Gene_symbol, Ensembl_ID, dominant_dir),
    by = c("Gene_symbol", "Ensembl_ID")
  ) %>%
  tidyr::pivot_wider(
    id_cols     = c(Gene_symbol, Ensembl_ID, dominant_dir),
    names_from  = dataset,
    values_from = c(model, coef_a, coef_b, coef_c, r2, delta_F4_F1),
    names_glue  = "{dataset}_{.value}"
  )

write.table(conserved, file.path(DIR_FITS, "conserved_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("conserved_genes.tsv guardado")
message("  up:   ", sum(conserved$dominant_dir == "up"))
message("  down: ", sum(conserved$dominant_dir == "down"))
