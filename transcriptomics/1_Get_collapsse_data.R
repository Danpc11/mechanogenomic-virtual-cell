##==========================================================================
##  Project: Liver fibrosis / NAFLD-NASH nuclear mechanics

##    SCRIPT FOR DOWNLOAD AND INTEGRATE DATASETS:

##  GSE130970  | Platform: Illumina HiSeq       |  Salmon estimated counts  double    Entrez IDs   78  samples
##  GSE135251  | Platform: Illumina NextSeq 500 |  Raw integer counts       integer   Ensembl IDs  216 samples
##  GSE162694  | Platform: Illumina RNA-seq     |  Raw integer counts       integer   Ensembl IDs  143 samples
## ==========================================================================

library(GEOquery)
library(data.table)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
# paths relative to this script
.file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
BASE  <- if (length(.file)) dirname(normalizePath(.file)) else getwd()

DIR_RAW <- file.path(BASE, "geo", "raw_downloads")
DIR_OBJ <- file.path(BASE, "geo", "R_objects")
dir.create(DIR_RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJ, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1800)

####################      DOWNLOAD      ######################################

## ---------------------------------------------------------------------------
##  GSE130970
##  Salmon tximport estimated counts | Entrez gene IDs | 78 samples
## ---------------------------------------------------------------------------

message("=== GSE130970 ===")
dir130 <- file.path(DIR_RAW, "GSE130970"); dir.create(dir130, showWarnings = FALSE)

gse130       <- getGEO("GSE130970", destdir = dir130, GSEMatrix = TRUE,
                       getGPL = FALSE, parseCharacteristics = TRUE)
GSE130970_META <- pData(phenoData(if (is.list(gse130)) gse130[[1]] else gse130))

getGEOSuppFiles("GSE130970", makeDirectory = FALSE, baseDir = dir130, fetch_files = TRUE)

# Load counts file
raw <- fread(list.files(dir130, pattern = "counts.*entrez.*\\.gz$", full.names = TRUE),
             data.table = FALSE)
rownames(raw)      <- as.character(raw[[1]])
GSE130970_COUNTS   <- as.matrix(raw[, -1])
storage.mode(GSE130970_COUNTS) <- "double"    # keep as double; Salmon estimates have decimals

# Map GSM
colnames(GSE130970_COUNTS) <- GSE130970_META$geo_accession[
  match(colnames(GSE130970_COUNTS), GSE130970_META$title)
]
# Confirm sample mapping
stopifnot(all(colnames(GSE130970_COUNTS) %in% GSE130970_META$geo_accession))
# Preview structure
str(GSE130970_COUNTS) 
#num [1:19585, 1:78] 15969 1899 100 2969 500 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:19585] "1" "10" "100" "1000" ...
#..$ : chr [1:78] "440349.1.X_1" "440350.1.X_1" "440351.1.X_4" "440352.1.X_4" ...

colnames(GSE130970_META)
head(GSE130970_COUNTS[1:4,1:4])
#440349.1.X_1 440350.1.X_1 440351.1.X_4 440352.1.X_4
#1      15968.6717  15623.28134  12255.15581  13328.01305
#10      1898.6296   1635.90151   1476.10785   1359.12886
#100      100.2519     72.35322     67.08624     70.67883
#1000    2969.2037   2997.20030   2547.50438   2625.51244

## ---------------------------------------------------------------------------
##  GSE135251
##  True raw integer counts | Ensembl gene IDs | 216 samples
## ---------------------------------------------------------------------------

message("\n=== GSE135251 ===")
dir135 <- file.path(DIR_RAW, "GSE135251"); dir.create(dir135, showWarnings = FALSE)

gse135       <- getGEO("GSE135251", destdir = dir135, GSEMatrix = TRUE,
                       getGPL = FALSE, parseCharacteristics = TRUE)
GSE135251_META <- pData(phenoData(if (is.list(gse135)) gse135[[1]] else gse135))

getGEOSuppFiles("GSE135251", makeDirectory = FALSE, baseDir = dir135, fetch_files = TRUE)

# Unpack .tar → 216 individual per-sample files
untar_dir <- file.path(dir135, "extracted"); dir.create(untar_dir, showWarnings = FALSE)
untar(list.files(dir135, pattern = "RAW\\.tar$", full.names = TRUE), exdir = untar_dir)

sample_files <- list.files(untar_dir, pattern = "\\.gz$", full.names = TRUE)
message("  ", length(sample_files), " sample files extracted")

# Read each file (2 columns: gene_id, count); drop STAR "__" summary rows
counts_list <- lapply(sample_files, function(f) {
  d <- fread(f, data.table = FALSE, header = FALSE)
  d <- d[!grepl("^__", d[[1]]), ]
  setNames(d[[2]], d[[1]])
})
sample_ids <- str_extract(basename(sample_files), "^GSM[0-9]+")

GSE135251_COUNTS <- matrix(unlist(counts_list),
                           nrow = length(names(counts_list[[1]])),
                           dimnames = list(names(counts_list[[1]]), sample_ids))
storage.mode(GSE135251_COUNTS) <- "integer"

# Align columns to metadata row order
shared <- intersect(rownames(GSE135251_META), colnames(GSE135251_COUNTS))
GSE135251_COUNTS <- GSE135251_COUNTS[, shared]
GSE135251_META   <- GSE135251_META[shared, ]
stopifnot(all(colnames(GSE135251_COUNTS) == rownames(GSE135251_META)))

str(GSE135251_COUNTS)
#int [1:64253, 1:216] 2565 0 605 315 92 96 41898 826 3373 360 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:64253] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" ...
#..$ : chr [1:216] "GSM3998167" "GSM3998168" "GSM3998169" "GSM3998170" ...

GSE135251_COUNTS[1:4, 1:3]
#                 GSM3998167 GSM3998168 GSM3998169
#ENSG00000000003       2565       2400       2391
#ENSG00000000005          0         14          0
#ENSG00000000419        605        525        709
#ENSG00000000457        315        330        329

colnames(GSE130970_META)

## ---------------------------------------------------------------------------
##  GSE162694
##  True raw integer counts | Ensembl gene IDs | 143 samples
## ---------------------------------------------------------------------------

message("\n=== GSE162694 ===")
dir162 <- file.path(DIR_RAW, "GSE162694"); dir.create(dir162, showWarnings = FALSE)

gse162       <- getGEO("GSE162694", destdir = dir162, GSEMatrix = TRUE,
                       getGPL = FALSE, parseCharacteristics = TRUE)
GSE162694_META <- pData(phenoData(if (is.list(gse162)) gse162[[1]] else gse162))

getGEOSuppFiles("GSE162694", makeDirectory = FALSE, baseDir = dir162, fetch_files = TRUE)

raw <- fread(list.files(dir162, pattern = "raw_counts.*\\.gz$", full.names = TRUE),
             data.table = FALSE)
rownames(raw)      <- as.character(raw[[1]])
GSE162694_COUNTS   <- as.matrix(raw[, -1])
storage.mode(GSE162694_COUNTS) <- "integer"

# Add sample_code column to metadata for easy joining later
GSE162694_META$sample_code <- str_extract(GSE162694_META$title, "\\S+$")

#map GSM
colnames(GSE162694_COUNTS) <- GSE162694_META$geo_accession[
  match(colnames(GSE162694_COUNTS), GSE162694_META$sample_code)]

# ensure
stopifnot(all(colnames(GSE162694_COUNTS) %in% GSE162694_META$geo_accession))


str(GSE162694_COUNTS)
#int [1:31683, 1:143] 2060 8 386 772 260 64 81897 975 4238 567 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:31683] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" ...
#..$ : chr [1:143] "GSM4957321" "GSM4957322" "GSM4957323" "GSM4957324" ...

GSE162694_COUNTS[1:4, 1:3]
#               GSM4957321 GSM4957322 GSM4957323
#ENSG00000000003     2060      1923       1487
#ENSG00000000005        8        10          7
#ENSG00000000419      386       362        260
#ENSG00000000457      772       911        640

colnames(GSE162694_META)

## ---------------------------------------------------------------------------
##  SAVE
## ---------------------------------------------------------------------------

for (nm in c("GSE130970_COUNTS", "GSE130970_META",
             "GSE135251_COUNTS", "GSE135251_META",
             "GSE162694_COUNTS", "GSE162694_META")) {
  saveRDS(get(nm), file.path(DIR_OBJ, paste0(nm, ".rds")))
  message("  saved → ", file.path(DIR_OBJ, paste0(nm, ".rds")))
}


####################      INTEGRATE      ######################################
# ===============================
#  COUNTS
# ===============================

# -------- GSE130970 --------
g130 <- GSE130970_COUNTS
map <- mapIds(
  org.Hs.eg.db,
  keys = rownames(g130),
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)
g130 <- g130[!is.na(map), ]
rownames(g130) <- map[!is.na(map)]

GSE130970_COUNTS <- rowsum(round(g130), group = rownames(g130))
storage.mode(GSE130970_COUNTS) <- "integer"

# -------- genes comunes --------
common_genes <- Reduce(intersect, list(
  rownames(GSE130970_COUNTS),
  rownames(GSE135251_COUNTS),
  rownames(GSE162694_COUNTS)
))

counts_master <- cbind(
  GSE130970_COUNTS[common_genes, ],
  GSE135251_COUNTS[common_genes, ],
  GSE162694_COUNTS[common_genes, ]
)

counts_master[1:5,1:5]
#                 GSM3758005 GSM3758006 GSM3758007 GSM3758008 GSM3758009
#ENSG00000000003       3410       2674       2577       2561       1750
#ENSG00000000005         14          2         21          3         20
#ENSG00000000419        654        504        528        515        576
#ENSG00000000457        970        898        835        921        948

dim(counts_master)
#[1] 16162   437

# ===============================
# METADATA – with true control annotation
# ===============================
#
# Background: fibrosis stage 0 (F0) does NOT always mean a healthy liver.
# ── GSE130970 (Hoang et al., Sci Rep 2019) ─────────────────────────────────
#   The cohort includes both healthy controls AND NAFLD patients spanning
#   the full histological spectrum
#   True controls = F0 samples with steatosis grade == 0 AND NAS == 0.
#   F0 samples with steatosis grade > 0 are NAFL patients (simple steatosis,
#   no fibrosis) and should NOT be used as true controls.
#
# ── GSE135251 (Govaere et al., Sci Transl Med 2020) ────────────────────────
#   206 biopsy-proven NAFLD patients (European NAFLD Registry) + 10 healthy
#   normal controls → 216 samples total.
#   True controls = samples where disease:ch1 matches "Control" AND stage == 0.
#   The extra stage == 0 guard prevents mislabeled samples from leaking into
#   Control_Normal (the source of the cross-tab anomaly seen previously).
#
# ── GSE162694 (Pantano et al., Sci Rep 2021) ────────────────────────────────
#   143 samples split into two distinct F0 groups (Table 1 of the paper):
#     • "Normal histology"  (n = 31): no NAFLD, no steatosis → TRUE controls
#     • "NAFLD fibrosis stage 0" (n = 35): steatosis present, no fibrosis
#   In the GEO metadata, fibrosis stage:ch1 == "normal" marks the true
#   controls, while fibrosis stage:ch1 == "0" marks NAFLD F0 patients.
#   true_control is derived from the raw text field (before numeric coercion),
#   so it is already guaranteed to be stage == 0 by construction.
# ───────────────────────────────────────────────────────────────────────────

## ── Verification helper: inspect disease labels before building metadata ──
message("\n--- GSE135251: unique values of disease:ch1 ---")
print(table(GSE135251_META$`disease:ch1`))
# Expected output includes a "Normal" or similar label for the 10 controls.
# Adjust the grepl() pattern in meta135 below if the label differs.

message("\n--- GSE162694: unique values of fibrosis stage:ch1 ---")
print(table(GSE162694_META$`fibrosis stage:ch1`))
# Expected: "normal" for 31 samples and "0","1","2","3","4" for the rest.

message("\n--- GSE130970: cross-tab fibrosis stage vs steatosis grade ---")
print(table(
  stage    = GSE130970_META$`fibrosis stage:ch1`,
  steatosis = GSE130970_META$`steatosis grade:ch1`
))
# Samples at stage 0 with steatosis 0 are the healthy controls.

# ── GSE130970 ──────────────────────────────────────────────────────────────
meta130 <- GSE130970_META %>%
  transmute(
    sample_id       = geo_accession,
    dataset         = "GSE130970",
    stage           = as.numeric(`fibrosis stage:ch1`),
    sex             = recode(`Sex:ch1`, "F" = "Female", "M" = "Male"),
    age             = as.numeric(`age at biopsy:ch1`),
    steatosis_grade = as.numeric(`steatosis grade:ch1`),
    nas_score       = as.numeric(`nafld activity score:ch1`),
    # FIX: true_control requires stage == 0 explicitly (redundant here but
    #      guards against unexpected metadata values in future re-downloads)
    true_control    = (as.numeric(`fibrosis stage:ch1`) == 0 &
                         as.numeric(`steatosis grade:ch1`) == 0 &
                         as.numeric(`nafld activity score:ch1`) == 0),
    disease         = NA_character_
  ) %>%
  mutate(
    histology_group = case_when(
      true_control                                       ~ "Control_Normal",
      stage == 0 & steatosis_grade == 0 & nas_score > 0 ~ NA_character_, # drop: 3 samples, no steatosis, NAS > 0
      stage == 0 & !true_control                        ~ "F0",  # 18 with steatosis, fibrosis baseline
      stage == 1                                        ~ "F1",
      stage == 2                                        ~ "F2",
      stage == 3                                        ~ "F3",
      stage == 4                                        ~ "F4",
      TRUE                                              ~ NA_character_
    )
  )

# ── GSE135251 ──────────────────────────────────────────────────────────────
meta135 <- GSE135251_META %>%
  transmute(
    sample_id       = geo_accession,
    dataset         = "GSE135251",
    stage           = as.numeric(`fibrosis stage:ch1`),
    sex             = NA_character_,
    age             = NA_real_,
    steatosis_grade = NA_real_,
    nas_score       = as.numeric(`nas score:ch1`),
    disease         = `disease:ch1`,
    # FIX: add stage == 0 guard so that any "Control"-labelled sample with
    #      stage > 0 (data inconsistency in GEO) is NOT assigned Control_Normal.
    #      Those anomalous samples fall through to their numeric stage (F1–F4).
    true_control    = (`disease:ch1` == "Control") & (as.numeric(`fibrosis stage:ch1`) == 0)
  ) %>%
  mutate(
    histology_group = case_when(
      true_control               ~ "Control_Normal",
      stage == 0 & !true_control ~ "F0",   # 38 confirmed NAFLD patients
      stage == 1                 ~ "F1",
      stage == 2                 ~ "F2",
      stage == 3                 ~ "F3",
      stage == 4                 ~ "F4",
      TRUE                       ~ NA_character_
    )
  )

meta135$sex <- ifelse(GSE135251_COUNTS["ENSG00000229807", ] >
                        GSE135251_COUNTS["ENSG00000129824", ] &
                        GSE135251_COUNTS["ENSG00000229807", ] >
                        GSE135251_COUNTS["ENSG00000012817", ], "Female", "Male")

# ── GSE162694 ──────────────────────────────────────────────────────────────
meta162 <- GSE162694_META %>%
  mutate(
    # true_control derived from raw text field before numeric coercion →
    # guaranteed stage == 0 by construction; no additional guard needed.
    true_control = str_detect(`fibrosis stage:ch1`, "(?i)normal"),
    stage_clean  = case_when(
      str_detect(`fibrosis stage:ch1`, "(?i)normal") ~ "0",
      TRUE ~ `fibrosis stage:ch1`
    )
  ) %>%
  transmute(
    sample_id       = geo_accession,
    dataset         = "GSE162694",
    stage           = as.numeric(stage_clean),
    sex             = `Sex:ch1`,
    age             = as.numeric(`age:ch1`),
    steatosis_grade = NA_real_,
    nas_score       = as.numeric(`nas score:ch1`),
    disease         = NA_character_,
    true_control    = true_control,
    histology_group = case_when(
      true_control               ~ "Control_Normal",
      stage == 0 & !true_control ~ "F0",   # 35 patients with steatosis, no fibrosis
      stage == 1                 ~ "F1",
      stage == 2                 ~ "F2",
      stage == 3                 ~ "F3",
      stage == 4                 ~ "F4",
      TRUE                       ~ NA_character_
    )
  )

metadata_master <- bind_rows(meta130, meta135, meta162)

metadata_master <- metadata_master[
  match(colnames(counts_master), metadata_master$sample_id), ]

#####CHECK####
stopifnot(all(colnames(counts_master) == metadata_master$sample_id))

dim(metadata_master)
#[1] 437   10   (sample_id, dataset, stage, sex, age,
#                steatosis_grade, nas_score, disease,
#                true_control, histology_group)

# ── Summary tables ──────────────────────────────────────────────────────────
message("\n--- Fibrosis stage × dataset ---")
print(table(metadata_master$dataset, metadata_master$stage))
#           0  1  2  3  4
#GSE130970 25 28  9 14  2
#GSE135251 46 48 54 54 14
#GSE162694 66 30 27  8 12

message("\n--- histology_group × dataset (use this to filter true controls) ---")
print(table(metadata_master$histology_group, metadata_master$dataset))
# 
# GSE130970 GSE135251 GSE162694
# Control_Normal         4        10        31
# F0                    21        38        35
# F1                    28        47        30
# F2                     9        53        27
# F3                    14        54         8
# F4                     2        14        12

# ── Validation: no Control_Normal sample should have stage > 0 ─────────────
message("\n--- Cross-check: histology_group × stage (must be 0 for Control_Normal) ---")
print(table(metadata_master$histology_group, metadata_master$stage))
# Expected: Control_Normal row must have counts ONLY in the stage == 0 column.
stopifnot(all(
  metadata_master$stage[metadata_master$histology_group == "Control_Normal" &
                          !is.na(metadata_master$histology_group)] == 0
))
message("  OK: all Control_Normal samples are at stage 0.")

message("\n--- true_control × dataset ---")
print(table(metadata_master$true_control, metadata_master$dataset))

#           GSE130970 GSE135251 GSE162694
# FALSE        74       206       112
# TRUE          4        10        31
# 

message("\n--- sex × dataset ---")
print(table(metadata_master$dataset, metadata_master$sex))
#           Female Male
#GSE130970     48   30
#GSE135251     93  123
#GSE162694    103   40

# ── Usage examples ───────────────────────────────────────────────────────────
# Filter to TRUE controls only:
#   ctrl_meta   <- metadata_master[metadata_master$true_control, ]
#   ctrl_counts <- counts_master[, ctrl_meta$sample_id]
#
# Filter to NAFLD F0 only (steatosis, no fibrosis — NOT controls):
#   f0_meta     <- metadata_master[metadata_master$histology_group == "F0", ]
#
# Filter to all samples EXCLUDING true controls:
#   nafld_meta  <- metadata_master[!metadata_master$true_control, ]

#################       SAVE      #############################
write.table(
  counts_master,
  file = file.path(DIR_OBJ, "1_1_counts_master.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA)

write.table(
  metadata_master,
  file = file.path(DIR_OBJ, "1_1_metadata_master.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

save.image(file.path(DIR_OBJ, "1_1_Data_Image.RData"))

#load(file.path(DIR_OBJ, "1_1_Data_Image.RData"))