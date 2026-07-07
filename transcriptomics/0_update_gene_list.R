# ==========================================================================
# Update gene list: add TEAD1, TEAD2, TEAD3, TEAD4, SRF + annotate Ensembl IDs
# ==========================================================================

# paths relative to this script
.file    <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
BASE     <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
DIR_DATA <- file.path(BASE, "data")

library(vroom)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --------------------------------------------------------------------------
# 1. LOAD original table
# --------------------------------------------------------------------------

genes_nucleo_all <- vroom(file.path(DIR_DATA, "genes_nucleo_all.tsv"))

# --------------------------------------------------------------------------
# 2. ADD TEAD/SRF rows
# --------------------------------------------------------------------------

new_rows <- tibble(
  Protein = c("TEAD1", "TEAD2", "TEAD3", "TEAD4", "SRF"),
  
  Gene = c("TEAD1", "TEAD2", "TEAD3", "TEAD4", "SRF"),
  
  Role_in_nuclear_mechanics = c(
    "Transcription factor mechanosensitive; mediates YAP/TAZ-dependent gene expression in response to nuclear deformation and ECM stiffness",
    
    "Core mechanotransducer transcription factor; binds YAP/TAZ to drive stiffness-induced transcriptional programs; regulates genes linked to ECM remodeling",
    
    "Transcription factor mechanosensitive; mediates YAP/TAZ-dependent gene expression in response to nuclear deformation and ECM stiffness; contributes to cytoskeleton–nucleus feedback",
    
    "Core mechanotransducer transcription factor; binds YAP/TAZ to drive stiffness-induced transcriptional programs; regulates genes linked to proliferation, survival, and ECM remodeling",
    
    "Serum Response Factor; integrates actin-cytoskeleton signals (via MRTF-A/B coactivators) to regulate nuclear mechanotransduction and cytoskeletal gene expression"
  ),
  
  General_function = c(
    "TEA domain transcription factor; DNA-binding partner of YAP/TAZ in Hippo pathway",
    
    "TEA domain transcription factor; key effector of Hippo signaling controlling cell fate and growth",
    
    "TEA domain transcription factor; DNA-binding partner of YAP/TAZ in Hippo pathway",
    
    "TEA domain transcription factor; key effector of Hippo signaling controlling cell fate and growth",
    
    "MADS-box transcription factor; master regulator of cytoskeletal and immediate-early genes; partner of MRTF-A/B"
  ),
  
  Signaling_pathways = c(
    "Hippo pathway (YAP/TAZ-TEAD axis); mechanotransduction; cytoskeleton-nucleus signaling",
    
    "Hippo pathway (YAP/TAZ-TEAD axis); mechanotransduction; ECM stiffness response",
    
    "Hippo pathway (YAP/TAZ-TEAD axis); mechanotransduction; cytoskeleton-nucleus signaling",
    
    "Hippo pathway (YAP/TAZ-TEAD axis); mechanotransduction; ECM stiffness response",
    
    "MRTF-SRF axis; actin dynamics; mechanotransduction; cytoskeleton-nucleus feedback"
  ),
  
  Epigenetics_tag = c("NO", "NO", "NO", "NO", "NO"),
  
  References_IDs = c(
    "10.1038/s41418-020-00643-5",
    "10.1016/j.cellsig.2026.112492",
    "10.1038/s41418-020-00643-5 : 10.3389/fphar.2016.00462",
    "10.1016/j.cellsig.2026.112492",
    "10.1083/jcb.200210130"
  )
)

genes_nucleo_all_updated <- bind_rows(
  genes_nucleo_all,
  new_rows
)

# --------------------------------------------------------------------------
# 3. FIX problematic symbols
# --------------------------------------------------------------------------

correction_map <- c(
  "H3F3A/H3F3B"            = "H3F3A",
  "MRTFA(MKL1)"            = "MKL1",
  "MRTFB(MKL2)"            = "MKL2",
  "PREP1(AKA PKNOX1)"      = "PKNOX1",
  "CTGF"                   = "CCN2",
  "KMT1A"                  = "SUV39H1",
  "KMT1B"                  = "SUV39H2",
  "SYNE1(DN-KASH context)" = "SYNE1"
)

genes_nucleo_all_updated <- genes_nucleo_all_updated %>%
  mutate(
    Gene_symbol = ifelse(
      Gene %in% names(correction_map),
      correction_map[Gene],
      Gene
    )
  )

# --------------------------------------------------------------------------
# 4. AUTOMATIC annotation
# --------------------------------------------------------------------------

genes_nucleo_all_updated$Ensembl_ID <- mapIds(
  org.Hs.eg.db,
  keys      = genes_nucleo_all_updated$Gene_symbol,
  keytype   = "SYMBOL",
  column    = "ENSEMBL",
  multiVals = "first"
) |> unname()

# --------------------------------------------------------------------------
# 5. MANUAL correction for unresolved genes
# --------------------------------------------------------------------------

genes_nucleo_all_updated$Ensembl_ID[
  genes_nucleo_all_updated$Gene_symbol == "H3F3A"
] <- "ENSG00000163041"

genes_nucleo_all_updated$Ensembl_ID[
  genes_nucleo_all_updated$Gene_symbol == "MKL1"
] <- "ENSG00000196588"  # MKL1/MRTFA; ENSG00000162493 was PDPN (wrong)

genes_nucleo_all_updated$Ensembl_ID[
  genes_nucleo_all_updated$Gene_symbol == "MKL2"
] <- "ENSG00000186260"  # MKL2/MRTFB; ENSG00000198959 was TGM2 (wrong)

genes_nucleo_all_updated$Ensembl_ID[
  genes_nucleo_all_updated$Gene_symbol == "MIR30C"
] <- "ENSG00000207715"

# --------------------------------------------------------------------------
# 7. SAVE
# --------------------------------------------------------------------------

genes_nucleo_all_updated <- genes_nucleo_all_updated %>%
  dplyr::distinct(Ensembl_ID, .keep_all = TRUE)

message("Genes únicos tras deduplicar: ", nrow(genes_nucleo_all_updated))

vroom_write(
  genes_nucleo_all_updated,
  file.path(DIR_DATA, "0_genes_nucleo_all_updated.tsv"),
  delim = "\t"
)
