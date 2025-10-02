# ==============================================================================
# DEFINITIVE DATA CORRECTION SCRIP
# ==============================================================================
# The SOLE PURPOSE of this script is to REGENERATE the two GEO .rds files correctly.
# This version incorporates the fix for the non-standard response categories in GSE91061.

# --- 1. SETUP & LOAD ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

library(GEOquery)
library(tidyverse)
library(data.table)

# --- 2. LOAD RAW DATA ---
geo_id <- "GSE91061" 

gse <- getGEO(geo_id, GSEMatrix = TRUE, getGPL = FALSE)
geo_clinical_raw <- pData(gse[[1]])
expression_file <- "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv"
if (!file.exists(expression_file)) { stop("Expression file not found.") }
geo_fpkm_matrix_raw_dt <- fread(expression_file)
setnames(geo_fpkm_matrix_raw_dt, "V1", "gene_symbol")
geo_fpkm_matrix_dt <- geo_fpkm_matrix_raw_dt[, lapply(.SD, mean), by = gene_symbol]
geo_fpkm_matrix <- geo_fpkm_matrix_dt %>% column_to_rownames("gene_symbol")

# --- 3. DEFINITIVE CLINICAL CLEANING ---
message("Applying the definitive cleaning logic to the clinical data...")
geo_clinical <- geo_clinical_raw %>%
  as_tibble() %>% 
  dplyr::select(
    patient_id = title, 
    timepoint_raw = `visit (pre or on treatment):ch1`, 
    response_raw = `response:ch1`
  ) %>%
  filter(timepoint_raw == "Pre") %>%
  
  # =========================================================================
# THE CRITICAL FIX IS HERE
# Based on our diagnostic, we now know that Responders (CR+PR) are
# coded as a single category "PRCR". We update the logic accordingly.
# =========================================================================
mutate(
  response_binary = case_when(
    response_raw == "PRCR" ~ "Responder",           # CORRECT: Look for "PRCR"
    response_raw %in% c("PD", "SD") ~ "Non-Responder", # CORRECT: This was already fine
    TRUE ~ NA_character_                             # CORRECT: This will exclude "UNK"
  )
) %>%
  # =========================================================================

dplyr::select(patient_id, response_binary) %>%
  filter(!is.na(response_binary)) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  column_to_rownames("patient_id")

# --- 4. FINAL ALIGNMENT ---
common_samples_geo <- intersect(rownames(geo_clinical), colnames(geo_fpkm_matrix))
geo_clinical_aligned <- geo_clinical[common_samples_geo, , drop = FALSE] 
geo_fpkm_matrix_aligned <- geo_fpkm_matrix[, common_samples_geo]

# --- 5. SAVE THE CORRECTED FILES ---
saveRDS(geo_clinical_aligned, file = "geo_gse91061_clinical_aligned.rds")
saveRDS(geo_fpkm_matrix_aligned, file = "geo_gse91061_fpkm_aligned.rds")

# --- 6. FINAL PROOF ---
message("\nCORRECTION COMPLETE. The .rds files have been regenerated correctly.")
message("Here is the final proof that the saved files are correct:")
final_clin_check <- readRDS("geo_gse91061_clinical_aligned.rds")
final_expr_check <- readRDS("geo_gse91061_fpkm_aligned.rds")
message("\nDimensions of the saved clinical data:")
print(dim(final_clin_check))
message("\nDimensions of the saved expression data:")
print(dim(final_expr_check))

# You can also check the response breakdown directly
message("\nBreakdown of patient responses in the final clinical file:")
print(table(final_clin_check$response_binary))

message("\nYou can now run the main analysis script again.")
