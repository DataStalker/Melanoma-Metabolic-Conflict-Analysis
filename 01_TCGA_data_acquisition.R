# ==============================================================================
# FULL SCRIPT
# TCGA-SKCM DATA ACQUISITION
# ==============================================================================

# --- 1. SETUP: INSTALL AND LOAD PACKAGES ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c(
  "TCGAbiolinks", 
  "tidyverse",    
  "SummarizedExperiment",
  "EDASeq" 
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    BiocManager::install(pkg)
  }
}

library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

message("All required packages are installed and loaded.")

# --- 2. DATA ACQUISITION: QUERY AND PREPARE ---

# Define the query (this will be fast as data is downloaded)
query_rnaseq <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Metastatic")
)

# Download data if needed (will skip if already present)
GDCdownload(query = query_rnaseq)

# Prepare the data into a SummarizedExperiment object
tcga_skcm_data <- GDCprepare(query_rnaseq)

message("Data prepared into a SummarizedExperiment object.")

# --- 4. DATA PROCESSING: EXPRESSION DATA (CORRECTED METHOD) ---

# =========================================================================
# Instead of calculating TPM from counts, we directly extract the
# pre-calculated 'tpm_unstrand' assay.
# =========================================================================
message("Extracting pre-calculated TPM matrix...")
rnaseq_tpm_matrix <- assay(tcga_skcm_data, "tpm_unstrand")
message("TPM matrix successfully extracted.")
# =========================================================================

# Clean the patient ID barcodes in the TPM matrix column names
colnames(rnaseq_tpm_matrix) <- str_sub(colnames(rnaseq_tpm_matrix), 1, 12)

# --- 5. DATA PROCESSING: CLINICAL DATA ---

clinical_data <- clinical_data_raw %>%
  select(
    patient_id = patient,
    os_time = days_to_death,
    os_status = vital_status,
    days_to_last_follow_up,
    age = age_at_index,
    gender,
    tumor_stage = ajcc_pathologic_stage
  ) %>%
  mutate(
    os_status_numeric = ifelse(os_status == "Dead", 1, 0),
    os_time_cleaned = ifelse(
      is.na(os_time), 
      days_to_last_follow_up, 
      os_time
    )
  ) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  
  # =========================================================================
# We remove the old, automatically-assigned row names first.
# =========================================================================
remove_rownames() %>%
  # =========================================================================

# Now this function will work as intended.
column_to_rownames("patient_id")

message("Clinical data successfully cleaned (error fixed).")

# --- 6. FINAL ALIGNMENT: ENSURING DATA INTEGRITY ---

common_patients <- intersect(rownames(clinical_data), colnames(rnaseq_tpm_matrix))

rnaseq_tpm_matrix_aligned <- rnaseq_tpm_matrix[, common_patients]
clinical_data_aligned <- clinical_data[common_patients, ]

# Final verification: This must return TRUE!
alignment_check <- all(rownames(clinical_data_aligned) == colnames(rnaseq_tpm_matrix_aligned))

if(alignment_check){
  message("SUCCESS: Expression and clinical data are perfectly aligned.")
} else {
  stop("CRITICAL ERROR: Data alignment failed. Do not proceed.")
}

# --- 7. SAVE FINAL OBJECTS ---
saveRDS(rnaseq_tpm_matrix_aligned, file = "tcga_skcm_tpm_aligned.rds")
saveRDS(clinical_data_aligned, file = "tcga_skcm_clinical_aligned.rds")

message(paste0(
  "Processing complete! Your final data files are saved in your working directory:\n",
  "- tcga_skcm_tpm_aligned.rds (", nrow(rnaseq_tpm_matrix_aligned), " genes, ", ncol(rnaseq_tpm_matrix_aligned), " patients)\n",
  "- tcga_skcm_clinical_aligned.rds (", nrow(clinical_data_aligned), " patients, ", ncol(clinical_data_aligned), " variables)\n",
  "You are now ready for Part 2: Acquiring the GEO validation data."
))