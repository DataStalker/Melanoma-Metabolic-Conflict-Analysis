# ==============================================================================
# FINAL & COMPLETE SCRIPT
# 
# ==============================================================================


# --- 1. SETUP: FULLY AUTOMATED PACKAGE INSTALLATION ---
message("--- Section 1: Setting up environment (automated) ---")

required_packages <- c(
  "here", "BiocManager", "AnnotationDbi", "org.Hs.eg.db", "GSVA", 
  "tidyverse", "survival", "survminer", "broom", "gtsummary", 
  "pheatmap", "cowplot", "cardx" 
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!pkg %in% c("AnnotationDbi", "org.Hs.eg.db", "GSVA", "BiocManager")) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

bioc_packages <- c("AnnotationDbi", "org.Hs.eg.db", "GSVA")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(here)
  library(GSVA)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(broom)
  library(pheatmap)
  library(cowplot)
})

output_dir <- here("publication_outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

message("All packages loaded. Outputs will be saved to 'publication_outputs'.")


# --- 2. DATA LOADING & 3. GENE SETS ---
message("\n--- Section 2 & 3: Loading Data & Defining Gene Sets ---")
tcga_expr <- readRDS(here("tcga_skcm_tpm_aligned.rds"))
tcga_clin <- readRDS(here("tcga_skcm_clinical_aligned.rds"))
geo_expr <- readRDS(here("geo_gse91061_fpkm_aligned.rds"))
geo_clin <- readRDS(here("geo_gse91061_clinical_aligned.rds"))
glutamine_genes <- c("NAGS", "CPS1", "OTC", "ASS1", "ASL", "ARG1", "ARG2", "GLS", "GLS2", "GLUD1", "GLUD2", "PYCR1", "PYCR2", "PYCR3", "ALDH18A1", "PRODH", "PRODH2", "SAT1", "SAT2", "SMOX", "SMYD3", "OAT")
ifng_genes <- c("PSMB10", "STAT1", "IFNG", "CXCL10", "CXCL9", "IDO1", "IFIT1", "IFIT3", "STAT2", "IRF1", "IRF7", "OAS1", "OAS2", "PSMB9", "TAP1", "GBP2", "WARS", "CMKLR1", "HLA-E", "IL2RG", "CD3D", "LAG3", "CD274")
gene_sets <- list(glutamine = glutamine_genes, ifng_response = ifng_genes)


# --- 4. CORE ANALYSIS ---
message("\n--- Section 4: Running core analysis engine ---")
# TCGA
ensembl_ids_core <- gsub("\\..*$", "", rownames(tcga_expr))
gene_symbols_tcga <- mapIds(org.Hs.eg.db, keys = ensembl_ids_core, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
tcga_expr_translated <- as.data.frame(tcga_expr) %>% mutate(gene_symbol = gene_symbols_tcga) %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>% summarise(across(everything(), mean), .groups = 'drop') %>% column_to_rownames("gene_symbol")
# GEO
entrez_ids_geo <- rownames(geo_expr)
gene_symbols_geo <- mapIds(org.Hs.eg.db, keys = entrez_ids_geo, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
geo_expr_translated <- as.data.frame(geo_expr) %>% mutate(gene_symbol = gene_symbols_geo) %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>% summarise(across(everything(), mean), .groups = 'drop') %>% column_to_rownames("gene_symbol")
# GSVA
tcga_gsvapar <- gsvaParam(exprData = as.matrix(log2(tcga_expr_translated + 1)), geneSets = gene_sets)
tcga_scores_matrix <- gsva(param = tcga_gsvapar, verbose = FALSE)
geo_gsvapar <- gsvaParam(exprData = as.matrix(log2(geo_expr_translated + 1)), geneSets = gene_sets)
geo_scores_matrix <- gsva(param = geo_gsvapar, verbose = FALSE)
# Master DFs
tcga_analysis_df <- as.data.frame(t(tcga_scores_matrix)) %>% rownames_to_column("patient_id") %>% inner_join(rownames_to_column(tcga_clin, "patient_id"), by = "patient_id") %>% mutate(gln_group = ifelse(glutamine > median(glutamine), "High_Gln", "Low_Gln"), ifng_group = ifelse(ifng_response > median(ifng_response), "High_IFNg", "Low_IFNg"), conflict_group = factor(paste(gln_group, ifng_group, sep = "/"), levels = c("Low_Gln/High_IFNg", "High_Gln/High_IFNg", "Low_Gln/Low_IFNg", "High_Gln/Low_IFNg")))
geo_analysis_df <- as.data.frame(t(geo_scores_matrix)) %>% rownames_to_column("patient_id") %>% inner_join(rownames_to_column(geo_clin, "patient_id"), by = "patient_id") %>% mutate(gln_group = ifelse(glutamine > median(glutamine), "High_Gln", "Low_Gln"), ifng_group = ifelse(ifng_response > median(ifng_response), "High_IFNg", "Low_IFNg"), conflict_group = factor(paste(gln_group, ifng_group, sep = "/"), levels = c("Low_Gln/High_IFNg", "High_Gln/High_IFNg", "Low_Gln/Low_IFNg", "High_Gln/Low_IFNg")), response_binary = factor(response_binary, levels = c("Responder", "Non-Responder")))
message("Core analysis complete.")


# --- 5. COHORT CHARACTERIZATION (SIMPLIFIED EXPORT) ---
message("\n--- Section 5: Characterizing cohorts and saving summary tables ---")
message(paste("\n*** FINAL COHORT SIZES FOR MANUSCRIPT ***"))
message(paste("TCGA Final Analysis Cohort: n =", nrow(tcga_analysis_df)))
print(table(tcga_analysis_df$conflict_group))
message(paste("\nGEO Final Analysis Cohort: n =", nrow(geo_analysis_df)))
print(table(geo_analysis_df$conflict_group, geo_analysis_df$response_binary))
message(paste("*** END OF SAMPLE COUNTS ***\n"))
tcga_summary <- tcga_analysis_df %>% group_by(conflict_group) %>% summarise(N = n(), Median_Age = median(age, na.rm = TRUE), Male_Count = sum(gender == "male", na.rm = TRUE), Female_Count = sum(gender == "female", na.rm = TRUE))
write.csv(tcga_summary, here(output_dir, "Table1_TCGA_Demographics.csv"), row.names = FALSE)
message("Saved: Table1_TCGA_Demographics.csv")
geo_summary <- geo_analysis_df %>% group_by(conflict_group, response_binary) %>% summarise(N = n(), .groups = 'drop') %>% pivot_wider(names_from = response_binary, values_from = N, values_fill = 0)
write.csv(geo_summary, here(output_dir, "Table2_GEO_Demographics.csv"), row.names = FALSE)
message("Saved: Table2_GEO_Demographics.csv")


# --- 6. DISCOVERY ANALYSIS (TCGA) ---
message("\n--- Section 6: Running Discovery Analysis on TCGA cohort ---")
# Figure 1: Kaplan-Meier Plot
surv_object_tcga <- Surv(time = as.numeric(tcga_analysis_df$os_time_cleaned), event = tcga_analysis_df$os_status_numeric)
surv_fit_tcga <- survfit(surv_object_tcga ~ conflict_group, data = tcga_analysis_df)
km_plot_tcga <- ggsurvplot(surv_fit_tcga, data = tcga_analysis_df, pval = TRUE, risk.table = TRUE, title = "Prognostic Value of the Metabolic Conflict Index (TCGA-SKCM)", xlab = "Time (Days)", legend.title = "Metabolic/Immune Group", legend.labs = c("Low Gln/High IFNg (Favorable)", "High Gln/High IFNg (Conflict)", "Low Gln/Low IFNg (Quiescent)", "High Gln/Low IFNg (Desert)"), palette = c("#00BA38", "#F8766D", "#619CFF", "#C77CFF"))
png(here(output_dir, "Figure1_TCGA_Survival_Plot.png"), width = 9, height = 7, units = "in", res = 300)
print(km_plot_tcga)
dev.off()
message("Saved: Figure1_TCGA_Survival_Plot.png")

# Table 3: Cox Proportional Hazards Models
cox_interaction_model <- coxph(surv_object_tcga ~ glutamine * ifng_response, data = tcga_analysis_df)
write.csv(tidy(cox_interaction_model, exponentiate = TRUE), here(output_dir, "Table3a_TCGA_Cox_Interaction_Model.csv"), row.names = FALSE)
message("Saved: Table3a_TCGA_Cox_Interaction_Model.csv")
tcga_for_mv_model <- tcga_analysis_df %>% mutate(simple_stage = str_extract(tumor_stage, "Stage [IVX]+"), simple_stage = case_when(simple_stage %in% c("Stage I", "Stage II") ~ "Early Stage (I/II)", simple_stage %in% c("Stage III", "Stage IV") ~ "Late Stage (III/IV)", TRUE ~ NA_character_), simple_stage = factor(simple_stage, levels = c("Early Stage (I/II)", "Late Stage (III/IV)"))) %>% filter(!is.na(simple_stage) & !is.na(age) & !is.na(gender))
surv_object_mv <- Surv(time = as.numeric(tcga_for_mv_model$os_time_cleaned), event = tcga_for_mv_model$os_status_numeric)
cox_multivariable_model <- coxph(surv_object_mv ~ glutamine * ifng_response + age + gender + simple_stage, data = tcga_for_mv_model)
write.csv(tidy(cox_multivariable_model, exponentiate = TRUE), here(output_dir, "Table3b_TCGA_Cox_Multivariable_Model.csv"), row.names = FALSE)
message("Saved: Table3b_TCGA_Cox_Multivariable_Model.csv")


# --- 7. VALIDATION ANALYSIS (GEO) ---
message("\n--- Section 7: Running Validation Analysis on GEO cohort ---")
chi_sq_test <- chisq.test(table(geo_analysis_df$conflict_group, geo_analysis_df$response_binary))
response_rates <- geo_analysis_df %>% count(conflict_group, response_binary) %>% group_by(conflict_group) %>% mutate(percentage = n / sum(n) * 100)
response_plot_geo <- ggplot(response_rates, aes(x = conflict_group, y = percentage, fill = response_binary)) + geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = c("Responder" = "#00BFC4", "Non-Responder" = "#F8766D"), name = "Response to Anti-PD-1") + labs(title = "Predictive Value of the Metabolic Conflict Index (GEO: GSE91061)", subtitle = paste("Chi-squared p-value =", format.pval(chi_sq_test$p.value, digits = 3)), x = "Metabolic/Immune Group", y = "Percentage of Patients") + theme_classic(base_size = 12) + scale_x_discrete(labels = function(x) str_replace(x, "/", "/\n")) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(here(output_dir, "Figure2_GEO_Response_Plot.png"), plot = response_plot_geo, width = 8, height = 6, dpi = 300)
message("Saved: Figure2_GEO_Response_Plot.png")
logistic_model <- glm(response_binary ~ glutamine * ifng_response, data = geo_analysis_df, family = "binomial")
write.csv(tidy(logistic_model, exponentiate = TRUE, conf.int = TRUE), here(output_dir, "Table4_GEO_Logistic_Regression_Model.csv"), row.names = FALSE)
message("Saved: Table4_GEO_Logistic_Regression_Model.csv")


# --- 8. SUPPORTING & BIOLOGICAL FIGURES ---
message("\n--- Section 8: Generating supporting and biological figures ---")
p_s1a <- ggplot(tcga_analysis_df, aes(x = glutamine)) + geom_density(fill="orange", alpha=0.7) + geom_vline(xintercept=median(tcga_analysis_df$glutamine), lty=2) + theme_classic() + labs(title="Glutamine Score Distribution (TCGA)")
p_s1b <- ggplot(tcga_analysis_df, aes(x = ifng_response)) + geom_density(fill="cyan", alpha=0.7) + geom_vline(xintercept=median(tcga_analysis_df$ifng_response), lty=2) + theme_classic() + labs(title="IFN-gamma Score Distribution (TCGA)")
ggsave(here(output_dir, "FigureS1_TCGA_Score_Distributions.png"), plot = plot_grid(p_s1a, p_s1b), width = 10, height = 5, dpi = 300)
message("Saved: FigureS1_TCGA_Score_Distributions.png")
p_s2 <- ggplot(tcga_analysis_df, aes(x = glutamine, y = ifng_response)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_classic() + labs(title = "Correlation of Signature Scores (TCGA)", subtitle = paste("Pearson correlation:", round(cor(tcga_analysis_df$glutamine, tcga_analysis_df$ifng_response), 3)))
ggsave(here(output_dir, "FigureS2_TCGA_Signature_Correlation.png"), plot = p_s2, width = 6, height = 5, dpi = 300)
message("Saved: FigureS2_TCGA_Signature_Correlation.png")

# Figure 3: Heatmap of Signature Gene Expression
genes_for_heatmap <- intersect(c(glutamine_genes, ifng_genes), rownames(tcga_expr_translated))
heatmap_expr <- tcga_expr_translated[genes_for_heatmap, tcga_analysis_df$patient_id]
# =========================================================================
# We explicitly call dplyr::select() to avoid confusion with other packages.
# =========================================================================
heatmap_annotation <- tcga_analysis_df %>%
  dplyr::select(patient_id, conflict_group) %>%
  column_to_rownames("patient_id")
# =========================================================================
gene_annotation <- data.frame(Pathway = c(rep("Glutamine", length(intersect(glutamine_genes, genes_for_heatmap))), rep("IFN-gamma", length(intersect(ifng_genes, genes_for_heatmap)))))
rownames(gene_annotation) <- genes_for_heatmap
sorted_patients <- tcga_analysis_df %>% arrange(conflict_group) %>% pull(patient_id)
heatmap_expr_sorted <- heatmap_expr[, sorted_patients]
heatmap_annotation_sorted <- heatmap_annotation[sorted_patients, , drop = FALSE]
png(here(output_dir, "Figure3_TCGA_Gene_Heatmap.png"), width = 10, height = 12, units = "in", res = 300)
pheatmap(heatmap_expr_sorted, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = heatmap_annotation_sorted, annotation_row = gene_annotation, scale = "row", main = "Expression of Signature Genes Across Patient Groups (TCGA)")
dev.off()
message("Saved: Figure3_TCGA_Gene_Heatmap.png")


# --- 9. FINAL SUMMARY ---
message("\n--- ANALYSIS COMPLETE ---")
message(paste("All outputs have been saved to the '", output_dir, "' folder.", sep=""))