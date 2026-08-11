# ==============================================================================
# Revised analysis for the melanoma glutamine-immune preprint
#
# Major changes from v1:
#   * uses the curated TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR)
#   * excludes non-positive survival times by a prespecified rule
#   * makes standardized continuous scores the primary analyses
#   * reports median-defined groups only as descriptive visualizations
#   * uses Fisher's exact test and correctly codes response = 1
#   * uses Firth penalized logistic regression and leave-one-out AUC
#   * adjusts survival and score correlation for ABSOLUTE tumor purity
#   * repeats key analyses with MSigDB-curated Reactome/Hallmark gene sets
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GSVA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(broom)
  library(logistf)
  library(pROC)
  library(cowplot)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_dir, "revision_outputs")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
data_dir <- file.path(output_dir, "data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 11, base_family = "Arial"))
okabe_ito <- c(
  "Low Gln / High IFNG" = "#0072B2",
  "High Gln / High IFNG" = "#E69F00",
  "Low Gln / Low IFNG" = "#009E73",
  "High Gln / Low IFNG" = "#CC79A7"
)

message("Loading expression and clinical inputs...")
tcga_expr_raw <- readRDS(file.path(project_dir, "tcga_skcm_tpm_aligned.rds"))
geo_expr_raw <- readRDS(file.path(project_dir, "geo_gse91061_fpkm_aligned.rds"))
geo_clin_raw <- readRDS(file.path(project_dir, "geo_gse91061_clinical_aligned.rds"))
cdr <- read.csv(
  file.path(project_dir, "external_data", "TCGA-CDR-SupplementalTableS1.csv"),
  check.names = FALSE
) %>%
  filter(type == "SKCM")

purity_raw <- read.delim(
  file.path(project_dir, "external_data", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"),
  check.names = FALSE
)

collapse_to_symbols <- function(expr, keytype) {
  keys <- sub("\\..*$", "", rownames(expr))
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = keys,
    column = "SYMBOL",
    keytype = keytype,
    multiVals = "first"
  )
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- as.matrix(expr[keep, , drop = FALSE])
  grouped_sum <- rowsum(mat, group = unname(symbols[keep]), reorder = TRUE)
  grouped_n <- as.numeric(table(factor(unname(symbols[keep]), levels = rownames(grouped_sum))))
  grouped_sum / grouped_n
}

message("Mapping identifiers and aggregating duplicate gene symbols...")
tcga_expr <- collapse_to_symbols(tcga_expr_raw, "ENSEMBL")
geo_expr <- collapse_to_symbols(geo_expr_raw, "ENTREZID")

original_glutamine <- c(
  "NAGS", "CPS1", "OTC", "ASS1", "ASL", "ARG1", "ARG2", "GLS", "GLS2",
  "GLUD1", "GLUD2", "PYCR1", "PYCR2", "PYCR3", "ALDH18A1", "PRODH",
  "PRODH2", "SAT1", "SAT2", "SMOX", "SMYD3", "OAT"
)
original_ifng <- c(
  "PSMB10", "STAT1", "IFNG", "CXCL10", "CXCL9", "IDO1", "IFIT1", "IFIT3",
  "STAT2", "IRF1", "IRF7", "OAS1", "OAS2", "PSMB9", "TAP1", "GBP2", "WARS",
  "CMKLR1", "HLA-E", "IL2RG", "CD3D", "LAG3", "CD274"
)

msigdb <- readRDS(file.path(project_dir, "external_data", "msigdbr_human_2026.1.rds"))
curated_glutamine <- msigdb %>%
  filter(gs_name == "REACTOME_GLUTAMATE_AND_GLUTAMINE_METABOLISM") %>%
  pull(gene_symbol) %>% unique()
curated_ifng <- msigdb %>%
  filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>%
  pull(gene_symbol) %>% unique()

gene_sets <- list(
  candidate_glutamine = original_glutamine,
  candidate_ifng = original_ifng,
  reactome_glutamine = curated_glutamine,
  hallmark_ifng = curated_ifng
)

score_gsva <- function(expr, sets) {
  par <- gsvaParam(
    exprData = log2(as.matrix(expr) + 1),
    geneSets = sets,
    kcdf = "Gaussian",
    minSize = 5,
    maxSize = Inf
  )
  gsva(par, verbose = FALSE)
}

message("Calculating GSVA scores...")
tcga_scores <- score_gsva(tcga_expr, gene_sets)
geo_scores <- score_gsva(geo_expr, gene_sets)

present_counts <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame(
    signature = nm,
    defined_genes = length(unique(gene_sets[[nm]])),
    tcga_mapped = sum(unique(gene_sets[[nm]]) %in% rownames(tcga_expr)),
    geo_mapped = sum(unique(gene_sets[[nm]]) %in% rownames(geo_expr))
  )
}))
write.csv(present_counts, file.path(table_dir, "TableS1_signature_mapping.csv"), row.names = FALSE)

make_score_df <- function(score_matrix) {
  as.data.frame(t(score_matrix)) %>%
    tibble::rownames_to_column("patient_id") %>%
    mutate(
      gln_score = candidate_glutamine,
      ifng_score = candidate_ifng,
      curated_gln_score = reactome_glutamine,
      curated_ifng_score = hallmark_ifng,
      gln_z = as.numeric(scale(gln_score)),
      ifng_z = as.numeric(scale(ifng_score)),
      curated_gln_z = as.numeric(scale(curated_gln_score)),
      curated_ifng_z = as.numeric(scale(curated_ifng_score)),
      gln_group = ifelse(gln_score > median(gln_score), "High Gln", "Low Gln"),
      ifng_group = ifelse(ifng_score > median(ifng_score), "High IFNG", "Low IFNG"),
      conflict_group = factor(
        paste(gln_group, ifng_group, sep = " / "),
        levels = names(okabe_ito)
      )
    )
}

tcga_score_df <- make_score_df(tcga_scores)
geo_score_df <- make_score_df(geo_scores)

stage_clean <- function(x) {
  case_when(
    grepl("^(Stage 0|Stage I($|A|B|C)|Stage II($|A|B|C)|I/II NOS)", x) ~ "Early (0/I/II)",
    grepl("^Stage III|^Stage IV", x) ~ "Late (III/IV)",
    TRUE ~ NA_character_
  )
}

tcga_clin <- cdr %>%
  transmute(
    patient_id = bcr_patient_barcode,
    age = age_at_initial_pathologic_diagnosis,
    age10 = age / 10,
    gender = factor(tools::toTitleCase(tolower(gender)), levels = c("Female", "Male")),
    stage_raw = ajcc_pathologic_tumor_stage,
    stage = factor(stage_clean(ajcc_pathologic_tumor_stage), levels = c("Early (0/I/II)", "Late (III/IV)")),
    os_event = as.numeric(OS),
    os_time = as.numeric(OS.time),
    pfi_event = as.numeric(PFI),
    pfi_time = as.numeric(PFI.time)
  )

purity <- purity_raw %>%
  mutate(
    patient_id = substr(array, 1, 12),
    sample_code = substr(array, 14, 15),
    sample_type = case_when(
      sample_code == "01" ~ "Primary",
      sample_code == "06" ~ "Metastatic",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(patient_id) %>%
  mutate(n_purity_records = n()) %>%
  ungroup() %>%
  filter(n_purity_records == 1) %>%
  transmute(
    patient_id,
    purity = as.numeric(purity),
    sample_type = factor(sample_type, levels = c("Primary", "Metastatic"))
  )

tcga <- tcga_score_df %>%
  inner_join(tcga_clin, by = "patient_id") %>%
  left_join(purity, by = "patient_id") %>%
  mutate(
    os_eligible = !is.na(os_event) & !is.na(os_time) & os_time > 0,
    mv_eligible = os_eligible & !is.na(age10) & !is.na(gender) & !is.na(stage),
    composition_eligible = mv_eligible & !is.na(purity) & !is.na(sample_type)
  )

geo <- geo_score_df %>%
  inner_join(tibble::rownames_to_column(as.data.frame(geo_clin_raw), "patient_id"), by = "patient_id") %>%
  mutate(
    response = factor(response_binary, levels = c("Non-Responder", "Responder")),
    response_num = as.integer(response == "Responder")
  )

write.csv(tcga, file.path(data_dir, "tcga_analysis_dataset.csv"), row.names = FALSE)
write.csv(geo, file.path(data_dir, "gse91061_analysis_dataset.csv"), row.names = FALSE)

message("Running TCGA correlation and survival analyses...")
cor_pearson <- cor.test(tcga$gln_z, tcga$ifng_z, method = "pearson")
cor_spearman <- cor.test(tcga$gln_z, tcga$ifng_z, method = "spearman", exact = FALSE)

partial_df <- tcga %>% filter(!is.na(purity), !is.na(sample_type))
gln_resid <- residuals(lm(gln_z ~ purity + sample_type, data = partial_df))
ifng_resid <- residuals(lm(ifng_z ~ purity + sample_type, data = partial_df))
partial_cor <- cor.test(gln_resid, ifng_resid, method = "pearson")

curated_cor <- cor.test(tcga$curated_gln_z, tcga$curated_ifng_z, method = "pearson")

tcga_os <- tcga %>% filter(os_eligible)
tcga_mv <- tcga %>% filter(mv_eligible)
tcga_comp <- tcga %>% filter(composition_eligible)

km_fit <- survfit(Surv(os_time, os_event) ~ conflict_group, data = tcga_os)
logrank <- survdiff(Surv(os_time, os_event) ~ conflict_group, data = tcga_os)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)

cox_unadjusted <- coxph(
  Surv(os_time, os_event) ~ gln_z * ifng_z,
  data = tcga_os,
  x = TRUE
)
cox_clinical_only <- coxph(
  Surv(os_time, os_event) ~ age10 + gender + stage,
  data = tcga_mv,
  x = TRUE
)
cox_clinical_ifng <- coxph(
  Surv(os_time, os_event) ~ ifng_z + age10 + gender + stage,
  data = tcga_mv,
  x = TRUE
)
cox_clinical_ifng_gln <- coxph(
  Surv(os_time, os_event) ~ ifng_z + gln_z + age10 + gender + stage,
  data = tcga_mv,
  x = TRUE
)
cox_clinical <- coxph(
  Surv(os_time, os_event) ~ gln_z * ifng_z + age10 + gender + stage,
  data = tcga_mv,
  x = TRUE
)
cox_composition <- coxph(
  Surv(os_time, os_event) ~ gln_z * ifng_z + age10 + gender + stage + purity + sample_type,
  data = tcga_comp,
  x = TRUE
)
cox_curated <- coxph(
  Surv(os_time, os_event) ~ curated_gln_z * curated_ifng_z + age10 + gender + stage,
  data = tcga_mv,
  x = TRUE
)

extract_lrt_p <- function(smaller, larger) {
  tab <- anova(smaller, larger, test = "LRT")
  as.numeric(tab[2, grep("Pr\\(", colnames(tab))[1]])
}

incremental_models <- list(
  "Clinical covariates only" = cox_clinical_only,
  "Clinical + IFNG" = cox_clinical_ifng,
  "Clinical + IFNG + glutamine" = cox_clinical_ifng_gln,
  "Clinical + IFNG + glutamine + interaction" = cox_clinical
)
incremental_results <- bind_rows(lapply(names(incremental_models), function(label) {
  model <- incremental_models[[label]]
  conc <- summary(model)$concordance
  data.frame(
    model = label,
    n = model$n,
    events = model$nevent,
    parameters = length(coef(model)),
    AIC = AIC(model),
    concordance = unname(conc[1]),
    concordance_se = unname(conc[2])
  )
})) %>%
  mutate(
    incremental_lrt_p = c(
      NA,
      extract_lrt_p(cox_clinical_only, cox_clinical_ifng),
      extract_lrt_p(cox_clinical_ifng, cox_clinical_ifng_gln),
      extract_lrt_p(cox_clinical_ifng_gln, cox_clinical)
    )
  )
write.csv(incremental_results,
          file.path(table_dir, "TableS4_incremental_model_comparison.csv"), row.names = FALSE)

# Sensitivity analysis for the mild IFNG proportional-hazards signal: estimate
# separate IFNG associations before and after three years using a counting-
# process formulation. Median follow-up extends beyond this cut point.
tcga_mv_split <- survSplit(
  Surv(os_time, os_event) ~ .,
  data = tcga_mv,
  cut = 1095,
  episode = "time_period",
  id = "split_id",
  start = "tstart",
  end = "os_time",
  event = "os_event"
)
tcga_mv_split$time_period <- factor(
  tcga_mv_split$time_period,
  levels = c(1, 2),
  labels = c("0-3 years", ">3 years")
)
cox_time_split <- coxph(
  Surv(tstart, os_time, os_event) ~ gln_z + ifng_z * time_period +
    gln_z:ifng_z + age10 + gender + stage + cluster(patient_id),
  data = tcga_mv_split,
  x = TRUE
)
time_split_results <- broom::tidy(cox_time_split, exponentiate = TRUE, conf.int = TRUE)
write.csv(time_split_results,
          file.path(table_dir, "TableS5_time_split_IFNG_sensitivity.csv"), row.names = FALSE)

tidy_cox <- function(model, label) {
  model_n <- model$n
  model_events <- model$nevent
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(model = label, n = model_n, events = model_events) %>%
    select(model, n, events, term, estimate, conf.low, conf.high, p.value)
}

cox_results <- bind_rows(
  tidy_cox(cox_unadjusted, "Unadjusted"),
  tidy_cox(cox_clinical, "Clinical covariates"),
  tidy_cox(cox_composition, "Clinical + purity/sample type"),
  tidy_cox(cox_curated, "Curated signatures + clinical covariates")
)
write.csv(cox_results, file.path(table_dir, "Table2_cox_models.csv"), row.names = FALSE)

ph_results <- bind_rows(lapply(list(
  "Unadjusted" = cox_unadjusted,
  "Clinical covariates" = cox_clinical,
  "Clinical + purity/sample type" = cox_composition,
  "Curated signatures + clinical covariates" = cox_curated
), function(model) {
  z <- cox.zph(model)
  data.frame(term = rownames(z$table), z$table, row.names = NULL)
}), .id = "model")
write.csv(ph_results, file.path(table_dir, "TableS2_proportional_hazards_tests.csv"), row.names = FALSE)

message("Running GEO exact and penalized response analyses...")
response_table <- table(geo$conflict_group, geo$response)
fisher_result <- fisher.test(response_table)

firth_additive <- logistf(response_num ~ gln_z + ifng_z, data = geo, pl = TRUE)
firth_interaction <- logistf(response_num ~ gln_z * ifng_z, data = geo, pl = TRUE)
firth_curated <- logistf(response_num ~ curated_gln_z + curated_ifng_z, data = geo, pl = TRUE)

tidy_firth <- function(model, label) {
  data.frame(
    model = label,
    term = names(model$coefficients),
    odds_ratio = exp(model$coefficients),
    conf.low = exp(model$ci.lower),
    conf.high = exp(model$ci.upper),
    p.value = model$prob,
    row.names = NULL
  )
}

firth_results <- bind_rows(
  tidy_firth(firth_additive, "Candidate signatures - additive"),
  tidy_firth(firth_interaction, "Candidate signatures - exploratory interaction"),
  tidy_firth(firth_curated, "Curated signatures - additive")
)
write.csv(firth_results, file.path(table_dir, "Table3_firth_logistic_models.csv"), row.names = FALSE)

loocv_prob <- rep(NA_real_, nrow(geo))
for (i in seq_len(nrow(geo))) {
  fit_i <- logistf(response_num ~ gln_z + ifng_z, data = geo[-i, ], pl = FALSE)
  loocv_prob[i] <- as.numeric(predict(fit_i, newdata = geo[i, , drop = FALSE], type = "response"))
}
geo$loocv_probability <- loocv_prob
loocv_roc <- roc(geo$response_num, geo$loocv_probability, quiet = TRUE, direction = "<")
loocv_auc <- as.numeric(auc(loocv_roc))
loocv_auc_ci <- as.numeric(ci.auc(loocv_roc, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE))
write.csv(geo, file.path(data_dir, "gse91061_analysis_dataset_with_loocv.csv"), row.names = FALSE)

missingness <- data.frame(
  analysis_stage = c(
    "Expression cohort",
    "Curated OS available with positive time",
    "Clinical-covariate Cox model",
    "Purity/sample-type sensitivity model"
  ),
  n = c(nrow(tcga), nrow(tcga_os), nrow(tcga_mv), nrow(tcga_comp)),
  events = c(NA, sum(tcga_os$os_event), sum(tcga_mv$os_event), sum(tcga_comp$os_event))
)
write.csv(missingness, file.path(table_dir, "Table1_cohort_flow.csv"), row.names = FALSE)

baseline <- tcga %>%
  group_by(conflict_group) %>%
  summarise(
    N = n(),
    age_median = median(age, na.rm = TRUE),
    age_q1 = quantile(age, 0.25, na.rm = TRUE),
    age_q3 = quantile(age, 0.75, na.rm = TRUE),
    male_n = sum(gender == "Male", na.rm = TRUE),
    female_n = sum(gender == "Female", na.rm = TRUE),
    early_stage_n = sum(stage == "Early (0/I/II)", na.rm = TRUE),
    late_stage_n = sum(stage == "Late (III/IV)", na.rm = TRUE),
    stage_missing_n = sum(is.na(stage)),
    primary_n = sum(sample_type == "Primary", na.rm = TRUE),
    metastatic_n = sum(sample_type == "Metastatic", na.rm = TRUE),
    sample_type_missing_n = sum(is.na(sample_type)),
    .groups = "drop"
  )
write.csv(baseline, file.path(table_dir, "TableS3_TCGA_baseline_by_group.csv"), row.names = FALSE)

response_summary <- geo %>%
  group_by(conflict_group) %>%
  summarise(
    total = n(),
    responders = sum(response_num),
    nonresponders = sum(1 - response_num),
    response_percent = 100 * mean(response_num),
    .groups = "drop"
  )
write.csv(response_summary, file.path(table_dir, "Table4_GEO_response_by_group.csv"), row.names = FALSE)

# ---- Figures -----------------------------------------------------------------

flow <- missingness %>% mutate(y = rev(seq_len(n())))
p_flow <- ggplot(flow) +
  geom_rect(aes(xmin = 0, xmax = 8, ymin = y - 0.32, ymax = y + 0.32),
            fill = "#EEF4F8", color = "#4A6572", linewidth = 0.4) +
  geom_text(aes(x = 4, y = y, label = paste0(analysis_stage, "\nn = ", n,
                                               ifelse(is.na(events), "", paste0(", events = ", events)))),
            size = 3.4, lineheight = 1.05) +
  geom_segment(data = flow[-nrow(flow), ],
               aes(x = 4, xend = 4, y = y - 0.34, yend = y - 0.66),
               arrow = arrow(length = unit(0.10, "in")), color = "#4A6572") +
  coord_cartesian(xlim = c(0, 8), ylim = c(0.5, 4.5), clip = "off") +
  labs(title = "TCGA-SKCM analysis flow") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
ggsave(file.path(figure_dir, "Figure1_TCGA_cohort_flow.png"), p_flow,
       width = 6.5, height = 5.2, dpi = 300, bg = "white")

p_cor <- ggplot(tcga, aes(gln_z, ifng_z, color = purity)) +
  geom_point(alpha = 0.65, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#333333", linewidth = 0.7) +
  scale_color_viridis_c(option = "C", na.value = "#BDBDBD", name = "ABSOLUTE\npurity") +
  labs(
    title = "Candidate glutamine and IFNG scores in TCGA-SKCM",
    subtitle = sprintf(
      "Pearson r = %.3f (p %s); purity/sample-type-adjusted r = %.3f (p %s)",
      unname(cor_pearson$estimate), format.pval(cor_pearson$p.value, digits = 2),
      unname(partial_cor$estimate), format.pval(partial_cor$p.value, digits = 2)
    ),
    x = "Candidate glutamine score (standardized)",
    y = "Candidate IFNG score (standardized)"
  ) +
  theme(legend.position = "right")
ggsave(file.path(figure_dir, "Figure2_TCGA_score_correlation.png"), p_cor,
       width = 7.2, height = 5.4, dpi = 300, bg = "white")

km_plot <- ggsurvplot(
  km_fit,
  data = tcga_os,
  risk.table = TRUE,
  pval = paste0("Log-rank p = ", format.pval(logrank_p, digits = 3)),
  conf.int = FALSE,
  palette = unname(okabe_ito),
  xlab = "Time from diagnosis (days)",
  ylab = "Overall survival probability",
  title = "Overall survival by median-defined metabolic-immune group",
  legend.title = "Descriptive group",
  legend.labs = names(okabe_ito),
  risk.table.height = 0.28,
  ggtheme = theme_classic(base_size = 11, base_family = "Arial"),
  tables.theme = theme_cleantable(base_family = "Arial")
)
png(file.path(figure_dir, "Figure3_TCGA_Kaplan_Meier.png"),
    width = 2400, height = 2100, res = 300)
print(km_plot)
dev.off()

p_response <- ggplot(response_summary, aes(conflict_group, response_percent, fill = conflict_group)) +
  geom_col(width = 0.68, color = "#333333", linewidth = 0.3) +
  geom_text(aes(label = paste0(responders, "/", total)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" / ", "\n", x, fixed = TRUE)) +
  scale_y_continuous(limits = c(0, max(response_summary$response_percent) + 10),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "Observed anti-PD-1 response by descriptive group",
    subtitle = paste0("Fisher's exact p = ", format.pval(fisher_result$p.value, digits = 3),
                      "; labels show responders/total"),
    x = NULL,
    y = "Responder proportion"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
ggsave(file.path(figure_dir, "Figure4_GEO_response.png"), p_response,
       width = 7.2, height = 5.3, dpi = 300, bg = "white")

forest_terms <- c("gln_z", "ifng_z", "gln_z:ifng_z")
forest_labels <- c(
  gln_z = "Glutamine score (per SD)",
  ifng_z = "IFNG score (per SD)",
  `gln_z:ifng_z` = "Score interaction"
)
forest <- cox_results %>%
  filter(model %in% c("Unadjusted", "Clinical covariates", "Clinical + purity/sample type"),
         term %in% forest_terms) %>%
  mutate(term_label = forest_labels[term], model = factor(model,
    levels = c("Unadjusted", "Clinical covariates", "Clinical + purity/sample type")))
p_forest <- ggplot(forest, aes(estimate, term_label, color = model)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#777777") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.55), height = 0.16) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  scale_x_log10() +
  scale_color_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  labs(
    title = "Association of candidate scores with overall survival",
    x = "Hazard ratio (log scale)", y = NULL, color = "Model"
  ) +
  theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "Figure5_TCGA_forest.png"), p_forest,
       width = 7.2, height = 4.8, dpi = 300, bg = "white")

p_density1 <- ggplot(tcga, aes(gln_z)) +
  geom_density(fill = "#56B4E9", alpha = 0.65) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "Candidate glutamine score", x = "Standardized score", y = "Density")
p_density2 <- ggplot(tcga, aes(ifng_z)) +
  geom_density(fill = "#E69F00", alpha = 0.65) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "Candidate IFNG score", x = "Standardized score", y = "Density")
ggsave(file.path(figure_dir, "FigureS1_score_distributions.png"),
       plot_grid(p_density1, p_density2, nrow = 1),
       width = 8.5, height = 4.0, dpi = 300, bg = "white")

summary_values <- data.frame(
  metric = c(
    "tcga_expression_n", "tcga_os_n", "tcga_os_events", "tcga_mv_n", "tcga_mv_events",
    "tcga_composition_n", "tcga_composition_events", "pearson_r", "pearson_p",
    "spearman_rho", "spearman_p", "partial_r", "partial_p", "curated_r", "curated_r_p",
    "logrank_p", "geo_n", "geo_responders", "fisher_p", "loocv_auc",
    "loocv_auc_ci_low", "loocv_auc_ci_high"
  ),
  value = c(
    nrow(tcga), nrow(tcga_os), sum(tcga_os$os_event), nrow(tcga_mv), sum(tcga_mv$os_event),
    nrow(tcga_comp), sum(tcga_comp$os_event), unname(cor_pearson$estimate), cor_pearson$p.value,
    unname(cor_spearman$estimate), cor_spearman$p.value, unname(partial_cor$estimate), partial_cor$p.value,
    unname(curated_cor$estimate), curated_cor$p.value, logrank_p, nrow(geo), sum(geo$response_num),
    fisher_result$p.value, loocv_auc, loocv_auc_ci[1], loocv_auc_ci[3]
  )
)
write.csv(summary_values, file.path(output_dir, "analysis_summary_values.csv"), row.names = FALSE)

saveRDS(list(
  gene_sets = gene_sets,
  tcga_scores = tcga_scores,
  geo_scores = geo_scores,
  tcga = tcga,
  geo = geo,
  cox_unadjusted = cox_unadjusted,
  cox_clinical_only = cox_clinical_only,
  cox_clinical_ifng = cox_clinical_ifng,
  cox_clinical_ifng_gln = cox_clinical_ifng_gln,
  cox_clinical = cox_clinical,
  cox_composition = cox_composition,
  cox_curated = cox_curated,
  cox_time_split = cox_time_split,
  firth_additive = firth_additive,
  firth_interaction = firth_interaction,
  firth_curated = firth_curated
), file.path(data_dir, "revised_analysis_objects.rds"))

sink(file.path(output_dir, "analysis_report.txt"))
cat("REVISED MELANOMA ANALYSIS REPORT\n")
cat("================================\n\n")
print(missingness)
cat("\nSignature mapping:\n")
print(present_counts)
cat("\nCandidate-score correlation:\n")
print(cor_pearson)
cat("\nPurity/sample-type-adjusted residual correlation:\n")
print(partial_cor)
cat("\nCurated-signature correlation:\n")
print(curated_cor)
cat("\nLog-rank test:\n")
print(logrank)
cat("p =", logrank_p, "\n")
cat("\nCox models:\n")
print(cox_results)
cat("\nIncremental model comparison:\n")
print(incremental_results)
cat("\nProportional hazards tests:\n")
print(ph_results)
cat("\nTime-split IFNG sensitivity model:\n")
print(time_split_results)
cat("\nFisher exact test:\n")
print(fisher_result)
cat("\nFirth logistic models:\n")
print(firth_results)
cat("\nLOOCV AUC:", loocv_auc, "95% bootstrap CI", loocv_auc_ci[1], loocv_auc_ci[3], "\n")
sink()

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Revised analysis completed successfully. Outputs: ", output_dir)
