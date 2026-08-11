# ==============================================================================
# Strengthened analysis for the melanoma glutamine–IFNG manuscript
#
# Adds:
#   1. independent disease-specific survival validation in GSE65904;
#   2. EPIC and MCP-counter cell-composition sensitivity analyses in TCGA-SKCM;
#   3. patient-level single-cell localization in GSE72056;
#   4. missingness/selection comparisons; and
#   5. response-effect precision and approximate power in GSE91061.
# ==============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && dir.exists(args[[1]])) {
  .libPaths(c(normalizePath(args[[1]], winslash = "/"), .libPaths()))
}

suppressPackageStartupMessages({
  library(Biobase)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GSVA)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(pROC)
  library(EPIC)
  library(MCPcounter)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_dir, "external_data", "strengthened")
output_dir <- file.path(project_dir, "strengthened_revision_outputs")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
data_dir <- file.path(output_dir, "data")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 11, base_family = "sans"))
blue <- "#0072B2"
orange <- "#E69F00"
green <- "#009E73"
purple <- "#CC79A7"

old <- readRDS(file.path(
  project_dir, "revision_outputs", "data", "revised_analysis_objects.rds"
))
gene_sets <- old$gene_sets
tcga_old <- old$tcga
geo_response <- old$geo

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, format(p, scientific = TRUE, digits = 3), sprintf("%.3f", p)))
}

cox_tidy <- function(fit, model_name) {
  s <- summary(fit)
  cf <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)
  data.frame(
    model = model_name,
    n = fit$n,
    events = fit$nevent,
    term = rownames(cf),
    hazard_ratio = ci$`exp(coef)`,
    conf.low = ci$`lower .95`,
    conf.high = ci$`upper .95`,
    p.value = cf$`Pr(>|z|)`,
    row.names = NULL
  )
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), p,
         width = width, height = height, dpi = 400, bg = "white")
}

collapse_tcga_symbols <- function(expr) {
  ids <- sub("\\..*$", "", rownames(expr))
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- as.matrix(expr[keep, , drop = FALSE])
  sym <- unname(symbols[keep])
  sums <- rowsum(mat, group = sym, reorder = TRUE)
  counts <- tabulate(match(sym, rownames(sums)), nbins = nrow(sums))
  sums / counts
}

score_gsva <- function(expr, sets) {
  par <- gsvaParam(
    exprData = as.matrix(expr),
    geneSets = sets,
    kcdf = "Gaussian",
    minSize = 5,
    maxSize = Inf
  )
  gsva(par, verbose = FALSE)
}

partial_correlation <- function(df, covars, label) {
  keep <- complete.cases(df[, c("gln_z", "ifng_z", covars), drop = FALSE])
  d <- df[keep, , drop = FALSE]
  rhs <- paste(covars, collapse = " + ")
  gln_r <- residuals(lm(as.formula(paste("gln_z ~", rhs)), data = d))
  ifng_r <- residuals(lm(as.formula(paste("ifng_z ~", rhs)), data = d))
  ct <- cor.test(gln_r, ifng_r, method = "pearson")
  data.frame(
    adjustment = label,
    n = nrow(d),
    correlation = unname(ct$estimate),
    conf.low = ct$conf.int[[1]],
    conf.high = ct$conf.int[[2]],
    p.value = ct$p.value
  )
}

message("1/5 Independent survival validation in GSE65904...")
gse <- readRDS(file.path(input_dir, "GSE65904_series_matrix.rds"))
gse_probe <- exprs(gse)
gse_fd <- fData(gse)
gse_pd <- pData(gse)

symbols <- trimws(as.character(gse_fd$Symbol))
symbols <- sub("[;,].*$", "", symbols)
keep <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
gse_mat <- gse_probe[keep, , drop = FALSE]
gse_symbols <- symbols[keep]
gse_gene_sum <- rowsum(gse_mat, group = gse_symbols, reorder = TRUE)
gse_gene_n <- tabulate(match(gse_symbols, rownames(gse_gene_sum)), nbins = nrow(gse_gene_sum))
gse_gene <- gse_gene_sum / gse_gene_n

gse_scores <- score_gsva(gse_gene, gene_sets)
gse_score_df <- as.data.frame(t(gse_scores)) %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(
    gln_z = as.numeric(scale(candidate_glutamine)),
    ifng_z = as.numeric(scale(candidate_ifng)),
    curated_gln_z = as.numeric(scale(reactome_glutamine)),
    curated_ifng_z = as.numeric(scale(hallmark_ifng))
  )

na_string <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "[Not Available]")] <- NA_character_
  x
}

gse_clin <- data.frame(
  sample_id = rownames(gse_pd),
  age = suppressWarnings(as.numeric(gse_pd[["age:ch1"]])),
  gender = factor(na_string(gse_pd[["gender:ch1"]])),
  stage = factor(na_string(gse_pd[["tumor stage:ch1"]])),
  tissue = factor(na_string(gse_pd[["tissue:ch1"]])),
  dss_time = suppressWarnings(as.numeric(gse_pd[["disease specific survival in days:ch1"]])),
  dss_event = suppressWarnings(as.numeric(gse_pd[["disease specific survival (1=death, 0=alive):ch1"]]))
) %>%
  mutate(
    age10 = age / 10,
    dss_eligible = dss_time > 0 & dss_event %in% c(0, 1)
  ) %>%
  left_join(gse_score_df, by = "sample_id")

gse_dss <- gse_clin %>% filter(dss_eligible)
gse_mv <- gse_dss %>%
  filter(complete.cases(age10, gender, stage, gln_z, ifng_z)) %>%
  droplevels()

gse_clinical <- coxph(Surv(dss_time, dss_event) ~ age10 + gender + stage,
                          data = gse_mv, x = TRUE)
gse_ifng <- update(gse_clinical, . ~ . + ifng_z)
gse_ifng_gln <- update(gse_ifng, . ~ . + gln_z)
gse_candidate <- coxph(
  Surv(dss_time, dss_event) ~ gln_z * ifng_z + age10 + gender + stage,
  data = gse_mv,
  x = TRUE
)
gse_curated <- coxph(
  Surv(dss_time, dss_event) ~ curated_gln_z * curated_ifng_z + age10 + gender + stage,
  data = gse_mv,
  x = TRUE
)

gse_cox_table <- bind_rows(
  cox_tidy(gse_candidate, "Candidate signatures + clinical covariates"),
  cox_tidy(gse_curated, "Curated signatures + clinical covariates")
)
write.csv(gse_cox_table, file.path(table_dir, "Table5_GSE65904_cox_models.csv"), row.names = FALSE)

gse_incremental <- bind_rows(lapply(list(
  "Clinical covariates" = gse_clinical,
  "Clinical + IFNG" = gse_ifng,
  "Clinical + IFNG + glutamine" = gse_ifng_gln,
  "Clinical + IFNG + glutamine + interaction" = gse_candidate
), function(fit) {
  data.frame(
    n = fit$n,
    events = fit$nevent,
    parameters = length(coef(fit)),
    AIC = AIC(fit),
    concordance = unname(summary(fit)$concordance[[1]])
  )
}), .id = "model")
lrt_p <- function(smaller, larger) {
  tab <- anova(smaller, larger, test = "LRT")
  as.numeric(tab[2, ncol(tab)])
}
gse_incremental$LRT_vs_previous <- c(
  NA,
  lrt_p(gse_clinical, gse_ifng),
  lrt_p(gse_ifng, gse_ifng_gln),
  lrt_p(gse_ifng_gln, gse_candidate)
)
write.csv(gse_incremental, file.path(table_dir, "TableS6_GSE65904_incremental_models.csv"), row.names = FALSE)

gse_ph <- as.data.frame(cox.zph(gse_candidate)$table) %>%
  tibble::rownames_to_column("term")
write.csv(gse_ph, file.path(table_dir, "TableS7_GSE65904_PH_tests.csv"), row.names = FALSE)

gse_mapping <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame(
    signature = nm,
    defined_genes = length(unique(gene_sets[[nm]])),
    mapped_genes = sum(unique(gene_sets[[nm]]) %in% rownames(gse_gene))
  )
}))
write.csv(gse_mapping, file.path(table_dir, "TableS8_GSE65904_signature_mapping.csv"), row.names = FALSE)

message("2/5 TCGA deconvolution with EPIC and MCP-counter...")
tcga_raw <- readRDS(file.path(project_dir, "tcga_skcm_tpm_aligned.rds"))
tcga_gene <- collapse_tcga_symbols(tcga_raw)
tcga_log <- log2(tcga_gene + 1)

epic_result <- EPIC::EPIC(bulk = tcga_gene, reference = "TRef")
epic_frac_all <- as.data.frame(epic_result$cellFractions)
if (!all(colnames(tcga_gene) %in% rownames(epic_frac_all)) &&
    all(colnames(tcga_gene) %in% colnames(epic_frac_all))) {
  epic_frac_all <- as.data.frame(t(epic_frac_all))
}
epic_frac_all$patient_id <- rownames(epic_frac_all)
epic_gof <- as.data.frame(epic_result$fit.gof) %>%
  tibble::rownames_to_column("patient_id")
epic_frac_audit <- left_join(epic_frac_all, epic_gof, by = "patient_id")
epic_frac <- epic_frac_all %>%
  filter(patient_id %in% epic_gof$patient_id[epic_gof$convergeCode == 0])

mcp_genes <- read.delim(
  file.path(input_dir, "MCPcounter_genes.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
mcp_est <- MCPcounter::MCPcounter.estimate(
  tcga_log,
  featuresType = "HUGO_symbols",
  genes = mcp_genes
)
mcp_df <- as.data.frame(t(mcp_est))
mcp_df$patient_id <- rownames(mcp_df)

write.csv(epic_frac_audit, file.path(data_dir, "TCGA_EPIC_cell_fractions_and_fit.csv"), row.names = FALSE)
write.csv(mcp_df, file.path(data_dir, "TCGA_MCPcounter_scores.csv"), row.names = FALSE)

epic_names <- setdiff(colnames(epic_frac), "patient_id")
immune_cols <- grep("Bcells|CD4|CD8|Macroph|NK", epic_names, value = TRUE, ignore.case = TRUE)
stromal_cols <- grep("CAF|Endothelial", epic_names, value = TRUE, ignore.case = TRUE)

tcga_epic <- tcga_old %>%
  left_join(epic_frac, by = "patient_id") %>%
  mutate(
    epic_immune = rowSums(across(all_of(immune_cols)), na.rm = FALSE),
    epic_stromal = rowSums(across(all_of(stromal_cols)), na.rm = FALSE),
    epic_immune_z = as.numeric(scale(epic_immune)),
    epic_stromal_z = as.numeric(scale(epic_stromal))
  )

mcp_names <- setdiff(colnames(mcp_df), "patient_id")
mcp_immune_cols <- grep(
  "T cells|CD8|Cytotoxic|B lineage|NK cells|Monocytic|dendritic|Neutrophils",
  mcp_names, value = TRUE, ignore.case = TRUE
)
mcp_stromal_cols <- grep("Endothelial|Fibroblasts", mcp_names, value = TRUE, ignore.case = TRUE)

mcp_pc_scores <- function(df, cols, stem) {
  x <- scale(as.matrix(df[, cols, drop = FALSE]))
  pc <- prcomp(x, center = FALSE, scale. = FALSE)$x[, 1]
  avg <- rowMeans(x)
  if (cor(pc, avg, use = "complete.obs") < 0) pc <- -pc
  data.frame(patient_id = df$patient_id, value = as.numeric(scale(pc))) %>%
    rename(!!stem := value)
}

mcp_immune_pc <- mcp_pc_scores(mcp_df, mcp_immune_cols, "mcp_immune_z")
mcp_stromal_pc <- mcp_pc_scores(mcp_df, mcp_stromal_cols, "mcp_stromal_z")
tcga_mcp <- tcga_old %>%
  left_join(mcp_immune_pc, by = "patient_id") %>%
  left_join(mcp_stromal_pc, by = "patient_id")

raw_cor <- cor.test(tcga_old$gln_z, tcga_old$ifng_z, method = "pearson")
cor_table <- bind_rows(
  data.frame(
    adjustment = "Unadjusted",
    n = sum(complete.cases(tcga_old$gln_z, tcga_old$ifng_z)),
    correlation = unname(raw_cor$estimate),
    conf.low = raw_cor$conf.int[[1]],
    conf.high = raw_cor$conf.int[[2]],
    p.value = raw_cor$p.value
  ),
  partial_correlation(tcga_old, c("purity", "sample_type"), "ABSOLUTE purity + sample type"),
  partial_correlation(
    tcga_epic,
    c("purity", "sample_type", "epic_immune_z", "epic_stromal_z"),
    "Purity + sample type + EPIC immune/stromal"
  ),
  partial_correlation(
    tcga_mcp,
    c("purity", "sample_type", "mcp_immune_z", "mcp_stromal_z"),
    "Purity + sample type + MCP-counter immune/stromal PCs"
  )
)
write.csv(cor_table, file.path(table_dir, "Table6_bulk_composition_partial_correlations.csv"), row.names = FALSE)

tcga_epic_mv <- tcga_epic %>%
  filter(
    os_eligible,
    complete.cases(gln_z, ifng_z, age10, gender, stage, purity, sample_type,
                   epic_immune_z, epic_stromal_z)
  ) %>%
  droplevels()
cox_epic <- coxph(
  Surv(os_time, os_event) ~ gln_z * ifng_z + age10 + gender + stage +
    purity + sample_type + epic_immune_z + epic_stromal_z,
  data = tcga_epic_mv,
  x = TRUE
)

tcga_mcp_mv <- tcga_mcp %>%
  filter(
    os_eligible,
    complete.cases(gln_z, ifng_z, age10, gender, stage, purity, sample_type,
                   mcp_immune_z, mcp_stromal_z)
  ) %>%
  droplevels()
cox_mcp <- coxph(
  Surv(os_time, os_event) ~ gln_z * ifng_z + age10 + gender + stage +
    purity + sample_type + mcp_immune_z + mcp_stromal_z,
  data = tcga_mcp_mv,
  x = TRUE
)
deconv_cox <- bind_rows(
  cox_tidy(cox_epic, "EPIC-adjusted"),
  cox_tidy(cox_mcp, "MCP-counter-adjusted")
)
write.csv(deconv_cox, file.path(table_dir, "Table7_deconvolution_adjusted_cox.csv"), row.names = FALSE)

composition_long <- bind_rows(
  data.frame(
    patient_id = tcga_epic$patient_id,
    method = "EPIC",
    immune = tcga_epic$epic_immune_z,
    stromal = tcga_epic$epic_stromal_z,
    gln_z = tcga_epic$gln_z,
    ifng_z = tcga_epic$ifng_z
  ),
  data.frame(
    patient_id = tcga_mcp$patient_id,
    method = "MCP-counter",
    immune = tcga_mcp$mcp_immune_z,
    stromal = tcga_mcp$mcp_stromal_z,
    gln_z = tcga_mcp$gln_z,
    ifng_z = tcga_mcp$ifng_z
  )
)
composition_cor <- composition_long %>%
  pivot_longer(c(gln_z, ifng_z), names_to = "signature", values_to = "score") %>%
  group_by(method, signature) %>%
  summarise(
    immune_r = cor(score, immune, use = "complete.obs"),
    stromal_r = cor(score, stromal, use = "complete.obs"),
    .groups = "drop"
  )
write.csv(composition_cor, file.path(table_dir, "TableS9_signature_deconvolution_correlations.csv"), row.names = FALSE)

message("3/5 Patient-level single-cell localization in GSE72056...")
sc_file <- file.path(input_dir, "GSE72056_melanoma_single_cell_revised_v2.txt.gz")
con <- gzfile(sc_file, open = "rt")
header_lines <- readLines(con, n = 4)
header_parts <- strsplit(header_lines, "\t", fixed = TRUE)
cell_ids <- header_parts[[1]][-1]
tumor_id <- header_parts[[2]][-1]
malignant_code <- suppressWarnings(as.integer(header_parts[[3]][-1]))
nonmalignant_code <- suppressWarnings(as.integer(header_parts[[4]][-1]))

target_genes <- unique(unlist(gene_sets, use.names = FALSE))
sc_rows <- list()
repeat {
  lines <- readLines(con, n = 200)
  if (!length(lines)) break
  for (line in lines) {
    tab_pos <- regexpr("\t", line, fixed = TRUE)[[1]]
    if (tab_pos < 1) next
    gene <- substr(line, 1, tab_pos - 1)
    if (gene %in% target_genes) {
      values <- scan(text = substr(line, tab_pos + 1, nchar(line)),
                     what = numeric(), sep = "\t", quiet = TRUE)
      if (length(values) == length(cell_ids)) sc_rows[[gene]] <- values
    }
  }
}
close(con)

sc_expr <- do.call(rbind, sc_rows)
colnames(sc_expr) <- cell_ids
rownames(sc_expr) <- names(sc_rows)
nonzero_sd <- apply(sc_expr, 1, sd, na.rm = TRUE) > 0
sc_expr <- sc_expr[nonzero_sd, , drop = FALSE]
sc_gene_z <- t(scale(t(sc_expr)))
sc_gene_z[!is.finite(sc_gene_z)] <- 0

score_cells <- function(set_name) {
  mapped <- intersect(gene_sets[[set_name]], rownames(sc_gene_z))
  colMeans(sc_gene_z[mapped, , drop = FALSE], na.rm = TRUE)
}
detect_cells <- function(set_name) {
  mapped <- intersect(gene_sets[[set_name]], rownames(sc_expr))
  colMeans(sc_expr[mapped, , drop = FALSE] > 0, na.rm = TRUE)
}

cell_type <- ifelse(
  malignant_code == 2,
  "Malignant",
  ifelse(
    malignant_code == 1,
    c("Unresolved", "T cell", "B cell", "Macrophage", "Endothelial", "CAF", "NK cell")[nonmalignant_code + 1],
    "Unresolved"
  )
)

sc_meta <- data.frame(
  cell_id = cell_ids,
  patient_id = paste0("Tumor ", tumor_id),
  malignant_code = malignant_code,
  nonmalignant_code = nonmalignant_code,
  cell_type = factor(
    cell_type,
    levels = c("Malignant", "T cell", "NK cell", "B cell", "Macrophage", "Endothelial", "CAF", "Unresolved")
  ),
  gln_score = score_cells("candidate_glutamine"),
  ifng_score = score_cells("candidate_ifng"),
  curated_gln_score = score_cells("reactome_glutamine"),
  curated_ifng_score = score_cells("hallmark_ifng"),
  gln_detected = detect_cells("candidate_glutamine"),
  ifng_detected = detect_cells("candidate_ifng")
) %>%
  filter(cell_type != "Unresolved")

sc_patient_type <- sc_meta %>%
  group_by(patient_id, cell_type) %>%
  summarise(
    n_cells = n(),
    gln_score = median(gln_score),
    ifng_score = median(ifng_score),
    curated_gln_score = median(curated_gln_score),
    curated_ifng_score = median(curated_ifng_score),
    gln_detected = median(gln_detected),
    ifng_detected = median(ifng_detected),
    .groups = "drop"
  )
write.csv(sc_patient_type, file.path(data_dir, "GSE72056_patient_celltype_scores.csv"), row.names = FALSE)

sc_counts <- sc_meta %>%
  count(cell_type, name = "cells") %>%
  left_join(
    sc_patient_type %>% count(cell_type, name = "patients"),
    by = "cell_type"
  )
write.csv(sc_counts, file.path(table_dir, "Table8_GSE72056_cell_counts.csv"), row.names = FALSE)

bootstrap_median_ci <- function(x, reps = 2000) {
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  b <- replicate(reps, median(sample(x, replace = TRUE)))
  unname(quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

paired_results <- list()
nonmal_types <- setdiff(levels(sc_meta$cell_type), c("Malignant", "Unresolved"))
for (score_name in c("gln_score", "ifng_score", "curated_gln_score", "curated_ifng_score")) {
  malignant_df <- sc_patient_type %>%
    filter(cell_type == "Malignant") %>%
    select(patient_id, malignant = all_of(score_name))
  for (ct in nonmal_types) {
    other_df <- sc_patient_type %>%
      filter(cell_type == ct) %>%
      select(patient_id, other = all_of(score_name))
    paired <- inner_join(malignant_df, other_df, by = "patient_id")
    if (nrow(paired) >= 5) {
      diff <- paired$malignant - paired$other
      ci <- bootstrap_median_ci(diff)
      wt <- wilcox.test(paired$malignant, paired$other, paired = TRUE, exact = FALSE)
      paired_results[[paste(score_name, ct)]] <- data.frame(
        signature = score_name,
        comparison = paste("Malignant vs", ct),
        paired_patients = nrow(paired),
        median_difference = median(diff),
        conf.low = ci[[1]],
        conf.high = ci[[2]],
        p.value = wt$p.value
      )
    }
  }
}
paired_table <- bind_rows(paired_results) %>%
  group_by(signature) %>%
  mutate(p.adjusted = p.adjust(p.value, method = "holm")) %>%
  ungroup()
write.csv(paired_table, file.path(table_dir, "Table9_GSE72056_paired_localization.csv"), row.names = FALSE)

immune_patient <- sc_patient_type %>%
  filter(cell_type %in% c("T cell", "NK cell", "B cell", "Macrophage")) %>%
  group_by(patient_id) %>%
  summarise(
    immune_ifng = weighted.mean(ifng_score, n_cells),
    immune_curated_ifng = weighted.mean(curated_ifng_score, n_cells),
    immune_cells = sum(n_cells),
    .groups = "drop"
  )
malignant_patient <- sc_patient_type %>%
  filter(cell_type == "Malignant") %>%
  select(patient_id, malignant_gln = gln_score,
         malignant_ifng = ifng_score,
         malignant_curated_gln = curated_gln_score,
         malignant_curated_ifng = curated_ifng_score,
         malignant_cells = n_cells)
cross_compartment <- inner_join(malignant_patient, immune_patient, by = "patient_id")

bootstrap_spearman <- function(x, y, reps = 3000) {
  n <- length(x)
  b <- replicate(reps, {
    i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  })
  unname(quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

cross_tests <- list(
  "Candidate: malignant glutamine vs immune IFNG" = c("malignant_gln", "immune_ifng"),
  "Candidate: malignant glutamine vs malignant IFNG" = c("malignant_gln", "malignant_ifng"),
  "Curated: malignant glutamine vs immune IFNG" = c("malignant_curated_gln", "immune_curated_ifng"),
  "Curated: malignant glutamine vs malignant IFNG" = c("malignant_curated_gln", "malignant_curated_ifng")
)
cross_table <- bind_rows(lapply(names(cross_tests), function(label) {
  vars <- cross_tests[[label]]
  x <- cross_compartment[[vars[[1]]]]
  y <- cross_compartment[[vars[[2]]]]
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  ci <- bootstrap_spearman(x, y)
  data.frame(
    comparison = label,
    patients = length(x),
    spearman_rho = unname(ct$estimate),
    conf.low = ci[[1]],
    conf.high = ci[[2]],
    p.value = ct$p.value
  )
}))
write.csv(cross_table, file.path(table_dir, "TableS10_GSE72056_cross_compartment_correlations.csv"), row.names = FALSE)

sc_mapping <- bind_rows(lapply(names(gene_sets), function(nm) {
  data.frame(
    signature = nm,
    defined_genes = length(unique(gene_sets[[nm]])),
    mapped_genes = sum(unique(gene_sets[[nm]]) %in% rownames(sc_expr))
  )
}))
write.csv(sc_mapping, file.path(table_dir, "TableS11_GSE72056_signature_mapping.csv"), row.names = FALSE)

message("4/5 Missingness and response precision analyses...")
summarize_continuous <- function(df, group, variable, comparison) {
  g <- df[[group]]
  x <- df[[variable]]
  a <- x[g]
  b <- x[!g]
  p <- tryCatch(wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.frame(
    comparison = comparison,
    variable = variable,
    included_n = sum(!is.na(a)),
    included_summary = sprintf("%.2f [%.2f, %.2f]", median(a, na.rm = TRUE), quantile(a, .25, na.rm = TRUE), quantile(a, .75, na.rm = TRUE)),
    excluded_n = sum(!is.na(b)),
    excluded_summary = sprintf("%.2f [%.2f, %.2f]", median(b, na.rm = TRUE), quantile(b, .25, na.rm = TRUE), quantile(b, .75, na.rm = TRUE)),
    p.value = p
  )
}

summarize_categorical <- function(df, group, variable, comparison) {
  g <- df[[group]]
  x <- addNA(factor(df[[variable]]), ifany = TRUE)
  tab <- table(x, g)
  p <- tryCatch(fisher.test(tab, simulate.p.value = any(dim(tab) > 2), B = 10000)$p.value,
                error = function(e) NA_real_)
  fmt <- function(v) paste(sprintf("%s: %d (%.1f%%)", names(v), v, 100 * v / sum(v)), collapse = "; ")
  data.frame(
    comparison = comparison,
    variable = variable,
    included_n = sum(g),
    included_summary = fmt(tab[, "TRUE"]),
    excluded_n = sum(!g),
    excluded_summary = fmt(tab[, "FALSE"]),
    p.value = p
  )
}

tcga_missing <- tcga_old %>%
  mutate(
    include_os = os_eligible,
    include_mv_within_os = ifelse(os_eligible, mv_eligible, NA)
  )
selection_table <- bind_rows(
  bind_rows(lapply(c("gln_z", "ifng_z", "age", "purity"), function(v) {
    summarize_continuous(tcga_missing, "include_os", v, "OS-eligible vs expression-only")
  })),
  bind_rows(lapply(c("gender", "stage", "sample_type"), function(v) {
    summarize_categorical(tcga_missing, "include_os", v, "OS-eligible vs expression-only")
  })),
  {
    d <- tcga_missing %>% filter(os_eligible) %>% mutate(include_mv = mv_eligible)
    bind_rows(
      bind_rows(lapply(c("gln_z", "ifng_z", "age", "purity", "os_time"), function(v) {
        summarize_continuous(d, "include_mv", v, "Primary multivariable vs OS-only")
      })),
      bind_rows(lapply(c("gender", "stage", "sample_type", "os_event"), function(v) {
        summarize_categorical(d, "include_mv", v, "Primary multivariable vs OS-only")
      }))
    )
  }
)
write.csv(selection_table, file.path(table_dir, "TableS12_TCGA_missingness_comparison.csv"), row.names = FALSE)

response_models <- read.csv(file.path(
  project_dir, "revision_outputs", "tables", "Table3_firth_logistic_models.csv"
))
response_precision <- response_models %>%
  filter(term != "(Intercept)") %>%
  transmute(
    model,
    term,
    odds_ratio,
    conf.low,
    conf.high,
    p.value,
    compatible_decrease = conf.low < 1,
    compatible_increase = conf.high > 1
  )
write.csv(response_precision, file.path(table_dir, "TableS13_GSE91061_effect_precision.csv"), row.names = FALSE)

logit_intercept <- function(beta, prevalence) {
  x <- qnorm(seq(0.00005, 0.99995, length.out = 20000))
  uniroot(function(a) mean(plogis(a + beta * x)) - prevalence,
          interval = c(-15, 15))$root
}
sim_power <- function(or_value, n = 49, prevalence = 10 / 49, reps = 1000) {
  beta <- log(or_value)
  alpha <- logit_intercept(beta, prevalence)
  hits <- replicate(reps, {
    x <- rnorm(n)
    y <- rbinom(n, 1, plogis(alpha + beta * x))
    if (length(unique(y)) < 2) return(NA_real_)
    fit <- suppressWarnings(glm(y ~ x, family = binomial()))
    null <- glm(y ~ 1, family = binomial())
    p <- pchisq(null$deviance - fit$deviance, df = 1, lower.tail = FALSE)
    p < 0.05
  })
  c(power = mean(hits, na.rm = TRUE), valid = sum(!is.na(hits)))
}

power_grid <- seq(1, 4, by = 0.25)
power_results <- t(vapply(power_grid, sim_power, numeric(2)))
power_table <- data.frame(
  odds_ratio_per_sd = power_grid,
  power = power_results[, "power"],
  valid_simulations = power_results[, "valid"]
)
detectable_or <- power_table$odds_ratio_per_sd[which(power_table$power >= 0.80)[1]]
if (!length(detectable_or)) detectable_or <- NA_real_
write.csv(power_table, file.path(table_dir, "TableS14_GSE91061_power_simulation.csv"), row.names = FALSE)

message("5/5 Figures and compact result summary...")
tcga_forest <- read.csv(file.path(project_dir, "revision_outputs", "tables", "Table2_cox_models.csv")) %>%
  filter(model == "Clinical covariates", term %in% c("gln_z", "ifng_z", "gln_z:ifng_z")) %>%
  transmute(
    cohort = "TCGA-SKCM (overall survival)",
    term,
    hazard_ratio = estimate,
    conf.low,
    conf.high,
    p.value
  )
gse_forest <- gse_cox_table %>%
  filter(model == "Candidate signatures + clinical covariates",
         term %in% c("gln_z", "ifng_z", "gln_z:ifng_z")) %>%
  transmute(cohort = "GSE65904 (disease-specific survival)", term,
            hazard_ratio, conf.low, conf.high, p.value)
forest <- bind_rows(tcga_forest, gse_forest) %>%
  mutate(
    term = recode(term,
                  gln_z = "Glutamine score",
                  ifng_z = "IFNG score",
                  `gln_z:ifng_z` = "Interaction"),
    term = factor(term, levels = c("Glutamine score", "IFNG score", "Interaction"))
  )
p_forest <- ggplot(forest, aes(hazard_ratio, term, color = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.55), height = 0.14, linewidth = 0.65) +
  geom_point(position = position_dodge(width = 0.55), size = 2.5) +
  scale_x_log10() +
  scale_color_manual(values = c(blue, orange)) +
  labs(
    title = "Continuous signatures across two melanoma survival cohorts",
    subtitle = "Hazard ratios per 1-SD increase; models adjusted for available clinical covariates",
    x = "Hazard ratio (log scale)", y = NULL, color = NULL
  ) +
  theme(legend.position = "bottom")
save_plot(p_forest, "Figure5_cross_cohort_survival_forest", 7.2, 4.4)

cor_plot_df <- cor_table %>%
  mutate(adjustment = factor(adjustment, levels = rev(adjustment)))
p_cor <- ggplot(cor_plot_df, aes(correlation, adjustment)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.16, color = blue) +
  geom_point(size = 2.7, color = blue) +
  scale_x_continuous(limits = c(-0.55, 0.1)) +
  labs(
    title = "Cell composition explains the bulk association",
    subtitle = "Pearson correlations between residualized continuous scores",
    x = "Correlation (95% CI)", y = NULL
  )
save_plot(p_cor, "Figure6_bulk_cell_composition", 7.8, 4.4)

sc_plot_df <- sc_patient_type %>%
  select(patient_id, cell_type, gln_score, ifng_score) %>%
  pivot_longer(c(gln_score, ifng_score), names_to = "signature", values_to = "score") %>%
  mutate(
    signature = recode(signature,
                       gln_score = "Candidate glutamine",
                       ifng_score = "Candidate IFNG")
  )
p_sc <- ggplot(sc_plot_df, aes(cell_type, score, color = cell_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.12) +
  geom_jitter(width = 0.13, height = 0, size = 1.8, alpha = 0.78) +
  facet_wrap(~signature, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c(
    "Malignant" = "#333333", "T cell" = blue, "NK cell" = "#56B4E9",
    "B cell" = green, "Macrophage" = orange, "Endothelial" = purple,
    "CAF" = "#D55E00"
  )) +
  labs(
    title = "Single-cell localization using the patient as the unit of analysis",
    subtitle = "Each point is a patient-specific median across cells of one annotated compartment",
    x = NULL, y = "Standardized signature score", color = NULL
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
save_plot(p_sc, "Figure7_single_cell_localization", 9.4, 4.8)

p_power <- ggplot(power_table, aes(odds_ratio_per_sd, power)) +
  geom_hline(yintercept = 0.8, linetype = 2, color = "grey55") +
  geom_line(color = purple, linewidth = 0.9) +
  geom_point(color = purple, size = 2) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "Approximate power in the 49-patient immunotherapy cohort",
    subtitle = "Monte Carlo logistic model; 10/49 response prevalence, two-sided alpha = 0.05",
    x = "True odds ratio per 1-SD predictor increase", y = "Power"
  )
save_plot(p_power, "FigureS5_GSE91061_power_curve", 6.5, 4.4)

cohort_flow <- data.frame(
  cohort = c("TCGA-SKCM", "GSE65904", "GSE91061", "GSE72056"),
  role = c("Discovery survival and deconvolution", "Independent survival validation",
           "Exploratory immunotherapy response", "Single-cell localization"),
  available = c(nrow(tcga_old), nrow(gse_clin), nrow(geo_response), length(unique(sc_meta$patient_id))),
  primary_analysis = c(sum(tcga_old$mv_eligible), nrow(gse_mv), nrow(geo_response), length(unique(sc_meta$patient_id))),
  events_or_responders = c(sum(tcga_old$os_event[tcga_old$mv_eligible]),
                           sum(gse_mv$dss_event), sum(geo_response$response_num), NA)
)
write.csv(cohort_flow, file.path(table_dir, "Table1_strengthened_cohort_flow.csv"), row.names = FALSE)

candidate_external <- gse_cox_table %>%
  filter(model == "Candidate signatures + clinical covariates")
summary_values <- data.frame(
  metric = c(
    "gse65904_available", "gse65904_dss_eligible", "gse65904_mv_n", "gse65904_mv_events",
    "gse65904_gln_hr", "gse65904_gln_p", "gse65904_ifng_hr", "gse65904_ifng_p",
    "gse65904_interaction_hr", "gse65904_interaction_p",
    "epic_partial_r", "epic_partial_p", "mcp_partial_r", "mcp_partial_p",
    "single_cell_cells", "single_cell_patients", "single_cell_mapped_gln", "single_cell_mapped_ifng",
    "response_detectable_or_80pct"
  ),
  value = c(
    nrow(gse_clin), nrow(gse_dss), nrow(gse_mv), sum(gse_mv$dss_event),
    candidate_external$hazard_ratio[candidate_external$term == "gln_z"],
    candidate_external$p.value[candidate_external$term == "gln_z"],
    candidate_external$hazard_ratio[candidate_external$term == "ifng_z"],
    candidate_external$p.value[candidate_external$term == "ifng_z"],
    candidate_external$hazard_ratio[candidate_external$term == "gln_z:ifng_z"],
    candidate_external$p.value[candidate_external$term == "gln_z:ifng_z"],
    cor_table$correlation[grepl("EPIC", cor_table$adjustment)],
    cor_table$p.value[grepl("EPIC", cor_table$adjustment)],
    cor_table$correlation[grepl("MCP-counter", cor_table$adjustment)],
    cor_table$p.value[grepl("MCP-counter", cor_table$adjustment)],
    nrow(sc_meta), length(unique(sc_meta$patient_id)),
    sc_mapping$mapped_genes[sc_mapping$signature == "candidate_glutamine"],
    sc_mapping$mapped_genes[sc_mapping$signature == "candidate_ifng"],
    detectable_or
  )
)
write.csv(summary_values, file.path(output_dir, "strengthened_analysis_summary.csv"), row.names = FALSE)

saveRDS(list(
  gse65904 = gse_clin,
  gse65904_scores = gse_scores,
  gse65904_models = list(
    clinical = gse_clinical, ifng = gse_ifng, ifng_gln = gse_ifng_gln,
    candidate = gse_candidate, curated = gse_curated
  ),
  deconvolution = list(
    epic = epic_frac_audit, mcp = mcp_df, correlations = cor_table,
    cox_epic = cox_epic, cox_mcp = cox_mcp
  ),
  single_cell = list(
    patient_celltype = sc_patient_type, paired = paired_table,
    cross_compartment = cross_table, mapping = sc_mapping
  ),
  missingness = selection_table,
  response_power = power_table
), file.path(data_dir, "strengthened_analysis_objects.rds"))

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Strengthened analysis complete: ", output_dir)
