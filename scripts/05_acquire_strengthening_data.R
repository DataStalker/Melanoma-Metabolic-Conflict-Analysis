# ==============================================================================
# Acquire public datasets used to strengthen the revised melanoma analysis.
# Sources: NCBI GEO GSE65904 (bulk microarray) and GSE72056 (single-cell RNA-seq)
# ==============================================================================

options(stringsAsFactors = FALSE, timeout = 1800)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && dir.exists(args[[1]])) {
  .libPaths(c(normalizePath(args[[1]], winslash = "/"), .libPaths()))
}

suppressPackageStartupMessages({
  library(GEOquery)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
external_dir <- file.path(project_dir, "external_data", "strengthened")
dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)

gse65904_rds <- file.path(external_dir, "GSE65904_series_matrix.rds")
if (!file.exists(gse65904_rds)) {
  message("Downloading GSE65904 series matrix and sample metadata...")
  gse65904 <- getGEO("GSE65904", GSEMatrix = TRUE, AnnotGPL = FALSE,
                     destdir = external_dir)
  saveRDS(gse65904[[1]], gse65904_rds)
}

gpl10558_rds <- file.path(external_dir, "GPL10558_platform.rds")
if (!file.exists(gpl10558_rds)) {
  message("Downloading GPL10558 platform annotation...")
  gpl10558 <- getGEO("GPL10558", AnnotGPL = FALSE, destdir = external_dir)
  saveRDS(gpl10558, gpl10558_rds)
}

gse72056_file <- file.path(
  external_dir,
  "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
)
if (!file.exists(gse72056_file)) {
  message("Downloading GSE72056 processed single-cell matrix...")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz",
    destfile = gse72056_file,
    mode = "wb",
    quiet = FALSE
  )
}

mcp_genes_file <- file.path(external_dir, "MCPcounter_genes.txt")
if (!file.exists(mcp_genes_file)) {
  message("Downloading the versioned MCP-counter gene signatures...")
  download.file(
    "https://raw.githubusercontent.com/ebecht/MCPcounter/b6eac73/Signatures/genes.txt",
    destfile = mcp_genes_file,
    mode = "wb",
    quiet = FALSE
  )
}

manifest <- data.frame(
  file = c(gse65904_rds, gpl10558_rds, gse72056_file, mcp_genes_file),
  bytes = file.info(c(gse65904_rds, gpl10558_rds, gse72056_file, mcp_genes_file))$size,
  md5 = tools::md5sum(c(gse65904_rds, gpl10558_rds, gse72056_file, mcp_genes_file)),
  row.names = NULL
)
write.csv(manifest, file.path(external_dir, "download_manifest.csv"), row.names = FALSE)
print(manifest)
