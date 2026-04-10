################################################################################
# MICROARRAY PREPROCESSING PIPELINE
# Study: Adrenocortical Adenoma (ACA) vs Adrenocortical Carcinoma (ACC)
# Datasets: GSE12368 | GSE10927
# Author: [Alexander Mwaimbe]
# Date: 2026-04-07
#
# PIPELINE OVERVIEW:
#   We will process each dataset individually:
#     1. Data download (CEL files preferred)
#     2. Raw QC (boxplot, RNA degradation, NUSE, RLE)
#     3. Optional outlier removal
#     4. RMA normalization
#     5. Post-RMA QC (boxplot, PCA, dendrogram)
#     6. Low-expression probe filtering
#     7. Probe → gene symbol annotation
#     8. Probe collapsing to gene-level
#     9. Metadata curation (ACA vs ACC)
#    10. Save outputs
#
# OUTPUT: Independent gene level expression matrices ready for downstream
#         merging, batch correction, and then we perform splitting (70/30%).
#
# ⚠️ DO NOT merge datasets in this script.
# ⚠️ DO NOT perform batch correction here.
# ⚠️ DO NOT perform DEG or ML here.
################################################################################


################################################################################
# SECTION 1: PACKAGE MANAGEMENT
################################################################################

# ---- 1.1 Install BiocManager if needed ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ---- 1.2 CRAN packages ----
cran_pkgs <- c("ggplot2", "data.table", "reshape2", "RColorBrewer",
               "pheatmap", "gridExtra", "dplyr", "stringr", "ragg")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# ---- 1.3 Bioconductor packages ----
bioc_pkgs <- c(
  "GEOquery",    # GEO data retrieval
  "affy",        # CEL file reading & RMA (Affy 1-channel arrays)
  "oligo",       # Alternative for newer Affy/NimbleGen arrays
  "affyPLM",     # NUSE & RLE QC plots
  "limma",       # Probe collapsing (avereps/collapseRows)
  "AnnotationDbi",
  "hgu133plus2.db",   # GSE12368 platform annotation
  "hgu133a.db",       # GSE10927 platform annotation (verify)
  "sva"               # Loaded for downstream compatibility (not used here)
)
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

# ---- 1.4 Load all libraries ----
suppressPackageStartupMessages({
  library(GEOquery)
  library(affy)
  library(oligo)
  library(affyPLM)
  library(limma)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(hgu133a.db)
  library(ggplot2)
  library(data.table)
  library(reshape2)
  library(RColorBrewer)
  library(pheatmap)
  library(dplyr)
  library(stringr)
  library(sva)
})

# ---- 1.5 Session info (for reproducibility) ----
# Saved to file at the end. Print to console now for reference.
cat("=== R Session Info ===\n")
print(sessionInfo())


################################################################################
# SECTION 2: DIRECTORY STRUCTURE
################################################################################
# Here please set your directory, you can modify it to the way you like 

dirs <- c(
  "project/data/GSE12368",
  "project/data/GSE10927",
  "project/results",
  "project/results/qc_plots/GSE12368",
  "project/results/qc_plots/GSE10927",
  "project/results/expression",
  "project/metadata",
  "project/scripts"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
cat("Directory structure created.\n")


################################################################################
# HELPER FUNCTIONS
# (defined once, used for both datasets)
################################################################################

#' Save a base-R plot to PDF, TIFF, and PNG
#'
#' @param expr  An unquoted expression that draws a base-R plot.
#' @param prefix Full file path prefix (no extension).
#' @param width  Width in inches.
#' @param height Height in inches.
save_base_plot <- function(expr, prefix, width = 10, height = 7) {
  # PDF
  pdf(paste0(prefix, ".pdf"), width = width, height = height)
  eval(expr)
  dev.off()
  # TIFF 700 dpi
  tiff(paste0(prefix, ".tiff"), width = width, height = height,
       units = "in", res = 700, compression = "lzw")
  eval(expr)
  dev.off()
  # PNG 700 dpi
  png(paste0(prefix, ".png"), width = width, height = height,
      units = "in", res = 700)
  eval(expr)
  dev.off()
  cat("  Saved:", prefix, "[PDF/TIFF/PNG]\n")
}

#' Save a ggplot object to PDF, TIFF, and PNG
#'
#' @param p      A ggplot object.
#' @param prefix Full file path prefix (no extension).
#' @param width  Width in inches.
#' @param height Height in inches.
save_gg_plot <- function(p, prefix, width = 10, height = 7) {
  ggplot2::ggsave(paste0(prefix, ".pdf"),  plot = p, width = width, height = height, dpi = 700)
  ggplot2::ggsave(paste0(prefix, ".tiff"), plot = p, width = width, height = height, dpi = 700, device = "tiff")
  ggplot2::ggsave(paste0(prefix, ".png"),  plot = p, width = width, height = height, dpi = 700)
  cat("  Saved:", prefix, "[PDF/TIFF/PNG]\n")
}


################################################################################
# ==============================================================================
#  DATASET 1: GSE12368
# ==============================================================================
################################################################################

cat("\n\n========================================\n")
cat("PROCESSING DATASET: GSE12368\n")
cat("========================================\n")

GEO_ID_1   <- "GSE12368"
DATA_DIR_1 <- "project/data/GSE12368"
QC_DIR_1   <- "project/results/qc_plots/GSE12368"

# Platform: GPL570 = Affymetrix Human Genome U133 Plus 2.0 Array
# Annotation DB: hgu133plus2.db

# ===========================================================================
# STEP 1 (GSE12368): DOWNLOAD RAW DATA
# ===========================================================================
# Make sure you only download ACA and ACC samples. Check the dataset on the GEO database to ascertain which samples are ACA and ACC. Do not download normal samples

cat("\n[GSE12368] Step 1: Downloading data...\n")

# GEOquery: attempt to fetch supplementary files (CEL files are in SOFT/suppl)
# getGEOSuppFiles() downloads all supplementary files for the accession.
# CEL.gz files will be extracted below.
if (!file.exists(file.path(DATA_DIR_1, "filelist.txt"))) {
  tryCatch({
    supp_files <- GEOquery::getGEOSuppFiles(GEO_ID_1, makeDirectory = FALSE,
                                             baseDir = DATA_DIR_1)
    cat("[GSE12368] Supplementary files downloaded:\n")
    print(rownames(supp_files))
  }, error = function(e) {
    cat("[GSE12368] WARNING: Could not download supplementary files.\n")
    cat("  Error:", conditionMessage(e), "\n")
    cat("  Falling back to series matrix. See FALLBACK NOTE below.\n")
  })
  # Write a marker file so we skip on re-runs
  writeLines("downloaded", file.path(DATA_DIR_1, "filelist.txt"))
}

# ---- Untar / decompress CEL archives ----
# GEO typically packages all CEL files in a .tar or individual .cel.gz files.
tar_files <- list.files(DATA_DIR_1, pattern = "\\.tar$|\\.TAR$", full.names = TRUE)
if (length(tar_files) > 0) {
  for (tf in tar_files) {
    cat("[GSE12368] Extracting tar:", tf, "\n")
    untar(tf, exdir = DATA_DIR_1)
  }
}
# Decompress .gz files (CEL.gz → CEL)
gz_files <- list.files(DATA_DIR_1, pattern = "\\.cel\\.gz$|\\.CEL\\.gz$",
                        full.names = TRUE, recursive = TRUE)
for (gz in gz_files) {
  out <- sub("\\.gz$", "", gz)
  if (!file.exists(out)) {
    cat("[GSE12368] Decompressing:", basename(gz), "\n")
    R.utils::gunzip(gz, destname = out, overwrite = FALSE, remove = FALSE)
  }
}

# ---- Locate CEL files ----
cel_files_1 <- list.files(DATA_DIR_1, pattern = "\\.cel$|\\.CEL$",
                            full.names = TRUE, recursive = TRUE)

if (length(cel_files_1) == 0) {
  # =========================================================================
  # FALLBACK NOTE (GSE12368):
  #   CEL files were not found. This can happen if:
  #     (a) The GEO submitter did not deposit raw CEL files, OR
  #     (b) The download failed due to connectivity / GEO server issues.
  #   LIMITATION: Without raw CEL files, true RMA cannot be performed.
  #   We fall back to the series matrix (pre-processed, log2 intensities).
  #   Downstream results should be interpreted with caution as batch-correction
  #   and normalization assumptions may differ between datasets.
  # =========================================================================
  cat("[GSE12368] ⚠️  CEL files not found. Falling back to series matrix.\n")

  gse1_raw <- GEOquery::getGEO(GEO_ID_1, destdir = DATA_DIR_1,
                                 GSEMatrix = TRUE, AnnotGPL = FALSE)
  eset1_raw <- gse1_raw[[1]]

  cat("[GSE12368] Series matrix loaded. Dimensions:", dim(exprs(eset1_raw)), "\n")
  cat("[GSE12368] ⚠️  RMA will NOT be performed on series matrix data.\n")
  cat("           RMA-equivalent step is skipped; values are used as-is.\n")
  USE_CEL_1 <- FALSE

} else {
  cat("[GSE12368] Found", length(cel_files_1), "CEL files. Proceeding with raw pipeline.\n")
  print(basename(cel_files_1))
  USE_CEL_1 <- TRUE
}


# ===========================================================================
# STEP 2 (GSE12368): RAW QC — BEFORE NORMALIZATION
# ===========================================================================
# PURPOSE:We Identify low-quality arrays BEFORE any normalization so that
# defective samples will not corrupt the normalization step.
#
# We will generate the following Plots:
#   (a) Raw intensity boxplot   — should have similar median raw intensities
#   (b) RNA degradation plot    — 5'-to-3' signal ratio; degraded RNA → slope
#   (c) NUSE plot               — Normalized Unscaled Standard Error per probe;
#                                  NUSE > 1.1 suggests poor hybridization
#   (d) RLE plot                — Relative Log Expression; ideally centered at 0
# ===========================================================================

cat("\n[GSE12368] Step 2: Raw QC...\n")

if (USE_CEL_1) {

  # Read CEL files with affy (for HG-U133 Plus 2.0 arrays)
  raw_affy_1 <- affy::ReadAffy(filenames = cel_files_1)

  # (a) Raw intensity boxplot
  # Interpretation: Boxes should be roughly aligned before normalization.
  # Severely shifted boxes indicate potential scanning or labeling issues.
  save_base_plot(
    quote({
      affy::boxplot(raw_affy_1,
                    main = "GSE12368 — Raw Intensity Boxplot (Pre-RMA)",
                    col  = RColorBrewer::brewer.pal(min(12, ncol(raw_affy_1)), "Set3"),
                    las  = 2, cex.axis = 0.6,
                    ylab = "Raw Probe Intensity (log2)")
    }),
    prefix = file.path(QC_DIR_1, "GSE12368_raw_boxplot")
  )

  # (b) RNA degradation plot
  # Interpretation: Degraded samples show a steep positive slope (3' bias).
  # Outlier samples with drastically different slope should be flagged.
  deg_1 <- affy::AffyRNAdeg(raw_affy_1)
  save_base_plot(
    quote({
      affy::plotAffyRNAdeg(deg_1)
      title(main = "GSE12368 — RNA Degradation Plot (Pre-RMA)")
      legend("topleft", legend = Biobase::sampleNames(raw_affy_1),
             col = 1:length(Biobase::sampleNames(raw_affy_1)),
             lty = 1, cex = 0.4, bty = "n")
    }),
    prefix = file.path(QC_DIR_1, "GSE12368_rna_degradation"),
    width = 12, height = 7
  )

  # (c) NUSE plot & (d) RLE plot
  # fitPLM fits a probe-level model; required for NUSE and RLE.
  # This step may take several minutes on large datasets.
  cat("[GSE12368] Fitting PLM (required for NUSE/RLE)...\n")
  plm_fit_1 <- affyPLM::fitPLM(raw_affy_1, normalize = FALSE, background = FALSE)

  # NUSE: values clustered near 1.0 = good; > 1.1 = poor hybridization
  save_base_plot(
    quote({
      affyPLM::NUSE(plm_fit_1,
                    main = "GSE12368 — NUSE Plot (Pre-RMA)",
                    las = 2, cex.axis = 0.6,
                    outline = TRUE)
      abline(h = 1.1, col = "red", lty = 2)
      text(0.5, 1.12, "Threshold = 1.1", col = "red", adj = 0, cex = 0.7)
    }),
    prefix = file.path(QC_DIR_1, "GSE12368_NUSE"),
    width = 12, height = 7
  )

  # RLE: values should be centered at 0; large spread or median offset = problem
  save_base_plot(
    quote({
      affyPLM::RLE(plm_fit_1,
                   main = "GSE12368 — RLE Plot (Pre-RMA)",
                   las = 2, cex.axis = 0.6,
                   outline = TRUE)
      abline(h = 0, col = "blue", lty = 2)
    }),
    prefix = file.path(QC_DIR_1, "GSE12368_RLE"),
    width = 12, height = 7
  )

  cat("[GSE12368] Raw QC plots saved to:", QC_DIR_1, "\n")

} else {
  cat("[GSE12368] Skipping raw QC (no CEL files; series matrix fallback active).\n")
}


# ===========================================================================
# STEP 3 (GSE12368): OPTIONAL BUT RUN IT — REMOVE POOR-QUALITY SAMPLES
# ===========================================================================
# ⚠️ USER ACTION REQUIRED:
#   After reviewing the QC plots above, identify any samples that:
#     - Have NUSE median > 1.1
#     - Have extreme RLE spread (IQR >> 0)
#     - Have RNA degradation slope dramatically different from others
#     - Show outlier intensity distributions in the boxplot
#
#   Enter those sample names below and set REMOVE_SAMPLES_1 accordingly.
# ===========================================================================

# --- USER: Edit this vector with bad sample names (CEL file basenames) ---
BAD_SAMPLES_1 <- c("GSM310462.CEL", "GSM310464.CEL", "GSM310478.CEL", "GSM310492.CEL","GSM310490.CEL", "GSM310482.CEL","GSM310484.CEL", "GSM310483.CEL",
                   "GSM310487.CEL","GSM310488.CEL","GSM310489.CEL","GSM310491.CEL")
# Example:
# BAD_SAMPLES_1 <- c("GSM310123.CEL", "GSM310199.CEL")

if (USE_CEL_1 && length(BAD_SAMPLES_1) > 0) {
  keep_idx_1  <- !basename(cel_files_1) %in% BAD_SAMPLES_1
  cel_files_1 <- cel_files_1[keep_idx_1]
  cat("[GSE12368] Removed", sum(!keep_idx_1), "poor-quality samples.\n")
  cat("[GSE12368] Remaining samples:", length(cel_files_1), "\n")
  # Re-read after removal
  raw_affy_1 <- affy::ReadAffy(filenames = cel_files_1)
} else {
  cat("[GSE12368] No samples removed (BAD_SAMPLES_1 is empty).\n")
}


# ===========================================================================
# STEP 4 (GSE12368): RMA NORMALIZATION
# ===========================================================================
# WHY RMA?
#   RMA (Robust Multi-array Average) is the gold standard for Affymetrix arrays.
#   It performs three operations:
#     1. Background correction  — removes optical noise using convolution model
#     2. Quantile normalization — forces identical empirical distributions across
#                                 arrays, removing array-to-array technical variation
#     3. Log2 transformation    — stabilizes variance and reduces skew
#     4. Probe summarization    — aggregates multiple probes per probe set (median polish)
#
# WHY BEFORE ANNOTATION?
#   RMA operates on raw probe-level data. Annotation is a post-hoc mapping
#   from probe set IDs → gene symbols, applied after the expression matrix
#   is in final form.
#
# ⚠️ RMA is performed ONLY on CEL files. Series matrix data is already
#    pre-processed and cannot be re-normalized with RMA.
# ===========================================================================

cat("\n[GSE12368] Step 4: RMA Normalization...\n")

if (USE_CEL_1) {
  eset1 <- affy::rma(raw_affy_1)
  cat("[GSE12368] RMA complete. Expression matrix dimensions:", dim(Biobase::exprs(eset1)), "\n")
} else {
  # Series matrix fallback: use pre-processed values as-is (already log2)
  eset1 <- eset1_raw
  cat("[GSE12368] Using pre-processed series matrix expression values.\n")
}


# ===========================================================================
# STEP 5 (GSE12368): POST-RMA QC
# ===========================================================================
# EXPECTED OUTCOMES after RMA:
#   - Boxplot: All arrays should have nearly identical distributions (flat boxes)
#   - PCA:     Samples should cluster loosely by biological group (ACA vs ACC),
#              not by technical factors. Outlier clusters may indicate batch effects.
#   - Dendrogram: Biological replicates within the same group should cluster together.
#
# DETECTING REMAINING OUTLIERS:
#   - An isolated sample in PCA far from all others = candidate for removal
#   - A sample that clusters with the wrong group in the dendrogram = verify metadata
# ===========================================================================

cat("\n[GSE12368] Step 5: Post-RMA QC...\n")

expr_mat_1 <- Biobase::exprs(eset1)
n_samples_1 <- ncol(expr_mat_1)

# (a) Post-RMA boxplot
p_box_post_1 <- {
  df_long <- reshape2::melt(expr_mat_1, varnames = c("probe", "sample"), value.name = "expr")
  ggplot(df_long, aes(x = sample, y = expr, fill = sample)) +
    geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
    theme_bw(base_size = 9) +
    labs(title = "GSE12368 — Post-RMA Normalized Intensity",
         x = "Sample", y = "log2 Intensity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
}
save_gg_plot(p_box_post_1,
             prefix = file.path(QC_DIR_1, "GSE12368_postRMA_boxplot"),
             width = max(10, n_samples_1 * 0.35), height = 6)

# (b) PCA plot
pca_res_1 <- prcomp(t(expr_mat_1), scale. = TRUE)
pca_df_1  <- data.frame(
  sample = rownames(pca_res_1$x),
  PC1    = pca_res_1$x[, 1],
  PC2    = pca_res_1$x[, 2]
)
var_explained_1 <- round(100 * pca_res_1$sdev^2 / sum(pca_res_1$sdev^2), 1)

p_pca_1 <- ggplot(pca_df_1, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3, colour = "#2196F3", alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
  theme_bw(base_size = 10) +
  labs(
    title = "GSE12368 — PCA (Post-RMA)",
    subtitle = "Colour by group after metadata curation in Step 9",
    x = paste0("PC1 (", var_explained_1[1], "% variance)"),
    y = paste0("PC2 (", var_explained_1[2], "% variance)")
  )

# Install ggrepel if needed
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

save_gg_plot(p_pca_1,
             prefix = file.path(QC_DIR_1, "GSE12368_postRMA_PCA"),
             width = 9, height = 7)

# (c) Hierarchical clustering dendrogram
dist_mat_1  <- dist(t(expr_mat_1), method = "euclidean")
hclust_1    <- hclust(dist_mat_1, method = "complete")

save_base_plot(
  quote({
    plot(hclust_1,
         main = "GSE12368 — Hierarchical Clustering (Post-RMA)",
         xlab = "", sub = "",
         cex = 0.6, hang = -1)
  }),
  prefix = file.path(QC_DIR_1, "GSE12368_postRMA_dendrogram"),
  width = max(12, n_samples_1 * 0.4), height = 7
)

cat("[GSE12368] Post-RMA QC plots saved.\n")


# ===========================================================================
# STEP 6 (GSE12368): PROBE FILTERING
# ===========================================================================
# Rationale:
#   Low-expression probes add noise and reduce statistical power.
#   Retaining ~50-75% of probes is typical for Affymetrix arrays.
#
# Thresholds used:
#   - Expression threshold : log2 intensity > 3.0 in >= 20% of samples
#     (3.0 ~ 8-fold above background; commonly used in literature)
#   - Variance threshold   : configurable (set to FALSE by default)
#
# Impact on WGCNA:
#   WGCNA is sensitive to low-variance genes (they contribute noise to the
#   adjacency matrix). Variance filtering before WGCNA is strongly recommended.
#
# Impact on ML:
#   High-dimensional matrices with uninformative features slow training and
#   may reduce classifier accuracy. Pre-filtering reduces feature space.
# ===========================================================================

cat("\n[GSE12368] Step 6: Probe Filtering...\n")

# --- Configurable thresholds ---
EXPR_THRESHOLD_1    <- 3.0    # log2 intensity cutoff
SAMPLE_FRACTION_1   <- 0.20   # must exceed threshold in this fraction of samples
APPLY_VAR_FILTER_1  <- FALSE  # set TRUE to also apply variance filter
VAR_QUANTILE_1      <- 0.25   # if APPLY_VAR_FILTER_1=TRUE, remove bottom quartile

# Low-expression filter
n_above_1    <- rowSums(expr_mat_1 > EXPR_THRESHOLD_1)
keep_expr_1  <- n_above_1 >= (SAMPLE_FRACTION_1 * n_samples_1)
expr_filt_1  <- expr_mat_1[keep_expr_1, ]
cat("[GSE12368] Probes before expression filter:", nrow(expr_mat_1), "\n")
cat("[GSE12368] Probes after  expression filter:", nrow(expr_filt_1), "\n")

# Optional variance filter
if (APPLY_VAR_FILTER_1) {
  probe_vars_1   <- apply(expr_filt_1, 1, var)
  var_cutoff_1   <- quantile(probe_vars_1, VAR_QUANTILE_1)
  expr_filt_1    <- expr_filt_1[probe_vars_1 >= var_cutoff_1, ]
  cat("[GSE12368] Probes after variance filter:", nrow(expr_filt_1), "\n")
} else {
  cat("[GSE12368] Variance filter not applied (APPLY_VAR_FILTER_1 = FALSE).\n")
}


# ===========================================================================
# STEP 7 (GSE12368): PROBE ANNOTATION
# ===========================================================================
# Platform: GPL570 → hgu133plus2.db
#
# Mapping: probe set ID → SYMBOL (HGNC gene symbol)
#   - Probes without a valid SYMBOL are removed (unannotated controls, obsolete)
#   - One probe set can map to multiple genes (removed)
#   - Multiple probe sets can map to the same gene (collapsed in Step 8)
#
# If using series matrix fallback:
#   fData(eset1) may already contain gene symbols. We prefer AnnotationDbi
#   for consistency and reproducibility.
# ===========================================================================

cat("\n[GSE12368] Step 7: Probe Annotation (GPL570 → hgu133plus2.db)...\n")

probe_ids_1 <- rownames(expr_filt_1)

# Map probe IDs to gene symbols via AnnotationDbi
anno_1 <- AnnotationDbi::select(
  hgu133plus2.db,
  keys    = probe_ids_1,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)

# Keep only 1-to-1 mappings (remove probes mapping to multiple genes)
# Explanation: Multi-gene probes are ambiguous and should not be attributed
# to any single gene.
anno_1 <- anno_1[!is.na(anno_1$SYMBOL), ]
dup_probes_1 <- anno_1$PROBEID[duplicated(anno_1$PROBEID)]
anno_1 <- anno_1[!anno_1$PROBEID %in% dup_probes_1, ]

cat("[GSE12368] Annotated probe sets:", nrow(anno_1), "\n")
cat("[GSE12368] Unique gene symbols:", length(unique(anno_1$SYMBOL)), "\n")

# Subset expression matrix to annotated probes only
expr_anno_1 <- expr_filt_1[rownames(expr_filt_1) %in% anno_1$PROBEID, ]

# Attach gene symbol as rowname-compatible column
anno_1 <- anno_1[match(rownames(expr_anno_1), anno_1$PROBEID), ]
rownames(anno_1) <- anno_1$PROBEID


# ===========================================================================
# STEP 8 (GSE12368): COLLAPSE PROBES TO GENE-LEVEL
# ===========================================================================
# When multiple probe sets map to the same gene symbol, we must select ONE
# representative value per gene.
#
# Available methods:
#   "mean"    — average all probe sets; default, robust, reduces noise
#   "median"  — less sensitive to outlier probe sets
#   "maxvar"  — keep the probe set with highest variance; preferred for WGCNA
#                because high-variance probes capture more biological signal
#
# Recommendation:
#   - For WGCNA  → use "maxvar" (preserves biological variation)
#   - For ML     → use "mean" or "maxvar" (both acceptable)
#   - "median"   → conservative choice, less common in literature
#
# We use "mean" as default; change COLLAPSE_METHOD_1 to suit your analysis.
# ===========================================================================

cat("\n[GSE12368] Step 8: Collapsing probes to gene-level...\n")

COLLAPSE_METHOD_1 <- "mean"  # options: "mean", "median", "maxvar"

# Add gene symbol as a column for limma::avereps or manual collapsing
expr_with_sym_1 <- data.frame(
  SYMBOL = anno_1[rownames(expr_anno_1), "SYMBOL"],
  expr_anno_1,
  check.names = FALSE
)

if (COLLAPSE_METHOD_1 == "mean") {
  # limma::avereps: computes mean of duplicate rows
  gene_expr_1 <- limma::avereps(expr_anno_1,
                                  ID = anno_1[rownames(expr_anno_1), "SYMBOL"])

} else if (COLLAPSE_METHOD_1 == "median") {
  # Manual median collapsing
  gene_expr_1 <- do.call(rbind, lapply(
    split(as.data.frame(expr_anno_1),
          anno_1[rownames(expr_anno_1), "SYMBOL"]),
    function(g) apply(g, 2, median)
  ))
  gene_expr_1 <- as.matrix(gene_expr_1)

} else if (COLLAPSE_METHOD_1 == "maxvar") {
  # Keep probe set with highest row variance per gene
  probe_var_1 <- apply(expr_anno_1, 1, var)
  sym_1       <- anno_1[rownames(expr_anno_1), "SYMBOL"]
  # For each gene, find the probe with max variance
  best_probe_1 <- tapply(seq_along(sym_1), sym_1, function(idx) {
    idx[which.max(probe_var_1[idx])]
  })
  gene_expr_1  <- expr_anno_1[unlist(best_probe_1), ]
  rownames(gene_expr_1) <- names(best_probe_1)
}

cat("[GSE12368] Gene-level matrix dimensions:", dim(gene_expr_1), "\n")
cat("[GSE12368] Genes:", nrow(gene_expr_1), " | Samples:", ncol(gene_expr_1), "\n")


# ===========================================================================
# STEP 9 (GSE12368): METADATA CURATION
# ===========================================================================
# ⚠️ CRITICAL — USER ACTION REQUIRED
#
# PURPOSE: Correctly label each sample as ACA or ACC.
# DO NOT use naive keyword classification (e.g., grepping "ACA" from title).
# Many GEO submissions use inconsistent or abbreviated naming.
#
# PROCEDURE:
#   1. Print the full pheno data table (below).
#   2. Manually inspect the 'title' and 'characteristics_ch1' columns.
#   3. Fill in the sample_annotation data.frame with correct labels.
#   4. Remove any normal adrenal cortex samples (not ACA or ACC).
#   5. Remove ambiguous samples (e.g., unclear diagnosis, mixed tumor).
#
# Labels to use:
#   "ACA" — adrenocortical adenoma (benign)
#   "ACC" — adrenocortical carcinoma (malignant)
#   NA    — to be excluded
# ===========================================================================

cat("\n[GSE12368] Step 9: Metadata Curation...\n")

# Fetch GEO metadata (phenotype data)
if (USE_CEL_1) {
  gse1_meta <- GEOquery::getGEO(GEO_ID_1, GSEMatrix = TRUE, AnnotGPL = FALSE,
                                  destdir = DATA_DIR_1)
  pheno_1 <- pData(gse1_meta[[1]])
} else {
  pheno_1 <- pData(eset1_raw)
}

cat("[GSE12368] Pheno data columns available:\n")
print(colnames(pheno_1))

cat("\n[GSE12368] Sample titles and characteristics:\n")
print(pheno_1[, c("title", grep("characteristics", colnames(pheno_1), value = TRUE))])

# ---------------------------------------------------------------------------
# ⚠️ USER: Fill in the manual annotation below.
# Use the printed table above to verify each sample.
# The sample_id values MUST match the column names of gene_expr_1.
# ---------------------------------------------------------------------------

# Example template — REPLACE WITH ACTUAL LABELS AFTER REVIEWING PHENO DATA
# !! This is a PLACEHOLDER. Do NOT use without manual verification !!
sample_annotation_1 <- data.frame(
  sample_id = colnames(gene_expr_1),
  title     = pheno_1[match(colnames(gene_expr_1),
                             rownames(pheno_1)), "title"],
  group     = NA_character_,   # <-- USER FILLS THIS IN
  dataset   = "GSE12368",
  stringsAsFactors = FALSE
)

# Assign ACC (Carcinoma)
acc_cel <- paste0("GSM3104", 75:86, ".CEL")
sample_annotation_1$group[sample_annotation_1$sample_id %in% acc_cel] <- "ACC"

# Assign ACA (Adenoma)
aca_cel <- paste0("GSM3104", 59:74, ".CEL")
sample_annotation_1$group[sample_annotation_1$sample_id %in% aca_cel] <- "ACA"

# Verify: You should see numbers for ACA and ACC now
table(sample_annotation_1$group, useNA = "always")

# ---------------------------------------------------------------------------
# USER EXAMPLE — uncomment and edit with real sample IDs and groups:
# ---------------------------------------------------------------------------
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310459"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310460"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310461"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310462"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310463"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310464"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310465"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310466"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310467"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310468"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310469"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310470"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310471"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310472"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310473"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310474"] <- "ACA"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310475"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310476"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310477"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310478"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310479"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310480"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310481"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310482"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310483"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310484"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310485"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310486"] <- "ACC"
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310487"] <- NA
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310488"] <- NA
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310489"] <- NA
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310490"] <- NA
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310491"] <- NA
sample_annotation_1$group[sample_annotation_1$sample_id == "GSM310492"] <- NA
# ...etc.
# ---------------------------------------------------------------------------

cat("\n[GSE12368] ⚠️  MANUAL ACTION REQUIRED:\n")
cat("  Open sample_annotation_1 in R and fill in the 'group' column.\n")
cat("  Use: ACA, ACC, or NA (for normal/ambiguous — will be excluded).\n")

# Remove samples with NA group (normal, ambiguous, unclassified)
keep_samples_1       <- !is.na(sample_annotation_1$group)
sample_annotation_1  <- sample_annotation_1[keep_samples_1, ]
gene_expr_final_1    <- gene_expr_1[, sample_annotation_1$sample_id, drop = FALSE]

cat("[GSE12368] Samples retained (ACA + ACC):", nrow(sample_annotation_1), "\n")
cat("[GSE12368] Group distribution:\n")
print(table(sample_annotation_1$group))


# ===========================================================================
# STEP 10 (GSE12368): SAVE OUTPUTS
# ===========================================================================

cat("\n[GSE12368] Step 10: Saving outputs...\n")

# Expression matrix (gene-level)
fwrite(
  data.table(gene = rownames(gene_expr_final_1), as.data.frame(gene_expr_final_1)),
  file = "project/results/expression/GSE12368_expr_gene_level.csv"
)
cat("[GSE12368] Expression matrix saved: GSE12368_expr_gene_level.csv\n")

# Metadata
fwrite(
  as.data.table(sample_annotation_1),
  file = "project/metadata/GSE12368_metadata.csv"
)
cat("[GSE12368] Metadata saved: GSE12368_metadata.csv\n")

# Annotation table (probe → gene mapping)
fwrite(
  as.data.table(anno_1),
  file = "project/metadata/GSE12368_probe_annotation.csv"
)

cat("[GSE12368] All outputs saved.\n")
cat("========================================\n")
cat("GSE12368 PREPROCESSING COMPLETE\n")
cat("========================================\n\n")


################################################################################
# ==============================================================================
#  DATASET 2: GSE10927
# ==============================================================================
################################################################################

cat("\n\n========================================\n")
cat("PROCESSING DATASET: GSE10927\n")
cat("========================================\n")

GEO_ID_2   <- "GSE10927"
DATA_DIR_2 <- "project/data/GSE10927"
QC_DIR_2   <- "project/results/qc_plots/GSE10927"

# Platform: VERIFY — likely GPL96 (Affymetrix HG-U133A) or GPL570
# ⚠️ USER: Confirm the platform for GSE10927 on GEO before running.
# If platform is HG-U133A, use: hgu133a.db
# If platform is HG-U133 Plus 2.0, use: hgu133plus2.db

PLATFORM_DB_2 <- hgu133a.db   # change if needed based on GEO platform page


# ===========================================================================
# STEP 1 (GSE10927): DOWNLOAD RAW DATA
# ===========================================================================

cat("\n[GSE10927] Step 1: Downloading data...\n")

if (!file.exists(file.path(DATA_DIR_2, "filelist.txt"))) {
  tryCatch({
    supp_files_2 <- GEOquery::getGEOSuppFiles(GEO_ID_2, makeDirectory = FALSE,
                                               baseDir = DATA_DIR_2)
    cat("[GSE10927] Supplementary files downloaded:\n")
    print(rownames(supp_files_2))
  }, error = function(e) {
    cat("[GSE10927] WARNING: Could not download supplementary files.\n")
    cat("  Error:", conditionMessage(e), "\n")
    cat("  Falling back to series matrix.\n")
  })
  writeLines("downloaded", file.path(DATA_DIR_2, "filelist.txt"))
}

# Untar / decompress
tar_files_2 <- list.files(DATA_DIR_2, pattern = "\\.tar$|\\.TAR$", full.names = TRUE)
if (length(tar_files_2) > 0) {
  for (tf in tar_files_2) {
    cat("[GSE10927] Extracting tar:", tf, "\n")
    untar(tf, exdir = DATA_DIR_2)
  }
}
gz_files_2 <- list.files(DATA_DIR_2, pattern = "\\.cel\\.gz$|\\.CEL\\.gz$",
                           full.names = TRUE, recursive = TRUE)
for (gz in gz_files_2) {
  out <- sub("\\.gz$", "", gz)
  if (!file.exists(out)) {
    cat("[GSE10927] Decompressing:", basename(gz), "\n")
    R.utils::gunzip(gz, destname = out, overwrite = FALSE, remove = FALSE)
  }
}

cel_files_2 <- list.files(DATA_DIR_2, pattern = "\\.cel$|\\.CEL$",
                            full.names = TRUE, recursive = TRUE)

if (length(cel_files_2) == 0) {
  # =========================================================================
  # FALLBACK NOTE (GSE10927): Same as GSE12368. See above.
  # =========================================================================
  cat("[GSE10927] ⚠️  CEL files not found. Falling back to series matrix.\n")

  gse2_raw <- GEOquery::getGEO(GEO_ID_2, destdir = DATA_DIR_2,
                                 GSEMatrix = TRUE, AnnotGPL = FALSE)
  eset2_raw <- gse2_raw[[1]]
  cat("[GSE10927] Series matrix loaded. Dimensions:", dim(exprs(eset2_raw)), "\n")
  USE_CEL_2 <- FALSE

} else {
  cat("[GSE10927] Found", length(cel_files_2), "CEL files.\n")
  print(basename(cel_files_2))
  USE_CEL_2 <- TRUE
}


# ===========================================================================
# STEP 2 (GSE10927): RAW QC
# ===========================================================================

cat("\n[GSE10927] Step 2: Raw QC...\n")

if (USE_CEL_2) {

  raw_affy_2 <- affy::ReadAffy(filenames = cel_files_2)

  # (a) Raw boxplot
  save_base_plot(
    quote({
      affy::boxplot(raw_affy_2,
                    main = "GSE10927 — Raw Intensity Boxplot (Pre-RMA)",
                    col  = RColorBrewer::brewer.pal(min(12, ncol(raw_affy_2)), "Paired"),
                    las  = 2, cex.axis = 0.6,
                    ylab = "Raw Probe Intensity (log2)")
    }),
    prefix = file.path(QC_DIR_2, "GSE10927_raw_boxplot")
  )

  # (b) RNA degradation
  deg_2 <- affy::AffyRNAdeg(raw_affy_2)
  save_base_plot(
    quote({
      affy::plotAffyRNAdeg(deg_2)
                            
      legend("topleft", legend = sampleNames(raw_affy_2),
             col = 1:length(sampleNames(raw_affy_2)),
             lty = 1, cex = 0.4, bty = "n")
    }),
    prefix = file.path(QC_DIR_2, "GSE10927_rna_degradation"),
    width = 12, height = 7
  )

  # (c) NUSE & (d) RLE
  cat("[GSE10927] Fitting PLM...\n")
  plm_fit_2 <- affyPLM::fitPLM(raw_affy_2, normalize = FALSE, background = FALSE)

  save_base_plot(
    quote({
      affyPLM::NUSE(plm_fit_2,
                    main = "GSE10927 — NUSE Plot (Pre-RMA)",
                    las = 2, cex.axis = 0.6, outline = TRUE)
      abline(h = 1.1, col = "red", lty = 2)
      text(0.5, 1.12, "Threshold = 1.1", col = "red", adj = 0, cex = 0.7)
    }),
    prefix = file.path(QC_DIR_2, "GSE10927_NUSE"),
    width = 12, height = 7
  )

  save_base_plot(
    quote({
      affyPLM::RLE(plm_fit_2,
                   main = "GSE10927 — RLE Plot (Pre-RMA)",
                   las = 2, cex.axis = 0.6, outline = TRUE)
      abline(h = 0, col = "blue", lty = 2)
    }),
    prefix = file.path(QC_DIR_2, "GSE10927_RLE"),
    width = 12, height = 7
  )

  cat("[GSE10927] Raw QC plots saved.\n")

} else {
  cat("[GSE10927] Skipping raw QC (series matrix fallback).\n")
}


# ===========================================================================
# STEP 3 (GSE10927): OPTIONAL — REMOVE POOR-QUALITY SAMPLES
# ===========================================================================
# ⚠️ USER ACTION REQUIRED — same instructions as GSE12368 Step 3

BAD_SAMPLES_2 <- c("GSM277101.CEL","GSM277104.CEL","GSM277105.CEL","GSM277107.CEL","GSM277119.CEL",
                   "GSM277121.CEL","GSM277122.CEL","GSM277123.CEL","GSM277126.CEL","GSM277128.CEL",
                   "GSM277127.CEL","GSM277129.CEL","GSM277132.CEL","GSM277133.CEL","GSM277134.CEL",
                   "GSM277135.CEL","GSM277137.CEL","GSM277143.CEL","GSM277145.CEL","GSM277144.CEL",
                   "GSM277154.CEL")
# Example:
# BAD_SAMPLES_2 <- c("GSM274001.CEL")

if (USE_CEL_2 && length(BAD_SAMPLES_2) > 0) {
  keep_idx_2  <- !basename(cel_files_2) %in% BAD_SAMPLES_2
  cel_files_2 <- cel_files_2[keep_idx_2]
  cat("[GSE10927] Removed", sum(!keep_idx_2), "poor-quality samples.\n")
  raw_affy_2  <- affy::ReadAffy(filenames = cel_files_2)
} else {
  cat("[GSE10927] No samples removed.\n")
}


# ===========================================================================
# STEP 4 (GSE10927): RMA NORMALIZATION
# ===========================================================================

cat("\n[GSE10927] Step 4: RMA Normalization...\n")

if (USE_CEL_2) {
  eset2 <- affy::rma(raw_affy_2)
  cat("[GSE10927] RMA complete. Dimensions:", dim(exprs(eset2)), "\n")
} else {
  eset2 <- eset2_raw
  cat("[GSE10927] Using pre-processed series matrix expression values.\n")
}


# ===========================================================================
# STEP 5 (GSE10927): POST-RMA QC
# ===========================================================================

cat("\n[GSE10927] Step 5: Post-RMA QC...\n")

expr_mat_2 <- exprs(eset2)
n_samples_2 <- ncol(expr_mat_2)

# (a) Boxplot
p_box_post_2 <- {
  df_long_2 <- reshape2::melt(expr_mat_2, varnames = c("probe", "sample"), value.name = "expr")
  ggplot(df_long_2, aes(x = sample, y = expr, fill = sample)) +
    geom_boxplot(outlier.size = 0.3, show.legend = FALSE) +
    theme_bw(base_size = 9) +
    labs(title = "GSE10927 — Post-RMA Normalized Intensity",
         x = "Sample", y = "log2 Intensity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
}
save_gg_plot(p_box_post_2,
             prefix = file.path(QC_DIR_2, "GSE10927_postRMA_boxplot"),
             width = max(10, n_samples_2 * 0.35), height = 6)

# (b) PCA
pca_res_2 <- prcomp(t(expr_mat_2), scale. = TRUE)
pca_df_2  <- data.frame(
  sample = rownames(pca_res_2$x),
  PC1    = pca_res_2$x[, 1],
  PC2    = pca_res_2$x[, 2]
)
var_explained_2 <- round(100 * pca_res_2$sdev^2 / sum(pca_res_2$sdev^2), 1)

p_pca_2 <- ggplot(pca_df_2, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3, colour = "#E91E63", alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
  theme_bw(base_size = 10) +
  labs(
    title = "GSE10927 — PCA (Post-RMA)",
    subtitle = "Colour by group after metadata curation in Step 9",
    x = paste0("PC1 (", var_explained_2[1], "% variance)"),
    y = paste0("PC2 (", var_explained_2[2], "% variance)")
  )
save_gg_plot(p_pca_2,
             prefix = file.path(QC_DIR_2, "GSE10927_postRMA_PCA"),
             width = 9, height = 7)

# (c) Dendrogram
dist_mat_2 <- dist(t(expr_mat_2), method = "euclidean")
hclust_2   <- hclust(dist_mat_2, method = "complete")

save_base_plot(
  quote({
    plot(hclust_2,
         main = "GSE10927 — Hierarchical Clustering (Post-RMA)",
         xlab = "", sub = "",
         cex = 0.6, hang = -1)
  }),
  prefix = file.path(QC_DIR_2, "GSE10927_postRMA_dendrogram"),
  width = max(12, n_samples_2 * 0.4), height = 7
)

cat("[GSE10927] Post-RMA QC plots saved.\n")


# ===========================================================================
# STEP 6 (GSE10927): PROBE FILTERING
# ===========================================================================

cat("\n[GSE10927] Step 6: Probe Filtering...\n")

EXPR_THRESHOLD_2   <- 3.0
SAMPLE_FRACTION_2  <- 0.20
APPLY_VAR_FILTER_2 <- FALSE
VAR_QUANTILE_2     <- 0.25

n_above_2    <- rowSums(expr_mat_2 > EXPR_THRESHOLD_2)
keep_expr_2  <- n_above_2 >= (SAMPLE_FRACTION_2 * n_samples_2)
expr_filt_2  <- expr_mat_2[keep_expr_2, ]
cat("[GSE10927] Probes before filter:", nrow(expr_mat_2), "\n")
cat("[GSE10927] Probes after  filter:", nrow(expr_filt_2), "\n")

if (APPLY_VAR_FILTER_2) {
  probe_vars_2 <- apply(expr_filt_2, 1, var)
  var_cutoff_2 <- quantile(probe_vars_2, VAR_QUANTILE_2)
  expr_filt_2  <- expr_filt_2[probe_vars_2 >= var_cutoff_2, ]
  cat("[GSE10927] Probes after variance filter:", nrow(expr_filt_2), "\n")
}


# ===========================================================================
# STEP 7 (GSE10927): PROBE ANNOTATION
# ===========================================================================
# ⚠️ PLATFORM NOTE:
#   If GSE10927 uses GPL96 (HG-U133A), probe IDs differ from GPL570.
#   Confirm on the GEO platform page and adjust PLATFORM_DB_2 at top of section.

cat("\n[GSE10927] Step 7: Probe Annotation...\n")

probe_ids_2 <- rownames(expr_filt_2)

anno_2 <- AnnotationDbi::select(
  PLATFORM_DB_2,
  keys    = probe_ids_2,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)

anno_2        <- anno_2[!is.na(anno_2$SYMBOL), ]
dup_probes_2  <- anno_2$PROBEID[duplicated(anno_2$PROBEID)]
anno_2        <- anno_2[!anno_2$PROBEID %in% dup_probes_2, ]

cat("[GSE10927] Annotated probe sets:", nrow(anno_2), "\n")
cat("[GSE10927] Unique gene symbols:", length(unique(anno_2$SYMBOL)), "\n")

expr_anno_2 <- expr_filt_2[rownames(expr_filt_2) %in% anno_2$PROBEID, ]
anno_2      <- anno_2[match(rownames(expr_anno_2), anno_2$PROBEID), ]
rownames(anno_2) <- anno_2$PROBEID


# ===========================================================================
# STEP 8 (GSE10927): COLLAPSE PROBES TO GENE-LEVEL
# ===========================================================================

cat("\n[GSE10927] Step 8: Collapsing probes to gene-level...\n")

COLLAPSE_METHOD_2 <- "mean"   # options: "mean", "median", "maxvar"

if (COLLAPSE_METHOD_2 == "mean") {
  gene_expr_2 <- limma::avereps(expr_anno_2,
                                  ID = anno_2[rownames(expr_anno_2), "SYMBOL"])

} else if (COLLAPSE_METHOD_2 == "median") {
  gene_expr_2 <- do.call(rbind, lapply(
    split(as.data.frame(expr_anno_2),
          anno_2[rownames(expr_anno_2), "SYMBOL"]),
    function(g) apply(g, 2, median)
  ))
  gene_expr_2 <- as.matrix(gene_expr_2)

} else if (COLLAPSE_METHOD_2 == "maxvar") {
  probe_var_2 <- apply(expr_anno_2, 1, var)
  sym_2       <- anno_2[rownames(expr_anno_2), "SYMBOL"]
  best_probe_2 <- tapply(seq_along(sym_2), sym_2, function(idx) {
    idx[which.max(probe_var_2[idx])]
  })
  gene_expr_2  <- expr_anno_2[unlist(best_probe_2), ]
  rownames(gene_expr_2) <- names(best_probe_2)
}

cat("[GSE10927] Gene-level matrix dimensions:", dim(gene_expr_2), "\n")


# ===========================================================================
# STEP 9 (GSE10927): METADATA CURATION
# ===========================================================================
# ⚠️ USER ACTION REQUIRED — same instructions as GSE12368 Step 9
# ===========================================================================

cat("\n[GSE10927] Step 9: Metadata Curation...\n")

if (USE_CEL_2) {
  gse2_meta <- GEOquery::getGEO(GEO_ID_2, GSEMatrix = TRUE, AnnotGPL = FALSE,
                                  destdir = DATA_DIR_2)
  pheno_2 <- pData(gse2_meta[[1]])
} else {
  pheno_2 <- pData(eset2_raw)
}

cat("[GSE10927] Pheno data columns available:\n")
print(colnames(pheno_2))

cat("\n[GSE10927] Sample titles and characteristics:\n")
print(pheno_2[, c("title", grep("characteristics", colnames(pheno_2), value = TRUE))])

sample_annotation_2 <- data.frame(
  sample_id = colnames(gene_expr_2),
  title     = pheno_2[match(colnames(gene_expr_2),
                             rownames(pheno_2)), "title"],
  group     = NA_character_,
  dataset   = "GSE10927",
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# USER EXAMPLE — uncomment and edit:
# ---------------------------------------------------------------------------
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277090.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277091.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277092.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277093.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277094.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277095.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277096.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277097.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277098.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277099.CEL"] <- NA
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277100.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277101.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277102.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277103.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277104.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277105.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277106.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277107.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277108.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277109.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277110.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277111.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277112.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277113.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277114.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277115.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277116.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277117.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277118.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277119.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277120.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277121.CEL"] <- "ACA"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277122.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277123.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277124.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277125.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277126.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277127.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277128.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277129.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277130.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277131.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277132.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277133.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277134.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277135.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277136.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277137.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277138.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277139.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277140.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277141.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277142.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277143.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277144.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277145.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277146.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277147.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277148.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277149.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277150.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277151.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277152.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277153.CEL"] <- "ACC"
sample_annotation_2$group[sample_annotation_2$sample_id == "GSM277154.CEL"] <- "ACC"
# ...etc.
# ---------------------------------------------------------------------------

cat("\n[GSE10927] ⚠️  MANUAL ACTION REQUIRED:\n")
cat("  Fill in the 'group' column of sample_annotation_2.\n")
cat("  Use: ACA, ACC, or NA (for normal/ambiguous — will be excluded).\n")

keep_samples_2      <- !is.na(sample_annotation_2$group)
sample_annotation_2 <- sample_annotation_2[keep_samples_2, ]
gene_expr_final_2   <- gene_expr_2[, sample_annotation_2$sample_id, drop = FALSE]

cat("[GSE10927] Samples retained:", nrow(sample_annotation_2), "\n")
cat("[GSE10927] Group distribution:\n")
print(table(sample_annotation_2$group))


# ===========================================================================
# STEP 10 (GSE10927): SAVE OUTPUTS
# ===========================================================================

cat("\n[GSE10927] Step 10: Saving outputs...\n")

fwrite(
  data.table(gene = rownames(gene_expr_final_2), as.data.frame(gene_expr_final_2)),
  file = "project/results/expression/GSE10927_expr_gene_level.csv"
)
cat("[GSE10927] Expression matrix saved: GSE10927_expr_gene_level.csv\n")

fwrite(
  as.data.table(sample_annotation_2),
  file = "project/metadata/GSE10927_metadata.csv"
)
cat("[GSE10927] Metadata saved: GSE10927_metadata.csv\n")

fwrite(
  as.data.table(anno_2),
  file = "project/metadata/GSE10927_probe_annotation.csv"
)

cat("[GSE10927] All outputs saved.\n")
cat("========================================\n")
cat("GSE10927 PREPROCESSING COMPLETE\n")
cat("========================================\n")


################################################################################
# FINAL: SESSION INFO + PIPELINE SUMMARY
################################################################################

cat("\n\n============================================================\n")
cat("PIPELINE COMPLETE — BOTH DATASETS PREPROCESSED INDEPENDENTLY\n")
cat("============================================================\n\n")

cat("Output files:\n")
cat("  project/results/expression/GSE12368_expr_gene_level.csv\n")
cat("  project/results/expression/GSE10927_expr_gene_level.csv\n")
cat("  project/metadata/GSE12368_metadata.csv\n")
cat("  project/metadata/GSE10927_metadata.csv\n")
cat("  project/results/qc_plots/GSE12368/  [PDF + TIFF 700dpi + PNG 700dpi]\n")
cat("  project/results/qc_plots/GSE10927/  [PDF + TIFF 700dpi + PNG 700dpi]\n\n")

cat("NEXT STEPS (NOT implemented here — separate scripts):\n")
cat("  1. Merge gene-level matrices on common gene symbols\n")
cat("  2. Batch correction (ComBat from sva package)\n")
cat("  3. Split into training (GSE12368) / validation (GSE10927)\n")
cat("  4. WGCNA or ML feature selection\n\n")

cat("⚠️  REMINDER: Manually verify ACA/ACC labels in:\n")
cat("     sample_annotation_1  (GSE12368)\n")
cat("     sample_annotation_2  (GSE10927)\n")
cat("  Re-run Steps 9-10 for each dataset after labeling.\n\n")

# Save session info
sink("project/results/session_info.txt")
cat("Pipeline run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
cat("Session info saved: project/results/session_info.txt\n")