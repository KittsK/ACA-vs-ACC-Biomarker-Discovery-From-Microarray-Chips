############################
# ACC GEO Merger + Stratified 70/30 Split
# Datasets  : GSE12368, GSE10927
# Platform  : GPL570 (hgu133plus2.db) — Always Confirm the Platform of the datasets before running 
# Comparison: Adrenocortical Adenoma (ACA) vs Adrenocortical Carcinoma (ACC)
# ─────────────────────────────────────────────────────
# This is the Pipeline we are going to run
#   1.  Download series matrix for each dataset
#   2.  Log2 check
#   3.  Low-expression probe filter
#   4.  Probe annotation → gene symbol (hgu133plus2.db)
#   5.  Collapse probes → gene-level (mean)
#   6.  Assign ACA/ACC by GSM number range
#   7.  Find common genes & merge
#   8.  Sanity checks
#   9.  PCA before batch correction
#   10. ComBat batch correction
#   11. PCA after batch correction
#   12. Stratified 70/30 split
#   13. Save all outputs
#
# OUTPUTS (all in the Directory that you chose):
#   expr_merged_corrected.csv     ← full batch-corrected matrix (83 samples)
#   metadata_final.csv            ← full metadata (83 rows, ACA/ACC labels)
#   expr_train_70.csv             ← training set expression (70%)
#   expr_test_30.csv              ← test set expression (30%)
#   metadata_train_70.csv         ← training set metadata
#   metadata_test_30.csv          ← test set metadata
#   PCA_before_after_batch.tiff   ← QC PCA plot (700 dpi)
#   PCA_before_after_batch.pdf
#   PCA_before_after_batch.png
#   session_info.txt
############################


# ══════════════════════════════════════════════════════
#  ⚙️  USER CONFIGURATION — only edit this block
# ══════════════════════════════════════════════════════

# ── Working directory ─────────────────────────────────
setwd("D:/ACC Project")

# ── Redirect ALL cache/temp from C: → D: ─────────────
# Prevents C: drive from filling up during GEO downloads.
LOCAL_CACHE <- "D:/ACC Project/cache"
dir.create(LOCAL_CACHE,                          recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(LOCAL_CACHE, "GEOquery"),   showWarnings = FALSE)
dir.create(file.path(LOCAL_CACHE, "BiocCache"),  showWarnings = FALSE)
dir.create(file.path(LOCAL_CACHE, "Rtmp"),       showWarnings = FALSE)

Sys.setenv(TMPDIR = file.path(LOCAL_CACHE, "Rtmp"))
Sys.setenv(TMP    = file.path(LOCAL_CACHE, "Rtmp"))
Sys.setenv(TEMP   = file.path(LOCAL_CACHE, "Rtmp"))
options(GEOquery.destdir = file.path(LOCAL_CACHE, "GEOquery"))

if (requireNamespace("BiocFileCache", quietly = TRUE)) {
  options(BiocFileCache.rdir = file.path(LOCAL_CACHE, "BiocCache"))
}

cat("Cache directories redirected to D: drive\n")

# ── GEO Dataset IDs ───────────────────────────────────
GEO_IDS <- c("GSE12368", "GSE10927")

# ── Platform map — both datasets are GPL570 ───────────
PLATFORM_MAP <- c(
  "GPL570"   = "hgu133plus2.db",
  "GPL571"   = "hgu133a.db",
  "GPL96"    = "hgu133a.db",
  "GPL97"    = "hgu133b.db",
  "GPL6244"  = "hugene10sttranscriptcluster.db",
  "GPL6947"  = "illuminaHumanv3.db",
  "GPL10558" = "illuminaHumanv4.db",
  "GPL6883"  = "illuminaHumanv3.db",
  "GPL13667" = "hgu219.db",
  "GPL8300"  = "hgu95av2.db"
)

# ── Sample group assignment by GSM number range ───────
# ACA = adrenocortical adenoma (benign)
# ACC = adrenocortical carcinoma (malignant)
#

#Please Confirm again from the website that the GSM numbers accurately conform to the groups, whether ACA or ACC
# GSE12368: ACA = GSM310459–GSM310474 (16 samples)
#           ACC = GSM310475–GSM310486 (12 samples)
# GSE10927: ACA = GSM277100–GSM277121 (22 samples)
#           ACC = GSM277122–GSM277154 (33 samples)

GSM_GROUPS <- list(
  GSE12368 = list(
    ACA = c(310459, 310474),
    ACC = c(310475, 310486)
  ),
  GSE10927 = list(
    ACA = c(277100, 277121),
    ACC = c(277122, 277154)
  )
)

# ── Expression filter thresholds ─────────────────────
# Applied per dataset BEFORE annotation.
# log2 intensity > 4.0 in >= 20% of samples (appropriate for RMA/series matrix)
EXPR_THRESHOLD  <- 4.0
SAMPLE_FRACTION <- 0.20

# ── 70/30 split settings ──────────────────────────────
TRAIN_FRACTION <- 0.70   # proportion for training set
RANDOM_SEED    <- 42     # for reproducibility of split

# ── Output directory ──────────────────────────────────
#Here you have to set your directory where you'll want the output files to go 
OUT_DIR <- "D:/ACC Project/results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_EXPR_FULL  <- file.path(OUT_DIR, "expr_merged_corrected.csv")
OUT_META_FULL  <- file.path(OUT_DIR, "metadata_final.csv")
OUT_EXPR_TRAIN <- file.path(OUT_DIR, "expr_train_70.csv")
OUT_EXPR_TEST  <- file.path(OUT_DIR, "expr_test_30.csv")
OUT_META_TRAIN <- file.path(OUT_DIR, "metadata_train_70.csv")
OUT_META_TEST  <- file.path(OUT_DIR, "metadata_test_30.csv")
OUT_PCA        <- file.path(OUT_DIR, "PCA_before_after_batch.tiff")
OUT_SESSION    <- file.path(OUT_DIR, "session_info.txt")

# ══════════════════════════════════════════════════════
#  End of user configuration
# ══════════════════════════════════════════════════════

set.seed(RANDOM_SEED)
options(timeout = 600)


############################
# 0. Install & load packages
############################

bioc_pkgs <- c("GEOquery", "AnnotationDbi", "hgu133plus2.db",
               "sva", "limma", "Biobase")
cran_pkgs <- c("ggplot2", "patchwork", "scales", "dplyr", "tibble")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(sva)
  library(limma)
  library(Biobase)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
  library(tibble)
})

cat("============================== ENVIRONMENT ==============================\n")
cat("R version    :", R.version$version.string, "\n")
cat("Working dir  :", getwd(), "\n")
cat("Random seed  :", RANDOM_SEED, "\n")
cat("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("GEO IDs      :", paste(GEO_IDS, collapse = " | "), "\n")
cat("Platform     : GPL570 (hgu133plus2.db) — both datasets\n")
cat("Comparison   : ACA (adenoma) vs ACC (carcinoma)\n")
cat("=========================================================================\n\n")


############################
# Helper — auto-detect platform
# (kept from original script for robustness)
############################

detect_platform_db <- function(gse_id, platform_map) {
  cat("  Detecting platform for", gse_id, "...\n")
  gse_tmp <- GEOquery::getGEO(
    gse_id, GSEMatrix = TRUE, getGPL = FALSE,
    destdir = file.path(LOCAL_CACHE, "GEOquery")
  )
  gpl_id <- Biobase::annotation(gse_tmp[[1]])
  cat("  Platform detected:", gpl_id, "\n")

  if (!gpl_id %in% names(platform_map)) {
    stop(paste0(
      "Platform '", gpl_id, "' not in PLATFORM_MAP.\n",
      "  → Add: \"", gpl_id, "\" = \"<db_package>\"\n",
      "  → Browse: https://bioconductor.org/packages/release/BiocViews.html#___AnnotationData"
    ))
  }

  db_pkg <- platform_map[[gpl_id]]
  cat("  Annotation package:", db_pkg, "\n")
  return(list(gpl_id = gpl_id, db_pkg = db_pkg))
}


############################
# Helper — ensure log2 scale
############################

ensure_log2 <- function(mat, gse_id) {
  rng <- range(mat, na.rm = TRUE)
  cat("[", gse_id, "] Expression range:", round(rng[1], 2),
      "to", round(rng[2], 2), "\n")
  if (rng[2] > 100) {
    cat("[", gse_id, "] → Linear scale detected. Applying log2(x + 1).\n")
    mat <- log2(mat + 1)
  } else {
    cat("[", gse_id, "] → Already log2 scale. No transformation needed.\n")
  }
  return(mat)
}


############################
# Helper — low-expression probe filter
############################

filter_probes <- function(mat, gse_id,
                           threshold = EXPR_THRESHOLD,
                           fraction  = SAMPLE_FRACTION) {
  n_samples <- ncol(mat)
  n_above   <- rowSums(mat > threshold)
  keep      <- n_above >= (fraction * n_samples)

  cat("[", gse_id, "] Probes before filter:", nrow(mat), "\n")
  cat("[", gse_id, "] Probes after  filter:", sum(keep),
      "(log2 >", threshold, "in >=",
      fraction * 100, "% of samples)\n")
  cat("[", gse_id, "] Probes removed      :",
      nrow(mat) - sum(keep), "\n")

  return(mat[keep, , drop = FALSE])
}


############################
# Helper — probe annotation & gene-level collapsing
#
# CHANGED from original:
#   Original used aggregate() with mean.
#   Now uses limma::avereps() — faster, consistent with
#   the upstream CEL preprocessing pipeline.
############################

annotate_and_aggregate <- function(mat, gse_id, db_pkg) {

  cat("[", gse_id, "] Annotating probes (", db_pkg, ")...\n")

  db_obj <- get(db_pkg)

  gene_symbols <- AnnotationDbi::mapIds(
    db_obj,
    keys      = rownames(mat),
    column    = "SYMBOL",
    keytype   = "PROBEID",
    multiVals = "first"
  )

  cat("[", gse_id, "] Probes with valid symbol:",
      sum(!is.na(gene_symbols) & gene_symbols != ""),
      "/", length(gene_symbols), "\n")

  # Remove unannotated probes
  keep_anno    <- !is.na(gene_symbols) & gene_symbols != ""
  mat_anno     <- mat[keep_anno, , drop = FALSE]
  syms_anno    <- gene_symbols[keep_anno]

  # Collapse duplicate probes per gene → mean (limma::avereps)
  mat_gene <- limma::avereps(mat_anno, ID = syms_anno)

  cat("[", gse_id, "] Unique genes after collapsing:", nrow(mat_gene), "\n")
  return(mat_gene)
}


############################
# Helper — assign ACA/ACC by GSM number range
#
# REPLACES original assign_condition() keyword approach.
# GSM number ranges are explicit and verified — no risk
# of keyword mismatches from inconsistent GEO titles.
############################

assign_condition_gsm <- function(gsm_ids, gse_id, gsm_groups) {
  # Extract numeric part: "GSM310459" → 310459
  gsm_nums <- as.integer(gsub("GSM", "", gsm_ids, ignore.case = TRUE))
  grp_def  <- gsm_groups[[gse_id]]

  condition <- dplyr::case_when(
    gsm_nums >= grp_def$ACA[1] & gsm_nums <= grp_def$ACA[2] ~ "ACA",
    gsm_nums >= grp_def$ACC[1] & gsm_nums <= grp_def$ACC[2] ~ "ACC",
    TRUE ~ "exclude"
  )
  return(condition)
}


############################
# Helper — PCA plot
#
# CHANGED from original:
#   Points filled by condition (ACA=blue, ACC=red)
#   Points outlined by dataset colour (green/orange)
#   Ellipses drawn per dataset (shows batch effect clearly)
############################

plot_pca <- function(mat, meta, title_label, seed) {
  set.seed(seed)

  pca    <- prcomp(t(mat), scale. = TRUE, center = TRUE)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$sample_id <- rownames(scores)
  scores <- merge(scores,
                  meta[, c("sample_id", "batch", "condition")],
                  by = "sample_id")

  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  batch_colors <- c(GSE12368 = "#1B7837", GSE10927 = "#E08214")
  cond_fills   <- c(ACA = "#2166AC", ACC = "#D6604D")
  cond_shapes  <- c(ACA = 21, ACC = 24)

  ggplot(scores, aes(x = PC1, y = PC2)) +

    # Dataset ellipses (dashed) — reveals batch clustering
    stat_ellipse(aes(colour = batch, fill = batch),
                 geom = "polygon", type = "norm",
                 level = 0.95, alpha = 0.08,
                 show.legend = FALSE) +
    stat_ellipse(aes(colour = batch),
                 type = "norm", level = 0.95,
                 linewidth = 0.9, linetype = "dashed",
                 show.legend = FALSE) +

    # Reference lines at origin
    geom_hline(yintercept = 0, colour = "grey80",
               linewidth = 0.3, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "grey80",
               linewidth = 0.3, linetype = "dotted") +

    # Points — filled by condition, outlined by dataset
    geom_point(aes(fill = condition,
                   colour = batch,
                   shape  = condition),
               size = 3.5, stroke = 0.9, alpha = 0.92) +

    scale_colour_manual(values = batch_colors, name = "Dataset") +
    scale_fill_manual(  values = cond_fills,   name = "Diagnosis") +
    scale_shape_manual( values = cond_shapes,  name = "Diagnosis") +

    labs(
      title = title_label,
      x     = paste0("PC1: ", var_exp[1], "% variance"),
      y     = paste0("PC2: ", var_exp[2], "% variance")
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title     = element_text(face = "bold"),
      panel.border     = element_rect(colour = "black", fill = NA,
                                       linewidth = 1.2),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.35,
                                       linetype = "dashed"),
      legend.position  = "bottom",
      legend.box       = "horizontal"
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 4)),
      fill   = guide_legend(override.aes = list(size = 4, shape = 21)),
      shape  = guide_legend(override.aes = list(size = 4))
    )
}


############################
# 1. Detect platforms
############################

cat("\n========== Platform Detection ==========\n")
platform_info <- lapply(GEO_IDS, detect_platform_db,
                         platform_map = PLATFORM_MAP)
names(platform_info) <- GEO_IDS

cat("\nPlatform summary:\n")
for (id in GEO_IDS) {
  cat(" ", id, "→", platform_info[[id]]$gpl_id,
      "→", platform_info[[id]]$db_pkg, "\n")
}


############################
# 2. Download & process each dataset INDEPENDENTLY
#    Each dataset is fully processed before anything is merged.
############################

# ══════════════════════════════════════════════════════
# SETTING THE EXPRESSION THRESHOLDS 
# ══════════════════════════════════════════════════════
#We will update the expression thresholds again based on the expression range. Note that both the GSE12368 and GSE10927 use different expression thresholds and this based on the different expression ranges between the two samples

EXPR_THRESHOLDS <- list(
  GSE12368 = 4.0,   # after log2(x+1) transform → range 0–17 → keeps ~98% ✅
  GSE10927 = 2.3    # already log2, compressed   → range 1.7–4.74 → keeps ~66.5% ✅
)
SAMPLE_FRACTION <- 0.20

# ══════════════════════════════════════════════════════
# UPDATE THE filter_probes() FUNCTION AGAIN
# — This will now accept per-dataset threshold argument
# ══════════════════════════════════════════════════════

filter_probes <- function(mat, gse_id,
                          threshold = 4.0,
                          fraction  = SAMPLE_FRACTION) {
  n_samples <- ncol(mat)
  n_above   <- rowSums(mat > threshold)
  keep      <- n_above >= (fraction * n_samples)
  
  cat("[", gse_id, "] Probes before filter:", nrow(mat), "\n")
  cat("[", gse_id, "] Probes after  filter:", sum(keep),
      "(log2 >", threshold, "in >=", fraction * 100, "% of samples)\n")
  cat("[", gse_id, "] Probes removed      :", nrow(mat) - sum(keep), "\n")
  
  return(mat[keep, , drop = FALSE])
}

# 2. Download & process each dataset INDEPENDENTLY
# Note that Each dataset is fully processed before anything is merged.
#
#    KEY CHANGE in this section:
#    filter_probes() now receives EXPR_THRESHOLDS[[gse_id]]
#    so each dataset uses its own appropriate threshold:
#      GSE12368 → 4.0  (standard RMA scale after log2 transform)
#      GSE10927 → 2.3  (compressed pre-normalised scale)
############################

expr_list <- list()
meta_list <- list()

for (gse_id in GEO_IDS) {
  
  cat("\n========== Processing", gse_id, "==========\n")
  
  db_pkg <- platform_info[[gse_id]]$db_pkg
  gpl_id <- platform_info[[gse_id]]$gpl_id
  
  # ── 2a. Download ────────────────────────────────────
  cat("[", gse_id, "] Downloading series matrix...\n")
  gse_data <- GEOquery::getGEO(
    gse_id, GSEMatrix = TRUE, getGPL = FALSE,
    destdir = file.path(LOCAL_CACHE, "GEOquery")
  )
  eset  <- gse_data[[1]]
  pheno <- Biobase::pData(eset)
  mat   <- Biobase::exprs(eset)
  
  cat("[", gse_id, "] Downloaded:", nrow(mat), "probes x",
      ncol(mat), "samples\n")
  
  # ── 2b. Log2 check ──────────────────────────────────
  # Detects linear scale (max > 100) and applies log2(x+1)
  # GSE12368: linear → transformed  | GSE10927: already log2 → unchanged
  mat <- ensure_log2(mat, gse_id)
  
  # ── 2c. Probe filter ────────────────────────────────
  # CHANGED: now passes dataset-specific threshold
  # GSE12368 uses 4.0 | GSE10927 uses 2.3
  # Both thresholds retain ~60-70% of probes on their respective scales
  mat <- filter_probes(mat, gse_id,
                       threshold = EXPR_THRESHOLDS[[gse_id]])
  
  # ── 2d. Annotate & collapse to gene level ───────────
  expr_agg <- annotate_and_aggregate(mat, gse_id, db_pkg)
  
  # ── 2e. Assign ACA / ACC by GSM range ───────────────
  gsm_ids       <- rownames(pheno)
  condition_vec <- assign_condition_gsm(gsm_ids, gse_id, GSM_GROUPS)
  
  meta_gse <- data.frame(
    sample_id = gsm_ids,
    title     = as.character(pheno$title),
    condition = condition_vec,
    batch     = gse_id,
    dataset   = gse_id,
    platform  = gpl_id,
    stringsAsFactors = FALSE
  )
  
  cat("[", gse_id, "] Condition breakdown (before exclusion):\n")
  print(table(meta_gse$condition))
  
  # Report excluded samples (outside defined GSM ranges)
  n_excluded <- sum(meta_gse$condition == "exclude")
  if (n_excluded > 0) {
    cat("[", gse_id, "] Excluding", n_excluded,
        "samples outside defined GSM ranges:\n")
    excl <- meta_gse[meta_gse$condition == "exclude",
                     c("sample_id", "title")]
    for (i in seq_len(nrow(excl)))
      cat("    -", excl$sample_id[i], "|", excl$title[i], "\n")
  }
  
  # Validate both classes present
  if (sum(meta_gse$condition == "ACA") == 0)
    warning(paste0("[", gse_id, "] ⚠️  No ACA samples found. Check GSM_GROUPS."))
  if (sum(meta_gse$condition == "ACC") == 0)
    warning(paste0("[", gse_id, "] ⚠️  No ACC samples found. Check GSM_GROUPS."))
  
  # ── 2f. Keep ACA + ACC only ─────────────────────────
  keep              <- meta_gse$condition %in% c("ACA", "ACC")
  meta_gse_filtered <- meta_gse[keep, ]
  expr_filtered     <- expr_agg[, meta_gse_filtered$sample_id, drop = FALSE]
  
  cat("[", gse_id, "] Retained — ACA:",
      sum(meta_gse_filtered$condition == "ACA"),
      "| ACC:", sum(meta_gse_filtered$condition == "ACC"), "\n")
  
  expr_list[[gse_id]] <- expr_filtered
  meta_list[[gse_id]] <- meta_gse_filtered
  
  cat("[", gse_id, "] ✓ Done.\n")
}

############################
# 3. Common genes & merge
############################

cat("\n========== Gene Intersection ==========\n")
for (id in GEO_IDS)
  cat(id, ":", nrow(expr_list[[id]]), "genes\n")

common_genes <- Reduce(intersect, lapply(expr_list, rownames))
cat("Common genes across both datasets:", length(common_genes), "\n")

# Both GPL570 — expect >10,000 common genes
if (length(common_genes) < 5000)
  stop(paste(
    "Only", length(common_genes), "common genes found.",
    "Expected >10,000 since both datasets use GPL570.",
    "Check annotation step."
  ))

cat("\n========== Merging Datasets ==========\n")
expr_sub      <- lapply(expr_list, function(m) as.matrix(m[common_genes, ]))
merged_mat    <- do.call(cbind, expr_sub)
meta_combined <- do.call(rbind, meta_list)
rownames(meta_combined) <- NULL

meta_ordered <- meta_combined[
  match(colnames(merged_mat), meta_combined$sample_id), ]
rownames(meta_ordered) <- meta_ordered$sample_id

cat("Merged matrix:", nrow(merged_mat), "genes x",
    ncol(merged_mat), "samples\n")
cat("\nSample breakdown:\n")
print(table(meta_ordered$condition, meta_ordered$dataset))


############################
# 4. Sanity checks
############################

cat("\n========== Sanity Checks ==========\n")
stopifnot(
  "Column/metadata order mismatch" =
    all(colnames(merged_mat) == meta_ordered$sample_id),
  "NA condition labels" =
    !any(is.na(meta_ordered$condition)),
  "NA batch labels" =
    !any(is.na(meta_ordered$batch)),
  "NA values in expression matrix" =
    !any(is.na(merged_mat))
)
cat("✅ All sanity checks passed.\n")


############################
# 5. PCA — before batch correction
############################

cat("\n========== PCA Before Batch Correction ==========\n")
p_before <- plot_pca(merged_mat, meta_ordered,
                     "Before Batch Correction", RANDOM_SEED)
cat("PCA (before) computed.\n")


############################
# 6. ComBat batch correction
############################

cat("\n========== ComBat Batch Correction ==========\n")

batch_vec     <- meta_ordered$batch
condition_vec <- meta_ordered$condition

if (length(unique(batch_vec)) < 2) {
  cat("Only 1 batch — skipping ComBat.\n")
  expr_corrected <- merged_mat
} else {
  # Protect biological variance by including condition in model
  mod <- model.matrix(~ condition_vec)
  set.seed(RANDOM_SEED)
  expr_corrected <- sva::ComBat(
    dat         = as.matrix(merged_mat),
    batch       = batch_vec,
    mod         = mod,
    par.prior   = TRUE,
    prior.plots = FALSE
  )
  cat("ComBat complete:", nrow(expr_corrected), "genes x",
      ncol(expr_corrected), "samples\n")
}


############################
# 7. PCA — after batch correction
############################

cat("\n========== PCA After Batch Correction ==========\n")
p_after <- plot_pca(expr_corrected, meta_ordered,
                    "After Batch Correction", RANDOM_SEED)
cat("PCA (after) computed.\n")


############################
# 8. STRATIFIED 70/30 SPLIT
#
# Stratified = the ACA/ACC ratio is preserved in BOTH
# the training and test sets.
#
# Method:
#   - Split ACA samples independently: 70% train, 30% test
#   - Split ACC samples independently: 70% train, 30% test
#   - Combine → balanced training set + balanced test set
#
# This will ensure neither set is accidentally ACA-heavy or ACC-heavy,
# which is critical for unbiased DEG and ML performance.
############################

cat("\n========== Stratified 70/30 Split ==========\n")

set.seed(RANDOM_SEED)

# ── Split ACA samples ────────────────────────────────
aca_ids   <- meta_ordered$sample_id[meta_ordered$condition == "ACA"]
n_aca     <- length(aca_ids)
n_aca_trn <- round(n_aca * TRAIN_FRACTION)

aca_train_idx <- sample(seq_len(n_aca), size = n_aca_trn, replace = FALSE)
aca_train_ids <- aca_ids[aca_train_idx]
aca_test_ids  <- aca_ids[-aca_train_idx]

# ── Split ACC samples ────────────────────────────────
acc_ids   <- meta_ordered$sample_id[meta_ordered$condition == "ACC"]
n_acc     <- length(acc_ids)
n_acc_trn <- round(n_acc * TRAIN_FRACTION)

acc_train_idx <- sample(seq_len(n_acc), size = n_acc_trn, replace = FALSE)
acc_train_ids <- acc_ids[acc_train_idx]
acc_test_ids  <- acc_ids[-acc_train_idx]

# ── Combine ──────────────────────────────────────────
train_ids <- c(aca_train_ids, acc_train_ids)
test_ids  <- c(aca_test_ids,  acc_test_ids)

# ── Subset expression matrices ───────────────────────
expr_train <- expr_corrected[, train_ids, drop = FALSE]
expr_test  <- expr_corrected[, test_ids,  drop = FALSE]

# ── Subset metadata ──────────────────────────────────
meta_train <- meta_ordered[train_ids, ]
meta_test  <- meta_ordered[test_ids,  ]

# ── Report ───────────────────────────────────────────
cat("\nTraining set (70%):\n")
cat("  Total samples :", ncol(expr_train), "\n")
cat("  ACA           :", sum(meta_train$condition == "ACA"), "\n")
cat("  ACC           :", sum(meta_train$condition == "ACC"), "\n")
cat("  ACA %         :",
    round(100 * sum(meta_train$condition == "ACA") / ncol(expr_train), 1), "%\n")

cat("\nTest set (30%):\n")
cat("  Total samples :", ncol(expr_test), "\n")
cat("  ACA           :", sum(meta_test$condition == "ACA"), "\n")
cat("  ACC           :", sum(meta_test$condition == "ACC"), "\n")
cat("  ACA %         :",
    round(100 * sum(meta_test$condition == "ACA") / ncol(expr_test), 1), "%\n")

cat("\nOriginal ACA % (full set):",
    round(100 * n_aca / (n_aca + n_acc), 1), "%\n")
cat("Training ACA %           :",
    round(100 * sum(meta_train$condition == "ACA") / ncol(expr_train), 1), "%\n")
cat("Test ACA %               :",
    round(100 * sum(meta_test$condition == "ACA") / ncol(expr_test), 1), "%\n")
cat("(These three should be approximately equal — confirms stratification)\n")

# ── Sanity check ─────────────────────────────────────
stopifnot(
  "Train + test samples do not sum to total" =
    (ncol(expr_train) + ncol(expr_test)) == ncol(expr_corrected),
  "Train column/metadata mismatch" =
    all(colnames(expr_train) == rownames(meta_train)),
  "Test column/metadata mismatch" =
    all(colnames(expr_test) == rownames(meta_test))
)
cat("✅ Split sanity checks passed.\n")


############################
# 9. Save all outputs
############################

cat("\n========== Saving Outputs ==========\n")

# ── PCA figure ───────────────────────────────────────
combined_pca <- (p_before | p_after) +
  patchwork::plot_annotation(
    title    = "GSE12368 & GSE10927 \u2014 Batch Correction QC",
    subtitle = paste0(
      "ACA (adenoma) vs ACC (carcinoma)  \u2022  ",
      "GPL570 (Affymetrix HG-U133 Plus 2.0)  \u2022  ",
      "ComBat batch correction  \u2022  ",
      "Generated: ", format(Sys.Date(), "%Y-%m-%d")
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40")
    )
  )

ggsave(OUT_PCA, plot = combined_pca,
       width = 16, height = 8, dpi = 700,
       device = "tiff", compression = "lzw")

ggsave(sub("\\.tiff$", ".pdf", OUT_PCA),
       plot = combined_pca, width = 16, height = 8,
       device = cairo_pdf)

ggsave(sub("\\.tiff$", ".png", OUT_PCA),
       plot = combined_pca, width = 16, height = 8, dpi = 700)

cat("✅ Saved: PCA plots (TIFF / PDF / PNG)\n")

# ── Full merged matrix ────────────────────────────────
write.csv(expr_corrected, file = OUT_EXPR_FULL, row.names = TRUE)
cat("✅ Saved:", OUT_EXPR_FULL, "\n")

write.csv(meta_ordered, file = OUT_META_FULL, row.names = FALSE)
cat("✅ Saved:", OUT_META_FULL, "\n")

# ── Training set ─────────────────────────────────────
write.csv(expr_train, file = OUT_EXPR_TRAIN, row.names = TRUE)
cat("✅ Saved:", OUT_EXPR_TRAIN, "\n")

write.csv(meta_train, file = OUT_META_TRAIN, row.names = FALSE)
cat("✅ Saved:", OUT_META_TRAIN, "\n")

# ── Test set ──────────────────────────────────────────
write.csv(expr_test, file = OUT_EXPR_TEST, row.names = TRUE)
cat("✅ Saved:", OUT_EXPR_TEST, "\n")

write.csv(meta_test, file = OUT_META_TEST, row.names = FALSE)
cat("✅ Saved:", OUT_META_TEST, "\n")


############################
# 10. Session info
############################

sink(OUT_SESSION)
cat("Script     : acc_geo_merger_split.R\n")
cat("Run time   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Random seed:", RANDOM_SEED, "\n")
cat("GEO IDs    : GSE12368 | GSE10927\n")
cat("Platform   : GPL570 (hgu133plus2.db)\n")
cat("Comparison : ACA vs ACC\n\n")
cat("GSM group ranges:\n")
cat("  GSE12368 ACA: GSM310459 - GSM310474 (16 samples)\n")
cat("  GSE12368 ACC: GSM310475 - GSM310486 (12 samples)\n")
cat("  GSE10927 ACA: GSM277100 - GSM277121 (22 samples)\n")
cat("  GSE10927 ACC: GSM277122 - GSM277154 (33 samples)\n\n")
cat("Split settings:\n")
cat("  Method        : Stratified (ACA/ACC ratio preserved)\n")
cat("  Train fraction:", TRAIN_FRACTION, "\n")
cat("  Training set  :", ncol(expr_train), "samples\n")
cat("  Test set      :", ncol(expr_test),  "samples\n\n")
print(sessionInfo())
sink()
cat("✅ Saved:", OUT_SESSION, "\n")


############################
# 11. Final summary
############################

cat("\n============================== DONE ==============================\n")
cat("\nOutput directory:", OUT_DIR, "\n\n")
cat("Files saved:\n")
cat("  Full merged matrix    →", basename(OUT_EXPR_FULL),  "\n")
cat("  Full metadata         →", basename(OUT_META_FULL),  "\n")
cat("  Training expression   →", basename(OUT_EXPR_TRAIN), "\n")
cat("  Training metadata     →", basename(OUT_META_TRAIN), "\n")
cat("  Test expression       →", basename(OUT_EXPR_TEST),  "\n")
cat("  Test metadata         →", basename(OUT_META_TEST),  "\n")
cat("  PCA plot (TIFF/PDF/PNG) → PCA_before_after_batch.*\n")
cat("  Session info          →", basename(OUT_SESSION),    "\n")

cat("\n--- Full merged matrix ---\n")
cat("  Genes         :", nrow(expr_corrected), "\n")
cat("  Total samples :", ncol(expr_corrected), "\n")
cat("  ACA           :", sum(meta_ordered$condition == "ACA"), "\n")
cat("  ACC           :", sum(meta_ordered$condition == "ACC"), "\n")

cat("\n--- Training set (70%) ---\n")
cat("  Total samples :", ncol(expr_train), "\n")
cat("  ACA           :", sum(meta_train$condition == "ACA"), "\n")
cat("  ACC           :", sum(meta_train$condition == "ACC"), "\n")

cat("\n--- Test set (30%) ---\n")
cat("  Total samples :", ncol(expr_test), "\n")
cat("  ACA           :", sum(meta_test$condition == "ACA"), "\n")
cat("  ACC           :", sum(meta_test$condition == "ACC"), "\n")

cat("\n--- Per dataset breakdown (full set) ---\n")
print(table(meta_ordered$condition, meta_ordered$dataset))

cat("\n⚠️  NEXT STEP:\n")
cat("  DEG analysis → use expr_train_70.csv + metadata_train_70.csv\n")
cat("  Validation   → use expr_test_30.csv  + metadata_test_30.csv\n")
cat("==================================================================\n")
