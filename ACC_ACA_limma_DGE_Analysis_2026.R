################################################################################
# LIMMA DIFFERENTIAL EXPRESSION ANALYSIS
# Project   : Adrenocortical Adenoma (ACA) vs Adrenocortical Carcinoma (ACC)
# Input     : 70% training split from merged & batch-corrected dataset
#             (GSE12368 + GSE10927)
# Contrast  : ACC vs ACA (positive logFC = higher in ACC)
# Thresholds: adj.P.Val < 0.05 & |logFC| > 1
#
# OUTPUTS (all in D:/ACC Project/results/DEG/):
#   tables/
#     Limma_results_all.csv           — all genes with stats
#     Limma_results_significant.csv   — DEGs only (FDR<0.05, |logFC|>1)
#     Limma_upregulated.csv           — genes UP in ACC vs ACA
#     Limma_downregulated.csv         — genes DOWN in ACC vs ACA
#     GO_BP_results.csv               — GO Biological Process enrichment
#     GO_CC_results.csv               — GO Cellular Component enrichment
#     GO_MF_results.csv               — GO Molecular Function enrichment
#     KEGG_results.csv                — KEGG pathway enrichment
#   plots/
#     PCA_plot.tif                    — PCA of training samples
#     MA_plot.tif                     — MA plot
#     Volcano_plot.tif                — Volcano plot (top 10 up + 10 down labelled)
#     Heatmap_top_DEGs.tif            — Heatmap of top 50 DEGs
#     GO_BP_dotplot.tif
#     GO_CC_dotplot.tif
#     GO_MF_dotplot.tif
#     KEGG_dotplot.tif
#
# PREREQUISITES:
#   expr_train_70.csv       — from acc_geo_merger_split.R
#   metadata_train_70.csv   — from acc_geo_merger_split.R
################################################################################


# =====================================================
# 0. WORKING DIRECTORY & SETUP
# =====================================================

setwd("D:/ACC Project")

# Reproducibility
set.seed(42)

# =====================================================
# 1. LOAD LIBRARIES
# =====================================================

bioc_pkgs <- c("limma", "ComplexHeatmap", "circlize",
               "clusterProfiler", "org.Hs.eg.db", "enrichplot")
cran_pkgs <- c("tidyverse", "ggrepel", "grid")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# =====================================================
# 2. USER INPUTS
# =====================================================

# Input files — 70% training split from merger script
expr_file     <- "D:/ACC Project/results/expr_train_70.csv"
metadata_file <- "D:/ACC Project/results/metadata_train_70.csv"

# Output directory
base_dir      <- "D:/ACC Project/results/DEG"

# DEG thresholds
PVAL_THRESH <- 0.05   # adjusted p-value cutoff (BH/FDR)
FC_THRESH   <- 1.0    # |log2 fold change| cutoff

# Volcano plot: top N labelled per direction
N_VOLCANO   <- 10     # top 10 up + top 10 down

# Heatmap: top N DEGs
N_HEATMAP   <- 50

# Condition colours — consistent with PCA plots
#   ACA = Blue  (#2166AC)
#   ACC = Red   (#D6604D)
COND_COLS <- c(ACA = "#2166AC", ACC = "#D6604D")

# =====================================================
# 3. LOAD DATA
# =====================================================

cat("\n========== Loading Data ==========\n")

expr_data <- read.csv(expr_file,     row.names = 1, check.names = FALSE)
metadata  <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

# Align columns — ensure expression matrix columns match metadata rows
common_samples <- intersect(colnames(expr_data), rownames(metadata))
expr_data      <- expr_data[, common_samples, drop = FALSE]
metadata       <- metadata[common_samples,  , drop = FALSE]

cat("✅ Data loaded and aligned:\n")
cat("   Genes          :", nrow(expr_data), "\n")
cat("   Samples        :", ncol(expr_data), "\n")
cat("   ACA samples    :", sum(metadata$condition == "ACA"), "\n")
cat("   ACC samples    :", sum(metadata$condition == "ACC"), "\n")
cat("   Expression range (log2):", round(min(expr_data), 2),
    "to", round(max(expr_data), 2), "\n")

# =====================================================
# 4. FILTER NON-CODING GENES
# =====================================================

cat("\n========== Filtering Non-Coding Genes ==========\n")

gene_names <- rownames(expr_data)

# Patterns identifying non-coding / unannotated genes
# These are retained from the original script and are appropriate
# for Affymetrix HG-U133 Plus 2.0 annotations
remove_patterns <- c(
  "^LOC", "^MT-", "^MIR", "^LINC", "^SNORD", "^SNORA", "^SCARNA",
  "^snoRNA", "^miR-", "^lnc", "^RN7", "^RMRP", "^RPPH1", "^TERC",
  "^XIST", "^TSIX", "^NEAT1", "^MALAT1", "^H19", "^MEG", "^PTEN-",
  "^RP[0-9]", "^AC[0-9]", "^AL[0-9]", "^AP[0-9]", "^CTD-", "^CTB-",
  "^CTC-", "^LL[0-9]", "^KB-", "^CH[0-9]", "^C1orf"
)

remove_pattern_combined <- paste(remove_patterns, collapse = "|")
keep_genes  <- !grepl(remove_pattern_combined, gene_names, ignore.case = TRUE)
expr_data   <- expr_data[keep_genes, ]

cat("✅ Non-coding genes removed :", sum(!keep_genes), "\n")
cat("✅ Protein-coding genes kept:", nrow(expr_data), "\n")

# =====================================================
# 5. PREPARE METADATA
# =====================================================

cat("\n========== Preparing Metadata ==========\n")

# Set factor levels: ACA = reference, ACC = comparison
# Positive logFC will mean HIGHER in ACC vs ACA
metadata$condition <- factor(metadata$condition, levels = c("ACA", "ACC"))

cat("Sample distribution in training set:\n")
print(table(metadata$condition))

# =====================================================
# 6. CREATE OUTPUT DIRECTORIES
# =====================================================

dir.create(base_dir,                            showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "plots"),        showWarnings = FALSE)
dir.create(file.path(base_dir, "tables"),       showWarnings = FALSE)

cat("✅ Output directories created:", base_dir, "\n")

# =====================================================
# 7. LIMMA DIFFERENTIAL EXPRESSION ANALYSIS
# =====================================================

cat("\n========== Running limma DEG Analysis ==========\n")
cat("Contrast: ACC vs ACA\n")
cat("Positive logFC = higher expression in ACC (carcinoma)\n")
cat("Negative logFC = higher expression in ACA (adenoma)\n\n")

# Design matrix — no intercept, one column per condition
design <- model.matrix(~ 0 + condition, data = metadata)
colnames(design) <- levels(metadata$condition)   # → "ACA", "ACC"

cat("Design matrix columns:", colnames(design), "\n")
cat("Design matrix dimensions:", dim(design), "\n")

# Fit linear model
fit <- lmFit(expr_data, design)

# Define contrast: ACC - ACA
# Positive logFC = upregulated in ACC (carcinoma)
# Negative logFC = upregulated in ACA (adenoma)
contrast.matrix <- makeContrasts(
  ACCvsACA = ACC - ACA,
  levels   = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract all results (sorted by adjusted p-value)
res <- topTable(fit2,
                coef          = "ACCvsACA",
                number        = Inf,
                adjust.method = "BH",
                sort.by       = "P")

cat("✅ limma analysis complete\n")
cat("   Log2FC range:", round(min(res$logFC), 2),
    "to", round(max(res$logFC), 2), "\n")

# =====================================================
# 8. PROCESS & CLASSIFY RESULTS
# =====================================================

cat("\n========== Processing Results ==========\n")

res$gene <- rownames(res)

# Classify each gene
res$regulation <- ifelse(
  !is.na(res$adj.P.Val) &
    res$adj.P.Val < PVAL_THRESH &
    abs(res$logFC) > FC_THRESH,
  ifelse(res$logFC > FC_THRESH, "Up", "Down"),
  "NS"
)

# Label significant genes (used in volcano plot)
res$label <- ifelse(res$regulation != "NS", res$gene, "")

# Subsets
res_sig  <- res[res$regulation != "NS", ]
res_up   <- res[res$regulation == "Up",   ]
res_down <- res[res$regulation == "Down", ]

cat("✅ DEG summary:\n")
cat("   Total DEGs              :", nrow(res_sig),  "\n")
cat("   Upregulated in ACC      :", nrow(res_up),   "\n")
cat("   Downregulated in ACC    :", nrow(res_down), "\n")
cat("   (Downregulated = higher in ACA)\n")

# ── Save result tables ──────────────────────────────
write.csv(res,
          file.path(base_dir, "tables", "Limma_ACCvsACA_all.csv"),
          row.names = FALSE)
write.csv(res_sig,
          file.path(base_dir, "tables", "Limma_ACCvsACA_significant.csv"),
          row.names = FALSE)
write.csv(res_up,
          file.path(base_dir, "tables", "Limma_ACCvsACA_upregulated.csv"),
          row.names = FALSE)
write.csv(res_down,
          file.path(base_dir, "tables", "Limma_ACCvsACA_downregulated.csv"),
          row.names = FALSE)

cat("✅ Result tables saved to:", file.path(base_dir, "tables"), "\n")

# =====================================================
# 9. PCA PLOT (training set)
# =====================================================

cat("\n========== PCA Plot ==========\n")

pca_res    <- prcomp(t(expr_data), scale. = TRUE)
pca_df     <- as.data.frame(pca_res$x[, 1:2])
pca_df$condition <- metadata$condition
percentVar <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

pca_plot <- ggplot(pca_df, aes(PC1, PC2,
                                color = condition,
                                fill  = condition)) +
  geom_point(size = 3.5, alpha = 0.88, shape = 21,
             colour = "white", stroke = 0.6) +
  geom_point(size = 3.5, alpha = 0.88,
             aes(shape = condition)) +
  stat_ellipse(geom = "polygon", level = 0.95,
               alpha = 0.12, colour = NA) +
  stat_ellipse(level = 0.95, linewidth = 1.1) +
  scale_color_manual(name   = "Diagnosis", values = COND_COLS) +
  scale_fill_manual( name   = "Diagnosis", values = COND_COLS) +
  scale_shape_manual(name   = "Diagnosis",
                     values = c(ACA = 21, ACC = 24)) +
  theme_bw(base_size = 14) +
  theme(
    legend.title          = element_text(size = 14, face = "bold"),
    legend.text           = element_text(size = 12, face = "bold"),
    axis.title            = element_text(size = 13, face = "bold"),
    axis.text             = element_text(size = 11, face = "bold"),
    panel.grid.major      = element_blank(),
    panel.grid.minor      = element_blank(),
    legend.justification  = c("right", "top"),
    legend.background     = element_rect(color = "black",
                                          linewidth = 1, fill = "white"),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.margin         = margin(6, 6, 6, 6)
  ) +
  labs(
    title    = "PCA — Training Set (70%)",
    subtitle = paste0("ACA n=", sum(metadata$condition == "ACA"),
                      "  |  ACC n=", sum(metadata$condition == "ACC")),
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )

ggsave(file.path(base_dir, "plots", "PCA_plot.tif"),
       pca_plot, width = 11, height = 7, dpi = 600)
cat("✅ PCA plot saved\n")

# =====================================================
# 10. MA PLOT
# =====================================================

cat("\n========== MA Plot ==========\n")

ma_plot <- ggplot(res, aes(AveExpr, logFC, color = regulation)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = c(-FC_THRESH, FC_THRESH),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0,
             linetype = "solid",  color = "black", linewidth = 0.3) +
  scale_color_manual(
    name   = "Regulation",
    values = c("Up" = "#D6604D", "Down" = "#2166AC", "NS" = "grey70")
    # Up   = RED  (higher in ACC — carcinoma)
    # Down = BLUE (higher in ACA — adenoma)
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.title          = element_text(face = "bold"),
    axis.title            = element_text(face = "bold"),
    axis.text             = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    panel.grid.minor      = element_blank(),
    legend.justification  = c("right", "top"),
    legend.background     = element_rect(color = "black",
                                          linewidth = 1, fill = "white"),
    legend.box.background = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "MA Plot — ACC vs ACA",
    x     = "Average expression (log2)",
    y     = paste0("Log2 fold change (ACC vs ACA)\n",
                   "Positive = higher in ACC (carcinoma)")
  )

ggsave(file.path(base_dir, "plots", "MA_plot.tif"),
       ma_plot, width = 8, height = 6, dpi = 600)
cat("✅ MA plot saved\n")

# =====================================================
# 11. VOLCANO PLOT
# =====================================================

cat("\n========== Volcano Plot ==========\n")

# Top 10 up + top 10 down (by adjusted p-value)
top_up   <- res %>%
  filter(regulation == "Up")   %>%
  arrange(adj.P.Val)           %>%
  slice_head(n = N_VOLCANO)

top_down <- res %>%
  filter(regulation == "Down") %>%
  arrange(adj.P.Val)           %>%
  slice_head(n = N_VOLCANO)

top_labelled <- bind_rows(top_up, top_down)

volcano_plot <- ggplot(res, aes(logFC, -log10(adj.P.Val),
                                 color = regulation)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(PVAL_THRESH),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(
    data          = top_labelled,
    aes(label     = label),
    size          = 3.5,
    fontface      = "bold",
    max.overlaps  = 25,
    box.padding   = 0.4,
    point.padding = 0.3,
    segment.size  = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_manual(
    name   = "Regulation",
    values = c("Up" = "#D6604D", "Down" = "#2166AC", "NS" = "grey70")
    # Up   = RED  → upregulated in ACC (carcinoma)
    # Down = BLUE → upregulated in ACA (adenoma)
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x           = element_text(face = "bold"),
    axis.text.y           = element_text(face = "bold"),
    axis.title.x          = element_text(face = "bold"),
    axis.title.y          = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.justification  = c("right", "top"),
    legend.background     = element_rect(color = "black",
                                          linewidth = 1, fill = "white"),
    legend.box.background = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title    = "Volcano Plot — ACC vs ACA",
    subtitle = paste0("DEGs: ", nrow(res_up), " up in ACC  |  ",
                      nrow(res_down), " up in ACA",
                      "  (FDR < ", PVAL_THRESH,
                      ", |log2FC| > ", FC_THRESH, ")"),
    x        = "Log2 Fold Change (ACC vs ACA)",
    y        = expression(-log[10](italic(p)-value~adjusted))
  )

ggsave(file.path(base_dir, "plots", "Volcano_plot.tif"),
       volcano_plot, width = 8, height = 7, dpi = 600)
cat("✅ Volcano plot saved\n")

# =====================================================
# 12. HEATMAP — TOP 50 DEGs
# =====================================================

cat("\n========== Heatmap ==========\n")

# Select top N DEGs by adjusted p-value (FDR < 0.01 first, fall back to 0.05)
top_degs <- res %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.01, abs(logFC) > FC_THRESH) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = N_HEATMAP)

if (nrow(top_degs) < N_HEATMAP) {
  cat("  Less than", N_HEATMAP, "DEGs at FDR<0.01.",
      "Relaxing to FDR<0.05...\n")
  top_degs <- res %>%
    filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > FC_THRESH) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = N_HEATMAP)
}

cat("  DEGs selected for heatmap:", nrow(top_degs), "\n")

if (nrow(top_degs) > 0) {

  heatmap_matrix <- as.matrix(
    expr_data[top_degs$gene, rownames(metadata), drop = FALSE]
  )

  # Z-score scaling per gene (row)
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

  # Remove any rows with NA or Inf after scaling
  valid_rows <- apply(heatmap_matrix_scaled, 1,
                       function(x) !all(is.na(x)) & !any(is.infinite(x)))
  heatmap_matrix_scaled <- heatmap_matrix_scaled[valid_rows, , drop = FALSE]

  # Condition annotation bar
  # ACA = Blue (#2166AC) | ACC = Red (#D6604D) — consistent with all other plots
  annotation_col    <- data.frame(
    Diagnosis = metadata$condition,
    row.names = rownames(metadata)
  )
  annotation_colors <- list(
    Diagnosis = c(ACA = "#2166AC", ACC = "#D6604D")
  )

  tiff(file.path(base_dir, "plots", "Heatmap_top_DEGs.tif"),
       width = 10, height = 8, units = "in", res = 600)

  # Colour scale: blue (low) → white (mid) → red (high)
  col_fun  <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D"))

  row_dend <- as.dendrogram(hclust(dist(heatmap_matrix_scaled),
                                    method = "complete"))
  col_dend <- as.dendrogram(hclust(dist(t(heatmap_matrix_scaled)),
                                    method = "complete"))

  # Top annotation bar
  ha <- HeatmapAnnotation(
    Diagnosis = annotation_col$Diagnosis,
    col       = annotation_colors,
    annotation_name_gp     = gpar(fontsize = 10, fontface = "bold"),
    show_annotation_name   = TRUE,
    annotation_legend_param = list(
      Diagnosis = list(
        title_gp  = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 10, fontface = "bold"),
        border    = "black"
      )
    )
  )

  ht <- Heatmap(
    heatmap_matrix_scaled,
    name            = "Z-score",
    col             = col_fun,
    top_annotation  = ha,
    cluster_rows    = row_dend,
    cluster_columns = col_dend,
    show_row_names    = TRUE,
    show_column_names = FALSE,
    row_names_side  = "right",
    row_names_gp    = gpar(fontsize = 8, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp         = gpar(fontsize = 10, fontface = "bold"),
      labels_gp        = gpar(fontsize = 10, fontface = "bold"),
      legend_direction = "vertical",
      legend_height    = unit(4, "cm"),
      border           = "black"
    ),
    column_title    = paste0("Top ", nrow(heatmap_matrix_scaled),
                              " DEGs — ACC vs ACA (Training Set)"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    border          = TRUE
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right",
    merge_legend           = FALSE
  )

  dev.off()
  cat("✅ Heatmap saved\n")

} else {
  cat("⚠️  No significant DEGs found for heatmap at current thresholds.\n")
  cat("   Consider relaxing: PVAL_THRESH or FC_THRESH\n")
}

# =====================================================
# 13. PREPARE GENES FOR ENRICHMENT
# =====================================================

cat("\n========== Preparing Genes for Enrichment ==========\n")

sig_genes_for_enrich <- res %>%
  filter(adj.P.Val < PVAL_THRESH, abs(logFC) > FC_THRESH)

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys     = as.character(sig_genes_for_enrich$gene),
  keytype  = "SYMBOL",
  column   = "ENTREZID",
  multiVals = "first"
)
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
entrez_ids <- unique(as.character(entrez_ids))

cat("Significant DEGs for enrichment:", nrow(sig_genes_for_enrich), "\n")
cat("Mapped to Entrez IDs           :", length(entrez_ids), "\n")

if (length(entrez_ids) < 10) {
  cat("⚠️  Very few genes for enrichment (<10).",
      "Enrichment results may be unreliable.\n")
}

# =====================================================
# 14. GO ENRICHMENT ANALYSIS
# =====================================================

cat("\n========== GO Enrichment Analysis ==========\n")

ego_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  write.csv(as.data.frame(ego_bp),
            file.path(base_dir, "tables", "GO_BP_results.csv"),
            row.names = FALSE)
  cat("✅ GO BP:", nrow(as.data.frame(ego_bp)), "enriched terms\n")
}
if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
  write.csv(as.data.frame(ego_cc),
            file.path(base_dir, "tables", "GO_CC_results.csv"),
            row.names = FALSE)
  cat("✅ GO CC:", nrow(as.data.frame(ego_cc)), "enriched terms\n")
}
if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
  write.csv(as.data.frame(ego_mf),
            file.path(base_dir, "tables", "GO_MF_results.csv"),
            row.names = FALSE)
  cat("✅ GO MF:", nrow(as.data.frame(ego_mf)), "enriched terms\n")
}

# =====================================================
# 15. KEGG ENRICHMENT ANALYSIS
# =====================================================

cat("\n========== KEGG Enrichment Analysis ==========\n")

ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.csv(as.data.frame(ekegg),
            file.path(base_dir, "tables", "KEGG_results.csv"),
            row.names = FALSE)
  cat("✅ KEGG:", nrow(as.data.frame(ekegg)), "enriched pathways\n")
}

# =====================================================
# 16. SHARED DOT PLOT HELPERS
# =====================================================

wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  }, USE.NAMES = FALSE)
}

dotplot_theme <- theme_bw(base_size = 10) +
  theme(
    axis.text.x        = element_text(size = 9,  face = "bold"),
    axis.text.y        = element_text(size = 8,  face = "bold"),
    axis.title         = element_text(size = 11, face = "bold"),
    legend.title       = element_text(size = 9,  face = "bold"),
    legend.text        = element_text(size = 8,  face = "bold"),
    strip.text         = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.background  = element_blank(),
    legend.key         = element_blank()
  )

dotplot_color <- scale_color_gradientn(
  colours = c("#7B2D8B", "#C0308C", "#E8637A", "#F4A460", "#FDE8D8"),
  name    = "p.adjust"
)

dotplot_size <- scale_size_continuous(
  name  = "GeneCount",
  range = c(2, 12)
)

# =====================================================
# 17. GO BP DOT PLOT
# =====================================================

cat("\n========== GO Dot Plots ==========\n")

if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  bp_df <- as.data.frame(ego_bp) %>%
    arrange(p.adjust) %>% slice_head(n = 20) %>%
    mutate(Description = wrap_text(Description, width = 50)) %>%
    arrange(Count) %>%
    mutate(Description = factor(Description, levels = Description))

  p_go_bp <- ggplot(bp_df, aes(x = Count, y = Description,
                                 color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    dotplot_color + dotplot_size + dotplot_theme +
    labs(title = "GO Biological Process — ACC vs ACA",
         x = "Gene Count", y = NULL)

  ggsave(file.path(base_dir, "plots", "GO_BP_dotplot.tif"),
         p_go_bp, width = 10, height = 9, dpi = 600)
  cat("✅ GO BP dot plot saved\n")
}

# =====================================================
# 18. GO CC DOT PLOT
# =====================================================

if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
  cc_df <- as.data.frame(ego_cc) %>%
    arrange(p.adjust) %>% slice_head(n = 20) %>%
    mutate(Description = wrap_text(Description, width = 50)) %>%
    arrange(Count) %>%
    mutate(Description = factor(Description, levels = Description))

  p_go_cc <- ggplot(cc_df, aes(x = Count, y = Description,
                                 color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    dotplot_color + dotplot_size + dotplot_theme +
    labs(title = "GO Cellular Component — ACC vs ACA",
         x = "Gene Count", y = NULL)

  ggsave(file.path(base_dir, "plots", "GO_CC_dotplot.tif"),
         p_go_cc, width = 10, height = 7, dpi = 600)
  cat("✅ GO CC dot plot saved\n")
}

# =====================================================
# 19. GO MF DOT PLOT
# =====================================================

if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
  mf_df <- as.data.frame(ego_mf) %>%
    arrange(p.adjust) %>% slice_head(n = 20) %>%
    mutate(Description = wrap_text(Description, width = 50)) %>%
    arrange(Count) %>%
    mutate(Description = factor(Description, levels = Description))

  p_go_mf <- ggplot(mf_df, aes(x = Count, y = Description,
                                 color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    dotplot_color + dotplot_size + dotplot_theme +
    labs(title = "GO Molecular Function — ACC vs ACA",
         x = "Gene Count", y = NULL)

  ggsave(file.path(base_dir, "plots", "GO_MF_dotplot.tif"),
         p_go_mf, width = 10, height = 7, dpi = 600)
  cat("✅ GO MF dot plot saved\n")
}

# =====================================================
# 20. KEGG DOT PLOT
# =====================================================

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  kegg_df <- as.data.frame(ekegg) %>%
    arrange(p.adjust) %>% slice_head(n = 20) %>%
    mutate(Description = wrap_text(Description, width = 50)) %>%
    arrange(Count) %>%
    mutate(Description = factor(Description, levels = Description))

  p_kegg <- ggplot(kegg_df, aes(x = Count, y = Description,
                                  color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    dotplot_color + dotplot_size + dotplot_theme +
    labs(title = "KEGG Pathway Enrichment — ACC vs ACA",
         x = "Gene Count", y = NULL)

  ggsave(file.path(base_dir, "plots", "KEGG_dotplot.tif"),
         p_kegg, width = 10, height = 9, dpi = 600)
  cat("✅ KEGG dot plot saved\n")
}

# =====================================================
# 21. FINAL SUMMARY
# =====================================================

cat("\n\n")
cat("╔════════════════════════════════════════════════════════╗\n")
cat("║          DEG ANALYSIS COMPLETE — ACC vs ACA            ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")

cat("INPUT:\n")
cat("   Expression file :", expr_file,     "\n")
cat("   Metadata file   :", metadata_file, "\n")
cat("   Genes analysed  :", nrow(expr_data), "\n")
cat("   ACA (training)  :", sum(metadata$condition == "ACA"), "samples\n")
cat("   ACC (training)  :", sum(metadata$condition == "ACC"), "samples\n\n")

cat("THRESHOLDS:\n")
cat("   adj.P.Val <", PVAL_THRESH, "  |logFC| >", FC_THRESH, "\n\n")

cat("RESULTS:\n")
cat("   Total DEGs              :", nrow(res_sig),  "\n")
cat("   Upregulated in ACC      :", nrow(res_up),   "\n")
cat("   Downregulated in ACC    :", nrow(res_down), "\n")
cat("   (i.e. upregulated in ACA:", nrow(res_down), ")\n\n")

cat("OUTPUT DIRECTORY:", base_dir, "\n")
cat("   tables/ → 4 DEG tables + GO/KEGG enrichment tables\n")
cat("   plots/  → PCA, MA, Volcano, Heatmap, GO/KEGG dot plots\n\n")

cat("NEXT STEP:\n")
cat("   Validate top DEGs on the 30% test set:\n")
cat("   → expr_test_30.csv + metadata_test_30.csv\n")
cat("══════════════════════════════════════════════════════════\n")
