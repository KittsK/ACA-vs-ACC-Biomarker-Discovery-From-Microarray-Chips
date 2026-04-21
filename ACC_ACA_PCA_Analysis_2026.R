################################################################################
# ACC PROJECT — PCA PLOT (Before & After Batch Correction)
# Datasets  : GSE12368 & GSE10927
# Comparison: ACA vs ACC
#
# FIXED:
#   - ACA = Circle  (shape 21) filled BLUE  (#2166AC)
#   - ACC = Triangle (shape 24) filled RED   (#D6604D)
#   - Dataset (batch) shown via OUTLINE colour only
#     GSE12368 = Green (#1B7837) | GSE10927 = Orange (#E08214)
#
# PREREQUISITES:
#   The following objects must exist in your R session:
#     merged_mat     — merged expression matrix (before ComBat)
#     expr_corrected — batch-corrected expression matrix (after ComBat)
#     meta_ordered   — metadata with sample_id, condition, batch columns
#
# OUTPUT (saved to D:/ACC Project/results/):
#   PCA_before_after_batch.tiff  (700 dpi)
#   PCA_before_after_batch.pdf
#   PCA_before_after_batch.png
################################################################################

# ── Packages ──────────────────────────────────────────────────────────────────
pkgs <- c("ggplot2", "patchwork", "ggrepel", "scales", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(dplyr)
})

# ── Output directory ──────────────────────────────────────────────────────────
OUT_DIR <- "D:/ACC Project/results"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

RANDOM_SEED <- 42

# ==============================================================================
# COLOUR & SHAPE DEFINITIONS
# ==============================================================================
#
# CONDITION (fill of point):
#   ACA → Circle   (shape 21) → Blue  #2166AC
#   ACC → Triangle (shape 24) → Red   #D6604D
#
# DATASET (outline/colour of point):
#   GSE12368 → Green  #1B7837
#   GSE10927 → Orange #E08214
#
# Using filled shapes (21, 24) allows fill= and colour= to be
# mapped INDEPENDENTLY — fill = condition, colour = dataset.
# This cleanly separates biological group from technical batch
# in a single point without needing two separate glyphs.

COND_FILLS  <- c(ACA = "#2166AC", ACC = "#D6604D")
COND_SHAPES <- c(ACA = 21,        ACC = 24)
BATCH_COLS  <- c(GSE12368 = "#1B7837", GSE10927 = "#E08214")


# ==============================================================================
# PCA FUNCTION
# ==============================================================================

#' Build one PCA panel
#'
#' @param mat          numeric matrix (genes × samples)
#' @param meta         data.frame with sample_id, condition, batch columns
#' @param title_label  string — panel title
#' @param seed         integer — random seed for reproducibility

build_pca_panel <- function(mat, meta, title_label, seed = RANDOM_SEED) {

  set.seed(seed)

  # ── Compute PCA ─────────────────────────────────────────────────────────────
  pca    <- prcomp(t(mat), scale. = TRUE, center = TRUE)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$sample_id <- rownames(scores)

  # Merge with metadata
  scores <- merge(
    scores,
    meta[, c("sample_id", "condition", "batch")],
    by = "sample_id"
  )

  # Variance explained
  var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  # Short label — last 3 digits of GSM number
  scores$label <- gsub("GSM(\\d+)", "\\1", scores$sample_id)
  scores$label <- substr(scores$label,
                          nchar(scores$label) - 2,
                          nchar(scores$label))

  # ── Plot ─────────────────────────────────────────────────────────────────────
  ggplot(scores, aes(x = PC1, y = PC2)) +

    # ── Dataset ellipses (dashed border, no fill) ──────────────────────────
    # Shows where each dataset clusters — reveals batch effect before
    # ComBat and confirms mixing after ComBat
    stat_ellipse(
      aes(group = batch, colour = batch),
      type        = "norm",
      level       = 0.95,
      linewidth   = 0.85,
      linetype    = "dashed",
      show.legend = FALSE
    ) +

    # ── Condition ellipses (solid, semi-transparent fill) ──────────────────
    # Shows ACA vs ACC biological separation
    stat_ellipse(
      aes(group = condition, fill = condition),
      geom        = "polygon",
      type        = "norm",
      level       = 0.95,
      alpha       = 0.10,
      colour      = NA,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(group = condition, colour = condition),
      type        = "norm",
      level       = 0.95,
      linewidth   = 0.9,
      linetype    = "solid",
      show.legend = FALSE
    ) +

    # ── Reference lines at origin ──────────────────────────────────────────
    geom_hline(yintercept = 0, colour = "grey82",
               linewidth = 0.3, linetype = "dotted") +
    geom_vline(xintercept = 0, colour = "grey82",
               linewidth = 0.3, linetype = "dotted") +

    # ── Points ────────────────────────────────────────────────────────────
    # Layer 1: white background ring — makes outline colour visible on any bg
    geom_point(
      aes(shape = condition),
      fill   = "white",
      colour = "white",
      size   = 4.5,
      stroke = 1.8,
      show.legend = FALSE
    ) +
    # Layer 2: dataset colour outline (stroke)
    geom_point(
      aes(shape = condition, colour = batch),
      fill   = "white",
      size   = 4.5,
      stroke = 1.5,
      show.legend = TRUE
    ) +
    # Layer 3: condition fill (inside the shape)
    # This is the KEY layer — fill= drives the interior colour
    # shape= drives circle vs triangle
    # colour= (from layer 2) drives the visible outline colour
    geom_point(
      aes(shape = condition, fill = condition),
      colour = NA,
      size   = 4.5,
      stroke = 0,
      alpha  = 0.85,
      show.legend = TRUE
    ) +

    # ── Sample labels ─────────────────────────────────────────────────────
    geom_text_repel(
      aes(label = label, colour = batch),
      size               = 2.5,
      fontface           = "plain",
      box.padding        = 0.35,
      point.padding      = 0.25,
      force              = 1.2,
      max.overlaps       = 40,
      segment.size       = 0.25,
      segment.alpha      = 0.6,
      min.segment.length = 0.15,
      show.legend        = FALSE
    ) +

    # ── Scales ────────────────────────────────────────────────────────────
    # fill  → condition interior colour (ACA=blue, ACC=red)
    # colour → combined scale: condition solid ellipse + dataset outline
    # shape → condition shape (ACA=circle, ACC=triangle)
    scale_fill_manual(
      values = COND_FILLS,
      name   = "Diagnosis",
      guide  = guide_legend(
        override.aes = list(
          shape  = c(21, 24),
          fill   = COND_FILLS,
          colour = "grey30",
          size   = 4.5,
          stroke = 0.8
        )
      )
    ) +
    scale_colour_manual(
      values = c(COND_FILLS, BATCH_COLS),
      breaks = names(BATCH_COLS),          # only show dataset in colour legend
      name   = "Dataset",
      guide  = guide_legend(
        override.aes = list(
          shape  = 21,
          fill   = "white",
          colour = BATCH_COLS,
          size   = 4.5,
          stroke = 1.5
        )
      )
    ) +
    scale_shape_manual(
      values = COND_SHAPES,
      guide  = "none"                      # shape already shown in fill legend
    ) +

    # ── Labels ────────────────────────────────────────────────────────────
    labs(
      title = title_label,
      x     = paste0("PC1: ", var_exp[1], "% variance"),
      y     = paste0("PC2: ", var_exp[2], "% variance")
    ) +

    # ── Theme ─────────────────────────────────────────────────────────────
    theme_classic(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 14,
                                       colour = "#1C2833"),
      axis.title       = element_text(face = "bold", size = 11),
      axis.text        = element_text(size = 9, colour = "grey30"),
      panel.border     = element_rect(colour = "black", fill = NA,
                                       linewidth = 1.2),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.35,
                                       linetype = "dashed"),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.key.size  = unit(0.55, "cm"),
      legend.background = element_rect(fill = "white", colour = "grey85",
                                        linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(10, 12, 6, 12)
    )
}


# ==============================================================================
# BUILD PANELS
# ==============================================================================

cat("Building PCA panels...\n")

p_before <- build_pca_panel(
  mat         = merged_mat,
  meta        = meta_ordered,
  title_label = "Before Batch Correction"
)

p_after <- build_pca_panel(
  mat         = expr_corrected,
  meta        = meta_ordered,
  title_label = "After Batch Correction"
)

cat("PCA panels built.\n")


# ==============================================================================
# ASSEMBLE COMBINED FIGURE
# ==============================================================================

combined_pca <- (p_before | p_after) +
  plot_annotation(
    title    = "GSE12368 & GSE10927 \u2014 Batch Correction QC",
    subtitle = paste0(
      "ACA (adenoma, n=", sum(meta_ordered$condition == "ACA"), ")  ",
      "\u25CF Blue  \u2022  ",
      "ACC (carcinoma, n=", sum(meta_ordered$condition == "ACC"), ")  ",
      "\u25B2 Red  \u2022  ",
      "Ellipses: 95% CI  \u2022  ",
      "Dashed = dataset boundary  |  Solid = diagnosis boundary"
    ),
    caption = paste0(
      "Platform: GPL570 (Affymetrix HG-U133 Plus 2.0)  \u2022  ",
      "ComBat batch correction (sva)  \u2022  ",
      "Point fill = Diagnosis  |  Point outline = Dataset  \u2022  ",
      "Generated: ", format(Sys.Date(), "%Y-%m-%d")
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5,
                                    colour = "#1C2833",
                                    margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey38",
                                    margin = margin(b = 6)),
      plot.caption  = element_text(size = 8,  hjust = 0.5, colour = "grey55",
                                    margin = margin(t = 8)),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin     = margin(14, 14, 10, 14)
    )
  )


# ==============================================================================
# SAVE
# ==============================================================================

prefix <- file.path(OUT_DIR, "PCA_before_after_batch")

# TIFF 700 dpi
ggsave(paste0(prefix, ".tiff"),
       plot        = combined_pca,
       width       = 16, height = 8,
       dpi         = 700,
       device      = "tiff",
       compression = "lzw")
cat("Saved:", paste0(prefix, ".tiff\n"))

# PDF
ggsave(paste0(prefix, ".pdf"),
       plot   = combined_pca,
       width  = 16, height = 8,
       device = cairo_pdf)
cat("Saved:", paste0(prefix, ".pdf\n"))

# PNG 700 dpi
ggsave(paste0(prefix, ".png"),
       plot  = combined_pca,
       width = 16, height = 8,
       dpi   = 700)
cat("Saved:", paste0(prefix, ".png\n"))

cat("\n============================== DONE ==============================\n")
cat("ACA = Circle  (shape 21) filled Blue  | outline = dataset colour\n")
cat("ACC = Triangle(shape 24) filled Red   | outline = dataset colour\n")
cat("GSE12368 outline = Green  (#1B7837)\n")
cat("GSE10927 outline = Orange (#E08214)\n")
cat("==================================================================\n")
