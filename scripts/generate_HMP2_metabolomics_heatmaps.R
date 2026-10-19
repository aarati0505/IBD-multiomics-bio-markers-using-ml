# ================================================================
# HEATMAPS — HMP2 Metabolomics (untargeted LC-MS, 546 samples)
# 8,419 significant features — show top 40 up + 40 down
#
# Feature IDs are internal mass spec codes (C18n_QI*, C8p_QI*)
# not named metabolites — this is normal for untargeted data
#
# Outputs:
#   figures/heatmaps/HMP2_metabolomics_heatmap.pdf
#   figures/heatmaps/summary/HMP2_metabolomics_summary_heatmap.pdf
# ================================================================

library(pheatmap)
library(RColorBrewer)
library(dplyr)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
FIG_DIR <- file.path(BASE, "figures", "heatmaps")
SUM_DIR <- file.path(BASE, "figures", "heatmaps", "summary")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUM_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------
# COLOUR SETUP
# ----------------------------------------------------------------
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

anno_colors <- list(
  Diagnosis = c(IBD    = "#D73027", nonIBD = "#4575B4"),
  Direction = c("Up in IBD"   = "#E31A1C",
                "Down in IBD" = "#1F78B4")
)

# ----------------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------------
cat("Loading HMP2 metabolomics processed data...\n")
rds      <- readRDS(file.path(BASE,
             "processed/HMP2_metabolomics/HMP2_metabolomics_processed.rds"))
mat      <- rds$normalized_matrix
col_data <- rds$col_data
sig      <- rds$sig_features
de       <- rds$de_results

cat("Matrix:", nrow(mat), "features x", ncol(mat), "samples\n")
cat("Significant features:", nrow(sig), "\n")

# ----------------------------------------------------------------
# SELECT TOP 40 UP + TOP 40 DOWN by |LFC|
# ----------------------------------------------------------------
up_feats   <- sig %>% filter(log2FC > 0) %>%
              arrange(desc(log2FC)) %>%
              slice_head(n = 40) %>% pull(feature)

down_feats <- sig %>% filter(log2FC < 0) %>%
              arrange(log2FC) %>%
              slice_head(n = 40) %>% pull(feature)

sel_feats  <- c(up_feats, down_feats)
cat("Showing:", length(up_feats), "up +", length(down_feats),
    "down =", length(sel_feats), "features\n")

# Direction annotation vector
dir_vec <- setNames(
  ifelse(sel_feats %in% up_feats, "Up in IBD", "Down in IBD"),
  sel_feats
)

# ----------------------------------------------------------------
# SUBSET AND ORDER MATRIX
# ----------------------------------------------------------------
avail     <- sel_feats[sel_feats %in% rownames(mat)]
mat_sub   <- mat[avail, , drop = FALSE]

col_data$Diagnosis <- col_data$diagnosis_binary
# Sort: IBD first, then nonIBD
samp_order         <- order(col_data$Diagnosis)
mat_sub            <- mat_sub[, samp_order, drop = FALSE]
col_data           <- col_data[samp_order, ]
rownames(col_data) <- col_data$sample_id

# Z-score scale rows
mat_z <- t(scale(t(mat_sub)))
mat_z[mat_z >  3] <-  3
mat_z[mat_z < -3] <- -3

# Column annotation
anno_col <- data.frame(
  Diagnosis = col_data$Diagnosis,
  row.names = col_data$sample_id
)

# Row annotation
anno_row <- data.frame(
  Direction = dir_vec[rownames(mat_z)],
  row.names = rownames(mat_z)
)

# ================================================================
# HEATMAP 1 — Per-patient (all 546 samples)
# ================================================================
cat("\nGenerating per-patient heatmap...\n")

pdf(file.path(FIG_DIR, "HMP2_metabolomics_heatmap.pdf"),
    width = 16, height = 12)
pheatmap(
  mat_z,
  color             = heatmap_colors,
  annotation_col    = anno_col,
  annotation_row    = anno_row,
  annotation_colors = anno_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 7,
  fontsize          = 9,
  main              = paste0(
    "HMP2 Stool Metabolomics — Untargeted LC-MS\n",
    "Top 40 Up + 40 Down Regulated Features (IBD vs nonIBD)\n",
    "(Feature IDs are internal mass-spec codes — unnamed metabolites)"
  ),
  border_color      = NA,
  treeheight_row    = 30,
  treeheight_col    = 0
)
dev.off()
cat("  Saved: HMP2_metabolomics_heatmap.pdf\n")

# ================================================================
# HEATMAP 2 — Summary (group averages)
# ================================================================
cat("\nGenerating summary heatmap...\n")

ibd_s  <- col_data$sample_id[col_data$Diagnosis == "IBD"]
ctrl_s <- col_data$sample_id[col_data$Diagnosis == "nonIBD"]

mean_ibd  <- rowMeans(mat_sub[, ibd_s,  drop = FALSE])
mean_ctrl <- rowMeans(mat_sub[, ctrl_s, drop = FALSE])

sum_mat <- cbind(IBD = mean_ibd, nonIBD = mean_ctrl)
sum_mat <- t(scale(t(sum_mat)))
sum_mat[sum_mat >  2] <-  2
sum_mat[sum_mat < -2] <- -2

# Add direction arrows
lfc_vec <- mean_ibd - mean_ctrl
rownames(sum_mat) <- paste0(
  ifelse(lfc_vec[rownames(sum_mat)] > 0, "\u25b2 ", "\u25bc "),
  rownames(sum_mat)
)

# Row annotation (strip arrows for matching)
clean_ids <- gsub("^\u25b2 |\u25bc ", "", rownames(sum_mat))
anno_row_sum <- data.frame(
  Direction = dir_vec[clean_ids],
  row.names = rownames(sum_mat)
)

anno_col_sum <- data.frame(
  Group     = c("IBD", "nonIBD"),
  row.names = c("IBD", "nonIBD")
)
anno_colors_sum <- list(
  Group     = c(IBD = "#D73027", nonIBD = "#4575B4"),
  Direction = c("Up in IBD" = "#E31A1C", "Down in IBD" = "#1F78B4")
)

pdf(file.path(SUM_DIR, "HMP2_metabolomics_summary_heatmap.pdf"),
    width = 6, height = 13)
pheatmap(
  sum_mat,
  color             = heatmap_colors,
  annotation_col    = anno_col_sum,
  annotation_row    = anno_row_sum,
  annotation_colors = anno_colors_sum,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 7,
  fontsize_col      = 11,
  fontsize          = 9,
  main              = paste0(
    "HMP2 Stool Metabolomics\n",
    "Group Average: IBD vs nonIBD\n",
    "(\u25b2 = higher in IBD, \u25bc = lower in IBD)"
  ),
  border_color      = "grey90",
  treeheight_row    = 30,
  treeheight_col    = 0,
  cellwidth         = 40,
  angle_col         = 0
)
dev.off()
cat("  Saved: HMP2_metabolomics_summary_heatmap.pdf\n")

# ================================================================
# SUMMARY STATS
# ================================================================
cat("\n=== DONE ===\n")
cat("Total sig features:", nrow(sig), "\n")
cat("Up in IBD:", sum(sig$log2FC > 0),
    "| Down in IBD:", sum(sig$log2FC < 0), "\n")
cat("\nTop 5 most upregulated features (largest LFC):\n")
print(sig %>% filter(log2FC > 0) %>% arrange(desc(log2FC)) %>%
        select(feature, log2FC, padj) %>% head(5) %>%
        mutate(log2FC = round(log2FC, 3), padj = signif(padj, 3)))
cat("\nTop 5 most downregulated features:\n")
print(sig %>% filter(log2FC < 0) %>% arrange(log2FC) %>%
        select(feature, log2FC, padj) %>% head(5) %>%
        mutate(log2FC = round(log2FC, 3), padj = signif(padj, 3)))
