# ================================================================
# METABOLOMICS HEATMAPS — MTBLS8713 Targeted LC-MS
# Up & Down regulated metabolites — IBD vs nonIBD
#
# Outputs:
#   figures/heatmaps/MTBLS8713_metabolomics_heatmap.pdf
#   figures/heatmaps/summary/MTBLS8713_summary_heatmap.pdf
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
# COLOUR PALETTES
# ----------------------------------------------------------------
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

anno_colors <- list(
  Diagnosis = c(IBD = "#D73027", nonIBD = "#4575B4"),
  Direction = c("Up in IBD" = "#E31A1C", "Down in IBD" = "#1F78B4")
)

# ----------------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------------
cat("Loading MTBLS8713 processed data...\n")
rds      <- readRDS(file.path(BASE,
             "processed/MTBLS8713_metabolomics/MTBLS8713_processed.rds"))
mat      <- rds$normalized_matrix
col_data <- rds$col_data
sig      <- rds$sig_features
de       <- rds$de_results

cat("Matrix:", nrow(mat), "metabolites x", ncol(mat), "samples\n")
cat("Significant metabolites:", nrow(sig), "\n")

# ----------------------------------------------------------------
# SELECT FEATURES
# Top 40 up + all 21 down (show all downregulated since only 21)
# ----------------------------------------------------------------
up_mets   <- sig %>% filter(log2FC > 0) %>%
             arrange(desc(log2FC)) %>%
             slice_head(n = 40) %>%
             pull(metabolite)

down_mets <- sig %>% filter(log2FC < 0) %>%
             arrange(log2FC) %>%
             pull(metabolite)

all_mets  <- c(up_mets, down_mets)
cat("Showing:", length(up_mets), "up +", length(down_mets),
    "down =", length(all_mets), "metabolites\n")

# Clean up metabolite names for display (replace _ with space)
clean_names <- gsub("_", " ", all_mets)

# Build direction annotation
dir_vec <- setNames(
  ifelse(all_mets %in% up_mets, "Up in IBD", "Down in IBD"),
  all_mets
)

# ----------------------------------------------------------------
# SUBSET MATRIX AND RENAME ROWS
# ----------------------------------------------------------------
avail_mets <- all_mets[all_mets %in% rownames(mat)]
mat_sub    <- mat[avail_mets, , drop = FALSE]
rownames(mat_sub) <- gsub("_", " ", rownames(mat_sub))
names(dir_vec)    <- gsub("_", " ", names(dir_vec))

# Order samples: IBD first, nonIBD second
col_data$Diagnosis <- col_data$diagnosis_binary
sample_order       <- order(col_data$Diagnosis)
mat_sub            <- mat_sub[, sample_order, drop = FALSE]
col_data           <- col_data[sample_order, ]
rownames(col_data) <- col_data$sample_id

# Z-score scale
mat_z <- t(scale(t(mat_sub)))
mat_z[mat_z >  3] <-  3
mat_z[mat_z < -3] <- -3

# ----------------------------------------------------------------
# ANNOTATION DATA FRAMES
# ----------------------------------------------------------------
anno_col <- data.frame(
  Diagnosis = col_data$Diagnosis,
  row.names = rownames(col_data)
)

anno_row <- data.frame(
  Direction = dir_vec[rownames(mat_z)],
  row.names = rownames(mat_z)
)

# ================================================================
# HEATMAP 1 — Per-patient (all 179 samples)
# ================================================================
cat("\nGenerating per-patient heatmap...\n")

pdf(file.path(FIG_DIR, "MTBLS8713_metabolomics_heatmap.pdf"),
    width = 14, height = 12)
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
  main              = "MTBLS8713 — Targeted LC-MS Metabolomics\nTop 40 Up + 21 Down Regulated Metabolites (IBD vs nonIBD)",
  border_color      = NA,
  treeheight_row    = 30,
  treeheight_col    = 0
)
dev.off()
cat("  Saved: MTBLS8713_metabolomics_heatmap.pdf\n")

# ================================================================
# HEATMAP 2 — Summary (IBD average vs nonIBD average)
# ================================================================
cat("\nGenerating summary heatmap (group averages)...\n")

ibd_s  <- col_data$sample_id[col_data$Diagnosis == "IBD"]
ctrl_s <- col_data$sample_id[col_data$Diagnosis == "nonIBD"]

mean_ibd  <- rowMeans(mat_sub[, ibd_s,  drop = FALSE])
mean_ctrl <- rowMeans(mat_sub[, ctrl_s, drop = FALSE])

summary_mat        <- cbind(IBD = mean_ibd, nonIBD = mean_ctrl)
summary_mat        <- t(scale(t(summary_mat)))
summary_mat[summary_mat >  2] <-  2
summary_mat[summary_mat < -2] <- -2

# Add direction arrows to row names
lfc_vec <- mean_ibd - mean_ctrl
rownames(summary_mat) <- paste0(
  ifelse(lfc_vec[rownames(summary_mat)] > 0, "\u25b2 ", "\u25bc "),
  rownames(summary_mat)
)

anno_col_sum <- data.frame(
  Group     = c("IBD", "nonIBD"),
  row.names = c("IBD", "nonIBD")
)
anno_colors_sum <- list(
  Group     = c(IBD = "#D73027", nonIBD = "#4575B4"),
  Direction = c("Up in IBD" = "#E31A1C", "Down in IBD" = "#1F78B4")
)

# Row annotation for summary
anno_row_sum <- data.frame(
  Direction = dir_vec[gsub("^\u25b2 |\u25bc ", "", rownames(summary_mat))],
  row.names = rownames(summary_mat)
)

pdf(file.path(SUM_DIR, "MTBLS8713_summary_heatmap.pdf"),
    width = 6, height = 13)
pheatmap(
  summary_mat,
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
  main              = "MTBLS8713 — Targeted LC-MS Metabolomics\nGroup Average: IBD vs nonIBD\n(\u25b2 = higher in IBD, \u25bc = lower in IBD)",
  border_color      = "grey90",
  treeheight_row    = 30,
  treeheight_col    = 0,
  cellwidth         = 40,
  angle_col         = 0
)
dev.off()
cat("  Saved: MTBLS8713_summary_heatmap.pdf\n")

# ================================================================
# SUMMARY
# ================================================================
cat("\n=== DONE ===\n")
cat("Metabolites shown:", length(avail_mets), "\n")
cat("  Up in IBD:", length(up_mets), "\n")
cat("  Down in IBD:", length(down_mets), "\n")
cat("\nTop 5 upregulated metabolites in IBD:\n")
print(sig %>% filter(log2FC > 0) %>%
        arrange(desc(log2FC)) %>%
        select(metabolite, log2FC, padj) %>%
        head(5) %>%
        mutate(metabolite = gsub("_"," ",metabolite),
               log2FC = round(log2FC,3),
               padj = signif(padj,3)))
cat("\nTop 5 downregulated metabolites in IBD:\n")
print(sig %>% filter(log2FC < 0) %>%
        arrange(log2FC) %>%
        select(metabolite, log2FC, padj) %>%
        head(5) %>%
        mutate(metabolite = gsub("_"," ",metabolite),
               log2FC = round(log2FC,3),
               padj = signif(padj,3)))
