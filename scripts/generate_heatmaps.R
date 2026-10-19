# ================================================================
# HEATMAPS — Up & Down Regulated Genes and Proteins
# IBD Multi-Omics Project
#
# Datasets covered:
#   1. GSE137344  — Ileum transcriptomics (RNA-seq)
#   2. GSE186507  — Blood transcriptomics (microarray)
#   3. GSE193677  — Biopsy transcriptomics (microarray)
#   4. PXD001608  — Colon tissue proteomics (MaxQuant)
#   5. HMP2 ECS   — Stool metaproteomics (enzyme functions)
#   6. HMP2 KOS   — Stool metaproteomics (KEGG orthologs)
#
# Output: figures/heatmaps/*.pdf (one per dataset)
# ================================================================

library(pheatmap)
library(RColorBrewer)
library(dplyr)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
FIG_DIR <- file.path(BASE, "figures", "heatmaps")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------
# COLOUR PALETTE — red (up in IBD) to blue (down in IBD)
# ----------------------------------------------------------------
heatmap_colors <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

# Annotation colours for IBD / nonIBD bar
anno_colors <- list(
  Diagnosis = c(IBD = "#D73027", nonIBD = "#4575B4")
)

# ----------------------------------------------------------------
# HELPER: draw one heatmap and save as PDF
#   direction_vec: named character vector mapping feature -> "Up in IBD"
#                  or "Down in IBD". If NULL no row annotation is drawn.
#   subgroup_col:  name of a column in col_data to use as a second
#                  column annotation bar (e.g. "SubDiagnosis" for CD/UC split)
# ----------------------------------------------------------------
draw_heatmap <- function(mat,
                         col_data,
                         title,
                         filename,
                         direction_vec  = NULL,
                         subgroup_col   = NULL,
                         subgroup_colors = NULL,
                         show_rownames  = TRUE,
                         fontsize_row   = 7,
                         width          = 14,
                         height         = 10) {

  # Order samples: by subgroup if provided, else IBD then nonIBD
  if (!is.null(subgroup_col) && subgroup_col %in% names(col_data)) {
    sample_order <- order(col_data[[subgroup_col]])
  } else {
    sample_order <- order(col_data$Diagnosis)
  }
  mat      <- mat[, sample_order, drop = FALSE]
  col_data <- col_data[sample_order, , drop = FALSE]

  # Column annotation bar (top) — always includes Diagnosis
  anno_col_df        <- data.frame(Diagnosis = col_data$Diagnosis,
                                   row.names = rownames(col_data))
  col_anno_colors    <- anno_colors   # copy base colour list

  # Add subgroup bar if requested
  if (!is.null(subgroup_col) && subgroup_col %in% names(col_data)) {
    anno_col_df[[subgroup_col]] <- col_data[[subgroup_col]]
    if (!is.null(subgroup_colors)) {
      col_anno_colors[[subgroup_col]] <- subgroup_colors
    }
  }

  # Row annotation bar (left) — Direction
  anno_row   <- NULL

  if (!is.null(direction_vec)) {
    dir_aligned <- direction_vec[rownames(mat)]
    dir_aligned[is.na(dir_aligned)] <- "Unknown"
    anno_row <- data.frame(Direction = dir_aligned,
                           row.names = rownames(mat))
    col_anno_colors[["Direction"]] <- c(
      "Up in IBD"   = "#E31A1C",
      "Down in IBD" = "#1F78B4",
      "Unknown"     = "grey80"
    )
  }

  # Z-score scale rows
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[mat_scaled >  3] <-  3
  mat_scaled[mat_scaled < -3] <- -3

  pdf(filename, width = width, height = height)
  pheatmap(
    mat_scaled,
    color             = heatmap_colors,
    annotation_col    = anno_col_df,
    annotation_row    = anno_row,
    annotation_colors = col_anno_colors,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = show_rownames,
    show_colnames     = FALSE,
    fontsize_row      = fontsize_row,
    fontsize          = 9,
    main              = title,
    border_color      = NA,
    treeheight_row    = 30,
    treeheight_col    = 0
  )
  dev.off()
  cat("  Saved:", basename(filename), "\n")
}

# ================================================================
# 1. GSE137344 — Ileum Transcriptomics
#    165 significant DEGs, |log2FC| > 1, padj < 0.05
#    Show top 40 up + top 40 down by |LFC| for readability
# ================================================================
cat("\n--- GSE137344: Ileum Transcriptomics ---\n")

rds137  <- readRDS(file.path(BASE, "processed/GSE137344/GSE137344_normalized.rds"))
sig137  <- read.csv(file.path(BASE, "processed/GSE137344/GSE137344_significant_DEGs.csv"),
                    stringsAsFactors = FALSE)

cat("  Significant DEGs:", nrow(sig137), "\n")

# Top 40 up + top 40 down by |LFC|
up137    <- sig137 %>% filter(log2FoldChange > 0) %>%
            arrange(desc(log2FoldChange)) %>% slice_head(n = 40) %>% pull(gene)
down137  <- sig137 %>% filter(log2FoldChange < 0) %>%
            arrange(log2FoldChange)       %>% slice_head(n = 40) %>% pull(gene)
genes137 <- c(up137, down137)
cat("  Showing:", length(up137), "up +", length(down137), "down =", length(genes137), "genes\n")

norm_mat137 <- rds137$normalized_matrix
avail137    <- genes137[genes137 %in% rownames(norm_mat137)]
mat137      <- norm_mat137[avail137, , drop = FALSE]

cd137           <- rds137$col_data
cd137$Diagnosis <- ifelse(cd137$diagnosis_binary == "IBD", "IBD", "nonIBD")
cd137$SubDiagnosis <- factor(cd137$diagnosis,
                              levels = c("CD", "UC", "nonIBD"))
rownames(cd137) <- cd137$sample_id

# Build direction vector for row annotation
dir137 <- setNames(
  ifelse(genes137 %in% up137, "Up in IBD", "Down in IBD"),
  genes137
)

draw_heatmap(
  mat             = mat137,
  col_data        = cd137,
  direction_vec   = dir137,
  subgroup_col    = "SubDiagnosis",
  subgroup_colors = c(CD = "#FF7F00", UC = "#984EA3", nonIBD = "#4575B4"),
  title           = "GSE137344 — Ileum Transcriptomics\nTop 40 Up + 40 Down Regulated Genes (CD | UC | nonIBD)",
  filename        = file.path(FIG_DIR, "GSE137344_heatmap.pdf"),
  fontsize_row    = 6,
  height          = 12
)

# ================================================================
# 2. GSE186507 — Blood Transcriptomics
#    146 significant DEGs, |logFC| > 0.5, adj.P.Val < 0.05
#    Show top 40 up + top 40 down
# ================================================================
cat("\n--- GSE186507: Blood Transcriptomics ---\n")

rds186  <- readRDS(file.path(BASE, "processed/GSE186507/GSE186507_normalized.rds"))
deg186  <- read.csv(file.path(BASE, "processed/GSE186507/GSE186507_DEG_results.csv"),
                    stringsAsFactors = FALSE)

sig186 <- deg186 %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 0.5) %>%
  arrange(adj.P.Val)

cat("  Significant DEGs (|LFC|>0.5, padj<0.05):", nrow(sig186), "\n")

up186    <- sig186 %>% filter(logFC > 0) %>%
            arrange(desc(logFC)) %>% slice_head(n = 40) %>% pull(gene)
down186  <- sig186 %>% filter(logFC < 0) %>%
            arrange(logFC)       %>% slice_head(n = 40) %>% pull(gene)
genes186 <- c(up186, down186)
cat("  Showing:", length(up186), "up +", length(down186), "down =", length(genes186), "genes\n")

norm_mat186 <- rds186$normalized_matrix
avail186    <- genes186[genes186 %in% rownames(norm_mat186)]
mat186      <- norm_mat186[avail186, , drop = FALSE]

cd186           <- rds186$col_data
cd186$Diagnosis <- ifelse(cd186$diagnosis_binary == "IBD", "IBD", "nonIBD")
rownames(cd186) <- cd186$sample_id

dir186 <- setNames(
  ifelse(genes186 %in% up186, "Up in IBD", "Down in IBD"),
  genes186
)

draw_heatmap(
  mat           = mat186,
  col_data      = cd186,
  direction_vec = dir186,
  title         = "GSE186507 — Blood Transcriptomics\nTop 40 Up + 40 Down Regulated Genes (IBD vs nonIBD)",
  filename      = file.path(FIG_DIR, "GSE186507_heatmap.pdf"),
  fontsize_row  = 6,
  width         = 16,
  height        = 12
)

# ================================================================
# 3. GSE193677 — Biopsy Transcriptomics
#    LFC range very small (max ~0.45) — use top 50 by padj
#    Split into up/down by logFC sign
# ================================================================
cat("\n--- GSE193677: Biopsy Transcriptomics ---\n")

rds193  <- readRDS(file.path(BASE, "processed/GSE193677/GSE193677_normalized.rds"))
deg193  <- read.csv(file.path(BASE, "processed/GSE193677/GSE193677_DEG_results.csv"),
                    stringsAsFactors = FALSE)

sig193 <- deg193 %>%
  filter(!is.na(P.Value)) %>%
  arrange(P.Value)

cat("  Genes with padj < 0.05:", sum(!is.na(deg193$adj.P.Val) & deg193$adj.P.Val < 0.05), "\n")
cat("  Using top 50 by raw p-value for heatmap (LFC range too small for BH correction)\n")

# Top 25 up + top 25 down among top 200 most significant by p-value
top200   <- deg193 %>% filter(!is.na(P.Value)) %>% arrange(P.Value) %>% slice_head(n = 200)
up193    <- top200 %>% filter(logFC > 0) %>% slice_head(n = 25) %>% pull(gene)
down193  <- top200 %>% filter(logFC < 0) %>% slice_head(n = 25) %>% pull(gene)
genes193 <- c(up193, down193)
cat("  Showing:", length(up193), "up +", length(down193), "down =", length(genes193), "genes\n")

norm_mat193 <- rds193$normalized_matrix
avail193    <- genes193[genes193 %in% rownames(norm_mat193)]
mat193      <- norm_mat193[avail193, , drop = FALSE]

cd193           <- rds193$col_data
cd193$Diagnosis <- ifelse(cd193$diagnosis_binary == "IBD", "IBD", "nonIBD")
rownames(cd193) <- cd193$sample_id

dir193 <- setNames(
  ifelse(genes193 %in% up193, "Up in IBD", "Down in IBD"),
  genes193
)

draw_heatmap(
  mat           = mat193,
  col_data      = cd193,
  direction_vec = dir193,
  title         = "GSE193677 — Biopsy Transcriptomics\nTop 25 Up + 25 Down Regulated Genes by p-value (IBD vs nonIBD)",
  filename      = file.path(FIG_DIR, "GSE193677_heatmap.pdf"),
  fontsize_row  = 7,
  width         = 18,
  height        = 10
)

# ================================================================
# 4. PXD001608 — Colon Tissue Proteomics
#    721 significant proteins — show top 40 up + 40 down by |LFC|
# ================================================================
cat("\n--- PXD001608: Colon Tissue Proteomics ---\n")

pxd_rds   <- readRDS(file.path(BASE,
              "processed/proteomics_MaxQuant/proteomics_MaxQuant_processed.rds"))
sig_prot  <- read.csv(file.path(BASE,
              "processed/proteomics_MaxQuant/proteomics_MaxQuant_sig_proteins.csv"),
              stringsAsFactors = FALSE)

cat("  Significant proteins (|LFC|>1, padj<0.05):", nrow(sig_prot), "\n")

up_prot   <- sig_prot %>% filter(log2FC > 0) %>%
             arrange(desc(log2FC)) %>% slice_head(n = 40) %>% pull(protein)
down_prot <- sig_prot %>% filter(log2FC < 0) %>%
             arrange(log2FC)       %>% slice_head(n = 40) %>% pull(protein)
prots_sel <- c(up_prot, down_prot)
cat("  Showing:", length(up_prot), "up +", length(down_prot), "down =",
    length(prots_sel), "proteins\n")

norm_pxd  <- pxd_rds$normalized_matrix
avail_pxd <- prots_sel[prots_sel %in% rownames(norm_pxd)]
mat_pxd   <- norm_pxd[avail_pxd, , drop = FALSE]

cd_pxd           <- pxd_rds$col_data
cd_pxd$Diagnosis <- ifelse(cd_pxd$diagnosis_binary == "IBD", "IBD", "nonIBD")

dir_pxd <- setNames(
  ifelse(prots_sel %in% up_prot, "Up in IBD", "Down in IBD"),
  prots_sel
)

draw_heatmap(
  mat           = mat_pxd,
  col_data      = cd_pxd,
  direction_vec = dir_pxd,
  title         = "PXD001608 — Colon Tissue Proteomics\nTop 40 Up + 40 Down Regulated Proteins (UC vs Ctrl)",
  filename      = file.path(FIG_DIR, "PXD001608_proteomics_heatmap.pdf"),
  fontsize_row  = 6,
  height        = 12
)

# ================================================================
# 5. HMP2 Metaproteomics — ECS (Enzyme Functions)
#    75 significant features — show all
# ================================================================
cat("\n--- HMP2 ECS: Stool Metaproteomics (Enzyme Functions) ---\n")

ecs_rds  <- readRDS(file.path(BASE,
             "processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_processed.rds"))
sig_ecs  <- read.csv(file.path(BASE,
             "processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_sig_features.csv"),
             stringsAsFactors = FALSE)

cat("  Significant ECS features:", nrow(sig_ecs), "\n")

# Shorten long feature names for display (keep EC number + short name)
shorten_ecs <- function(x) {
  # Extract EC number part e.g. "1.1.1.1"
  ec <- regmatches(x, regexpr("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+", x))
  # Extract name part (after ": ")
  nm <- sub("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+: ", "", x)
  nm <- substr(nm, 1, 35)
  if (length(ec) > 0) paste0(ec, ": ", nm) else nm
}

up_ecs   <- sig_ecs %>% filter(log2FC > 0) %>%
            arrange(desc(log2FC)) %>% pull(feature)
down_ecs <- sig_ecs %>% filter(log2FC < 0) %>%
            arrange(log2FC)       %>% pull(feature)
feats_ecs <- c(up_ecs, down_ecs)

norm_ecs  <- ecs_rds$normalized_matrix
avail_ecs <- feats_ecs[feats_ecs %in% rownames(norm_ecs)]
mat_ecs   <- norm_ecs[avail_ecs, , drop = FALSE]
rownames(mat_ecs) <- sapply(rownames(mat_ecs), shorten_ecs)

cd_ecs           <- ecs_rds$col_data
cd_ecs$Diagnosis <- ifelse(cd_ecs$diagnosis_binary == "IBD", "IBD", "nonIBD")

dir_ecs <- setNames(
  ifelse(feats_ecs %in% up_ecs, "Up in IBD", "Down in IBD"),
  sapply(feats_ecs, shorten_ecs)
)

draw_heatmap(
  mat           = mat_ecs,
  col_data      = cd_ecs,
  direction_vec = dir_ecs,
  title         = "HMP2 Stool Metaproteomics — Enzyme Functions (ECS)\nSignificant Features Up & Down (IBD vs nonIBD)",
  filename      = file.path(FIG_DIR, "HMP2_ECS_heatmap.pdf"),
  fontsize_row  = 5.5,
  width         = 16,
  height        = 14
)

# ================================================================
# 6. HMP2 Metaproteomics — KOS (KEGG Orthologs)
#    81 significant features — show all
# ================================================================
cat("\n--- HMP2 KOS: Stool Metaproteomics (KEGG Orthologs) ---\n")

kos_rds  <- readRDS(file.path(BASE,
             "processed/HMP2_metaproteomics/HMP2_metaproteomics_kos_processed.rds"))
sig_kos  <- read.csv(file.path(BASE,
             "processed/HMP2_metaproteomics/HMP2_metaproteomics_kos_sig_features.csv"),
             stringsAsFactors = FALSE)

cat("  Significant KOS features:", nrow(sig_kos), "\n")

# Shorten KO names: keep KO ID + first 35 chars of name
shorten_ko <- function(x) {
  ko  <- regmatches(x, regexpr("^K[0-9]+", x))
  nm  <- sub("^K[0-9]+: ", "", x)
  # Remove EC numbers in brackets to save space
  nm  <- sub("\\s*\\[EC:[^\\]]+\\]", "", nm)
  nm  <- substr(nm, 1, 35)
  if (length(ko) > 0) paste0(ko, ": ", nm) else nm
}

up_kos   <- sig_kos %>% filter(log2FC > 0) %>%
            arrange(desc(log2FC)) %>% pull(feature)
down_kos <- sig_kos %>% filter(log2FC < 0) %>%
            arrange(log2FC)       %>% pull(feature)
feats_kos <- c(up_kos, down_kos)

norm_kos  <- kos_rds$normalized_matrix
avail_kos <- feats_kos[feats_kos %in% rownames(norm_kos)]
mat_kos   <- norm_kos[avail_kos, , drop = FALSE]
rownames(mat_kos) <- sapply(rownames(mat_kos), shorten_ko)

cd_kos           <- kos_rds$col_data
cd_kos$Diagnosis <- ifelse(cd_kos$diagnosis_binary == "IBD", "IBD", "nonIBD")

dir_kos <- setNames(
  ifelse(feats_kos %in% up_kos, "Up in IBD", "Down in IBD"),
  sapply(feats_kos, shorten_ko)
)

draw_heatmap(
  mat           = mat_kos,
  col_data      = cd_kos,
  direction_vec = dir_kos,
  title         = "HMP2 Stool Metaproteomics — KEGG Orthologs (KOS)\nSignificant Features Up & Down (IBD vs nonIBD)",
  filename      = file.path(FIG_DIR, "HMP2_KOS_heatmap.pdf"),
  fontsize_row  = 5.5,
  width         = 16,
  height        = 16
)

# ================================================================
# DONE
# ================================================================
cat("\n================================================================\n")
cat("All heatmaps saved to:", FIG_DIR, "\n")
cat("Files:\n")
for (f in list.files(FIG_DIR, pattern = "\\.pdf$")) {
  cat(" ", f, "\n")
}
