# ================================================================
# SUMMARY HEATMAPS — Group Averages (IBD vs nonIBD)
# Instead of one column per patient, shows:
#   Column 1 = mean expression across all IBD samples
#   Column 2 = mean expression across all nonIBD samples
#
# This makes the up/down pattern immediately visible
# and is much cleaner for presentations.
#
# Datasets:
#   1. GSE137344  — Ileum transcriptomics
#   2. GSE186507  — Blood transcriptomics
#   3. GSE193677  — Biopsy transcriptomics
#   4. PXD001608  — Colon tissue proteomics
#   5. HMP2 ECS   — Stool metaproteomics (enzyme functions)
#   6. HMP2 KOS   — Stool metaproteomics (KEGG orthologs)
#
# Output: figures/heatmaps/summary/
# ================================================================

library(pheatmap)
library(RColorBrewer)
library(dplyr)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "figures", "heatmaps", "summary")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------
# COLOUR PALETTE
# ----------------------------------------------------------------
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

anno_colors <- list(
  Group = c(IBD = "#D73027", nonIBD = "#4575B4")
)

# ----------------------------------------------------------------
# HELPER: build a 2-column summary matrix (IBD mean vs nonIBD mean)
# then draw and save heatmap
# ----------------------------------------------------------------
summary_heatmap <- function(mat,            # features x samples matrix
                             diagnosis_vec,  # named vector: sample -> "IBD"/"nonIBD"
                             features,       # which features (rows) to show
                             title,
                             filename,
                             up_label   = "Higher in IBD →",
                             down_label = "← Lower in IBD",
                             fontsize_row = 9,
                             width  = 5,
                             height = 10) {

  # Keep only features present in matrix
  features <- features[features %in% rownames(mat)]

  ibd_cols  <- names(diagnosis_vec)[diagnosis_vec == "IBD"]
  ctrl_cols <- names(diagnosis_vec)[diagnosis_vec == "nonIBD"]

  ibd_cols  <- ibd_cols[ibd_cols   %in% colnames(mat)]
  ctrl_cols <- ctrl_cols[ctrl_cols %in% colnames(mat)]

  # Compute group means
  mean_ibd  <- rowMeans(mat[features, ibd_cols,  drop = FALSE], na.rm = TRUE)
  mean_ctrl <- rowMeans(mat[features, ctrl_cols, drop = FALSE], na.rm = TRUE)

  # Build 2-column summary matrix
  summary_mat <- cbind(IBD = mean_ibd, nonIBD = mean_ctrl)

  # Z-score rows so colour scale is relative (same as per-patient heatmaps)
  summary_mat <- t(scale(t(summary_mat)))
  summary_mat[summary_mat >  2] <-  2
  summary_mat[summary_mat < -2] <- -2

  # Annotation for the 2 columns
  anno_col <- data.frame(
    Group = c("IBD", "nonIBD"),
    row.names = c("IBD", "nonIBD")
  )

  # Add up/down arrows to row names for clarity
  log2fc    <- mean_ibd - mean_ctrl
  row_labels <- paste0(
    ifelse(log2fc > 0, "▲ ", "▼ "),
    features
  )
  rownames(summary_mat) <- row_labels

  pdf(filename, width = width, height = height)
  pheatmap(
    summary_mat,
    color             = heatmap_colors,
    annotation_col    = anno_col,
    annotation_colors = anno_colors,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = fontsize_row,
    fontsize_col      = 11,
    fontsize          = 9,
    main              = title,
    border_color      = "grey90",
    treeheight_row    = 30,
    treeheight_col    = 0,
    cellwidth         = 40,
    angle_col         = 0
  )
  dev.off()
  cat("  Saved:", basename(filename), "\n")
}

# ================================================================
# 1. GSE137344 — Ileum Transcriptomics
# ================================================================
cat("\n--- GSE137344: Ileum Transcriptomics ---\n")

rds137 <- readRDS(file.path(BASE, "processed/GSE137344/GSE137344_normalized.rds"))
sig137 <- read.csv(file.path(BASE, "processed/GSE137344/GSE137344_significant_DEGs.csv"),
                   stringsAsFactors = FALSE)

# Top 25 up + 25 down by fold change
up137   <- sig137 %>% filter(log2FoldChange > 0) %>%
           arrange(desc(log2FoldChange)) %>% slice_head(n = 25) %>% pull(gene)
down137 <- sig137 %>% filter(log2FoldChange < 0) %>%
           arrange(log2FoldChange)       %>% slice_head(n = 25) %>% pull(gene)
feats137 <- c(up137, down137)

diag137 <- setNames(rds137$col_data$diagnosis_binary,
                    rds137$col_data$sample_id)

cat("  Showing:", length(up137), "up +", length(down137), "down\n")

summary_heatmap(
  mat           = rds137$normalized_matrix,
  diagnosis_vec = diag137,
  features      = feats137,
  title         = "GSE137344 — Ileum Transcriptomics\nGroup Average: IBD vs nonIBD\n(▲ = higher in IBD, ▼ = lower in IBD)",
  filename      = file.path(OUT_DIR, "GSE137344_summary_heatmap.pdf"),
  fontsize_row  = 8,
  height        = 13
)

# ================================================================
# 2. GSE186507 — Blood Transcriptomics
# ================================================================
cat("\n--- GSE186507: Blood Transcriptomics ---\n")

rds186 <- readRDS(file.path(BASE, "processed/GSE186507/GSE186507_normalized.rds"))
deg186 <- read.csv(file.path(BASE, "processed/GSE186507/GSE186507_DEG_results.csv"),
                   stringsAsFactors = FALSE)

sig186 <- deg186 %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 0.5) %>%
  arrange(adj.P.Val)

up186   <- sig186 %>% filter(logFC > 0) %>%
           arrange(desc(logFC)) %>% slice_head(n = 25) %>% pull(gene)
down186 <- sig186 %>% filter(logFC < 0) %>%
           arrange(logFC)       %>% slice_head(n = 25) %>% pull(gene)
feats186 <- c(up186, down186)

diag186 <- setNames(rds186$col_data$diagnosis_binary,
                    rds186$col_data$sample_id)

cat("  Showing:", length(up186), "up +", length(down186), "down\n")

summary_heatmap(
  mat           = rds186$normalized_matrix,
  diagnosis_vec = diag186,
  features      = feats186,
  title         = "GSE186507 — Blood Transcriptomics\nGroup Average: IBD vs nonIBD\n(▲ = higher in IBD, ▼ = lower in IBD)",
  filename      = file.path(OUT_DIR, "GSE186507_summary_heatmap.pdf"),
  fontsize_row  = 8,
  height        = 13
)

# ================================================================
# 3. GSE193677 — Biopsy Transcriptomics
# ================================================================
cat("\n--- GSE193677: Biopsy Transcriptomics ---\n")

rds193 <- readRDS(file.path(BASE, "processed/GSE193677/GSE193677_normalized.rds"))
deg193 <- read.csv(file.path(BASE, "processed/GSE193677/GSE193677_DEG_results.csv"),
                   stringsAsFactors = FALSE)

# Top 25 up + 25 down by raw p-value (BH correction yields only 3 genes)
top_by_p <- deg193 %>% filter(!is.na(P.Value)) %>% arrange(P.Value) %>% slice_head(n = 200)
up193    <- top_by_p %>% filter(logFC > 0) %>% slice_head(n = 25) %>% pull(gene)
down193  <- top_by_p %>% filter(logFC < 0) %>% slice_head(n = 25) %>% pull(gene)
feats193 <- c(up193, down193)

diag193 <- setNames(rds193$col_data$diagnosis_binary,
                    rds193$col_data$sample_id)

cat("  Showing:", length(up193), "up +", length(down193), "down\n")

summary_heatmap(
  mat           = rds193$normalized_matrix,
  diagnosis_vec = diag193,
  features      = feats193,
  title         = "GSE193677 — Biopsy Transcriptomics\nGroup Average: IBD vs nonIBD\n(▲ = higher in IBD, ▼ = lower in IBD)",
  filename      = file.path(OUT_DIR, "GSE193677_summary_heatmap.pdf"),
  fontsize_row  = 8,
  height        = 13
)

# ================================================================
# 4. PXD001608 — Colon Tissue Proteomics
# ================================================================
cat("\n--- PXD001608: Colon Tissue Proteomics ---\n")

pxd_rds  <- readRDS(file.path(BASE,
             "processed/proteomics_MaxQuant/proteomics_MaxQuant_processed.rds"))
sig_prot <- read.csv(file.path(BASE,
             "processed/proteomics_MaxQuant/proteomics_MaxQuant_sig_proteins.csv"),
             stringsAsFactors = FALSE)

up_prot   <- sig_prot %>% filter(log2FC > 0) %>%
             arrange(desc(log2FC)) %>% slice_head(n = 25) %>% pull(protein)
down_prot <- sig_prot %>% filter(log2FC < 0) %>%
             arrange(log2FC)       %>% slice_head(n = 25) %>% pull(protein)
feats_pxd <- c(up_prot, down_prot)

diag_pxd <- setNames(pxd_rds$col_data$diagnosis_binary,
                     pxd_rds$col_data$sample_id)

cat("  Showing:", length(up_prot), "up +", length(down_prot), "down\n")

summary_heatmap(
  mat           = pxd_rds$normalized_matrix,
  diagnosis_vec = diag_pxd,
  features      = feats_pxd,
  title         = "PXD001608 — Colon Tissue Proteomics\nGroup Average: IBD (UC) vs nonIBD\n(▲ = higher in IBD, ▼ = lower in IBD)",
  filename      = file.path(OUT_DIR, "PXD001608_summary_heatmap.pdf"),
  fontsize_row  = 8,
  height        = 13
)

# ================================================================
# 5. HMP2 ECS — Stool Metaproteomics (Enzyme Functions)
# ================================================================
cat("\n--- HMP2 ECS: Stool Metaproteomics ---\n")

ecs_rds <- readRDS(file.path(BASE,
            "processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_processed.rds"))
sig_ecs <- read.csv(file.path(BASE,
            "processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_sig_features.csv"),
            stringsAsFactors = FALSE)

# Shorten EC names for display
shorten_ecs <- function(x) {
  ec <- regmatches(x, regexpr("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+", x))
  nm <- sub("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+:\\s*", "", x)
  nm <- substr(nm, 1, 40)
  if (length(ec) > 0 && nchar(ec) > 0) paste0(ec, ": ", nm) else nm
}

up_ecs   <- sig_ecs %>% filter(log2FC > 0) %>%
            arrange(desc(log2FC)) %>% pull(feature)
down_ecs <- sig_ecs %>% filter(log2FC < 0) %>%
            arrange(log2FC)       %>%
            slice_head(n = 30)   %>% pull(feature)
feats_ecs <- c(up_ecs, down_ecs)

diag_ecs <- setNames(ecs_rds$col_data$diagnosis_binary,
                     ecs_rds$col_data$sample_id)

# Rename rows for display
norm_ecs_renamed <- ecs_rds$normalized_matrix
rownames(norm_ecs_renamed) <- sapply(rownames(norm_ecs_renamed), shorten_ecs)
feats_ecs_short  <- sapply(feats_ecs, shorten_ecs)

cat("  Showing:", length(up_ecs), "up +", length(down_ecs), "down\n")

summary_heatmap(
  mat           = norm_ecs_renamed,
  diagnosis_vec = diag_ecs,
  features      = feats_ecs_short,
  title         = "HMP2 Stool Metaproteomics — Enzyme Functions (ECS)\nGroup Average: IBD vs nonIBD\n(▲ = higher in IBD, ▼ = lower in IBD)",
  filename      = file.path(OUT_DIR, "HMP2_ECS_summary_heatmap.pdf"),
  fontsize_row  = 7,
  height        = 14
)

# ================================================================
# 6. HMP2 KOS — Stool Metaproteomics (KEGG Orthologs)
# ================================================================
cat("\n--- HMP2 KOS: Stool Metaproteomics (KEGG Orthologs) ---\n")

kos_rds <- readRDS(file.path(BASE,
            "processed/HMP2_metaproteomics/HMP2_metaproteomics_kos_processed.rds"))
sig_kos <- read.csv(file.path(BASE,
            "processed/HMP2_metaproteomics/HMP2_metaproteomics_kos_sig_features.csv"),
            stringsAsFactors = FALSE)

shorten_ko <- function(x) {
  ko <- regmatches(x, regexpr("^K[0-9]+", x))
  nm <- sub("^K[0-9]+:\\s*", "", x)
  nm <- sub("\\s*\\[EC:[^\\]]+\\]", "", nm)
  nm <- substr(nm, 1, 38)
  if (length(ko) > 0 && nchar(ko) > 0) paste0(ko, ": ", nm) else nm
}

up_kos   <- sig_kos %>% filter(log2FC > 0) %>%
            arrange(desc(log2FC)) %>% pull(feature)
down_kos <- sig_kos %>% filter(log2FC < 0) %>%
            arrange(log2FC)       %>%
            slice_head(n = 30)   %>% pull(feature)
feats_kos <- c(up_kos, down_kos)

diag_kos <- setNames(kos_rds$col_data$diagnosis_binary,
                     kos_rds$col_data$sample_id)

norm_kos_renamed <- kos_rds$normalized_matrix
rownames(norm_kos_renamed) <- sapply(rownames(norm_kos_renamed), shorten_ko)
feats_kos_short  <- sapply(feats_kos, shorten_ko)

cat("  Showing:", length(up_kos), "up +", length(down_kos), "down\n")

summary_heatmap(
  mat           = norm_kos_renamed,
  diagnosis_vec = diag_kos,
  features      = feats_kos_short,
  title         = "HMP2 Stool Metaproteomics — KEGG Orthologs (KOS)\nGroup Average: IBD vs nonIBD\n(▲ = higher in IBD, ▼ = lower in IBD)",
  filename      = file.path(OUT_DIR, "HMP2_KOS_summary_heatmap.pdf"),
  fontsize_row  = 7,
  height        = 14
)

# ================================================================
# DONE
# ================================================================
cat("\n================================================================\n")
cat("All summary heatmaps saved to:", OUT_DIR, "\n")
for (f in list.files(OUT_DIR, pattern = "\\.pdf$")) cat(" ", f, "\n")
