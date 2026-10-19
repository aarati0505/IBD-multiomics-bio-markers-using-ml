# ================================================================
# PREPROCESSING: HMP2 Stool Metaproteomics (ECS + KOS)
# Source  : IBDMDB / HMP2 longitudinal cohort
# Samples : 450 stool samples — CD, UC, nonIBD
# Files   : HMP2_proteomics_ecs.tsv  (910 enzyme functions x 450 samples)
#           HMP2_proteomics_kos.tsv (1823 KEGG orthologs  x 450 samples)
# Outputs : processed/HMP2_metaproteomics/
#             HMP2_metaproteomics_ecs_processed.rds
#             HMP2_metaproteomics_kos_processed.rds
#             HMP2_metaproteomics_ecs_ml_matrix.csv
#             HMP2_metaproteomics_kos_ml_matrix.csv
#             HMP2_metaproteomics_ecs_DE_results.csv
#             HMP2_metaproteomics_kos_DE_results.csv
#             HMP2_metaproteomics_col_data.csv
# ================================================================

library(dplyr)
library(data.table)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "HMP2_metaproteomics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("HMP2 Stool Metaproteomics Preprocessing\n")
cat("ECS (enzyme functions) + KOS (KEGG orthologs)\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# HELPER: full preprocessing pipeline for one PCL-style TSV
# ----------------------------------------------------------------
preprocess_hmp2_mpx <- function(tsv_path, id_col, label, meta_df,
                                 ibd_col = "diagnosis_binary",
                                 out_dir  = OUT_DIR) {

  cat("\n================================================================\n")
  cat("Processing:", label, "\n")
  cat("================================================================\n")

  # -- 1. LOAD ---------------------------------------------------
  cat("Step 1: Loading", basename(tsv_path), "...\n")
  dat <- fread(tsv_path, sep = "\t", header = TRUE, data.table = FALSE)
  cat("  Dimensions:", nrow(dat), "features x", ncol(dat), "columns\n")
  cat("  ID column:", colnames(dat)[1], "\n")

  # Feature IDs are in column 1; rest are sample columns
  feature_ids  <- dat[[1]]
  sample_cols  <- colnames(dat)[-1]

  # -- 2. REMOVE UNGROUPED / UNMAPPED ----------------------------
  cat("Step 2: Removing UNGROUPED / UNMAPPED rows...\n")
  keep_rows <- !feature_ids %in% c("UNGROUPED", "UNMAPPED", "")
  dat         <- dat[keep_rows, ]
  feature_ids <- feature_ids[keep_rows]
  cat("  After removal:", nrow(dat), "features\n")

  # -- 3. BUILD NUMERIC MATRIX -----------------------------------
  cat("Step 3: Building numeric matrix...\n")
  mat <- as.matrix(dat[, sample_cols, drop = FALSE])
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- make.unique(feature_ids)
  cat("  Matrix:", nrow(mat), "x", ncol(mat), "\n")

  # -- 4. ATTACH METADATA ----------------------------------------
  cat("Step 4: Attaching sample metadata...\n")

  # Filter metadata to proteomics rows only, keep one row per External ID
  prot_meta <- meta_df %>%
    filter(data_type == "proteomics") %>%
    distinct(`External ID`, .keep_all = TRUE)

  # Match sample columns to metadata
  matched <- prot_meta[match(sample_cols, prot_meta$`External ID`), ]
  has_meta <- !is.na(matched$`External ID`)
  cat("  Samples with metadata:", sum(has_meta), "/", length(sample_cols), "\n")

  col_data <- data.frame(
    sample_id        = sample_cols,
    participant_id   = matched$`Participant ID`,
    diagnosis        = matched$diagnosis,
    diagnosis_binary = matched[[ibd_col]],
    week_num         = matched$week_num,
    sex              = matched$sex,
    site_name        = matched$site_name,
    row.names        = sample_cols,
    stringsAsFactors = FALSE
  )

  # Drop samples without metadata
  valid_samples <- sample_cols[has_meta]
  col_data      <- col_data[has_meta, ]
  mat           <- mat[, valid_samples, drop = FALSE]
  cat("  Retained:", ncol(mat), "samples with complete metadata\n")

  cat("  Diagnosis breakdown:\n")
  print(table(col_data$diagnosis_binary))

  # -- 5. REPLACE ZEROS WITH NA ----------------------------------
  cat("Step 5: Replacing zeros with NA...\n")
  mat[mat == 0] <- NA
  miss_pct <- round(mean(is.na(mat)) * 100, 1)
  cat("  Missing (zero-replaced):", miss_pct, "%\n")

  # -- 6. FILTER: keep features detected in >30% of samples ------
  #    (lower threshold than tissue data — stool is sparser)
  cat("Step 6: Filtering by detection rate (>30% in at least one group)...\n")
  cat("  Before:", nrow(mat), "features\n")

  ibd_s  <- col_data$sample_id[col_data$diagnosis_binary == "IBD"]
  ctrl_s <- col_data$sample_id[col_data$diagnosis_binary == "nonIBD"]

  det_ibd  <- rowMeans(!is.na(mat[, ibd_s,  drop = FALSE]))
  det_ctrl <- rowMeans(!is.na(mat[, ctrl_s, drop = FALSE]))
  keep     <- (det_ibd >= 0.30) | (det_ctrl >= 0.30)
  mat      <- mat[keep, ]
  cat("  After:", nrow(mat), "features\n")

  # -- 7. LOG2 TRANSFORM (counts + 1 pseudocount) ----------------
  cat("Step 7: Log2(count + 1) transforming...\n")
  mat_log2 <- log2(mat + 1)
  mat_log2[is.na(mat_log2)] <- 0   # undetected features -> 0 after log
  cat("  Range:", round(range(mat_log2, na.rm = TRUE), 2), "\n")

  # -- 8. RELATIVE ABUNDANCE NORMALISATION (TSS per sample) ------
  #    Divide each sample by its column sum, then re-log
  #    This corrects for library depth differences between samples
  cat("Step 8: Total-sum scaling (TSS) per sample...\n")
  raw_mat   <- mat      # keep raw for imputation reference
  col_sums  <- colSums(mat, na.rm = TRUE)
  mat_tss   <- sweep(mat, 2, col_sums, "/") * 1e6  # CPM-style
  mat_log2  <- log2(mat_tss + 1)
  mat_log2[is.na(mat_log2)] <- 0
  cat("  Normalisation complete\n")

  # -- 9. DIFFERENTIAL ANALYSIS: Wilcoxon per feature -------------
  #    Using Wilcoxon (non-parametric) because count data is skewed
  cat("Step 9: Wilcoxon rank-sum test (IBD vs nonIBD)...\n")

  mean_ibd  <- rowMeans(mat_log2[, ibd_s,  drop = FALSE])
  mean_ctrl <- rowMeans(mat_log2[, ctrl_s, drop = FALSE])

  pvals <- apply(mat_log2, 1, function(row) {
    ibd_v  <- as.numeric(row[ibd_s])
    ctrl_v <- as.numeric(row[ctrl_s])
    if (sum(ibd_v != 0) < 3 || sum(ctrl_v != 0) < 3) return(NA)
    tryCatch(wilcox.test(ibd_v, ctrl_v, exact = FALSE)$p.value,
             error = function(e) NA)
  })

  de_results <- data.frame(
    feature   = rownames(mat_log2),
    mean_IBD  = mean_ibd,
    mean_ctrl = mean_ctrl,
    log2FC    = mean_ibd - mean_ctrl,
    pval      = pvals,
    stringsAsFactors = FALSE
  )
  de_results$padj <- p.adjust(de_results$pval, method = "BH")
  de_results       <- de_results %>% arrange(padj)

  sig_features <- de_results %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FC) > 0.5)

  cat("  Significant features (padj<0.05, |LFC|>0.5):", nrow(sig_features), "\n")
  if (nrow(sig_features) > 0) {
    cat("  Top 10:\n")
    print(head(sig_features[, c("feature", "log2FC", "pval", "padj")], 10))
  }

  # -- 10. ML MATRIX: top 100 by padj ---------------------------
  cat("Step 10: Building ML matrix (top 100 features)...\n")
  top100  <- de_results %>% filter(!is.na(padj)) %>% slice_head(n = 100) %>% pull(feature)
  ml_mat  <- t(mat_log2[top100, ])   # samples x features
  cat("  ML matrix:", nrow(ml_mat), "samples x", ncol(ml_mat), "features\n")

  # -- 11. SAVE --------------------------------------------------
  cat("Step 11: Saving outputs...\n")

  tag <- tolower(gsub("[^A-Za-z0-9]", "_", label))

  saveRDS(list(
    normalized_matrix = mat_log2,
    ml_matrix         = ml_mat,
    col_data          = col_data,
    de_results        = de_results,
    sig_features      = sig_features,
    dataset_info      = list(
      name          = paste0("HMP2_", label),
      data_type     = "metaproteomics",
      feature_type  = label,
      n_features    = nrow(mat_log2),
      n_samples     = ncol(mat_log2),
      n_IBD         = sum(col_data$diagnosis_binary == "IBD"),
      n_nonIBD      = sum(col_data$diagnosis_binary == "nonIBD"),
      n_sig         = nrow(sig_features)
    )
  ), file.path(out_dir, paste0("HMP2_metaproteomics_", tag, "_processed.rds")))

  write.csv(as.data.frame(ml_mat),
            file.path(out_dir, paste0("HMP2_metaproteomics_", tag, "_ml_matrix.csv")))
  write.csv(de_results,
            file.path(out_dir, paste0("HMP2_metaproteomics_", tag, "_DE_results.csv")),
            row.names = FALSE)
  write.csv(sig_features,
            file.path(out_dir, paste0("HMP2_metaproteomics_", tag, "_sig_features.csv")),
            row.names = FALSE)

  cat("  Saved:", paste0("HMP2_metaproteomics_", tag, "_processed.rds"), "\n")
  cat("  Saved:", paste0("HMP2_metaproteomics_", tag, "_ml_matrix.csv"), "\n")
  cat("  Saved:", paste0("HMP2_metaproteomics_", tag, "_DE_results.csv"), "\n")
  cat("  Saved:", paste0("HMP2_metaproteomics_", tag, "_sig_features.csv"), "\n")

  invisible(list(
    col_data     = col_data,
    de_results   = de_results,
    sig_features = sig_features,
    ml_mat       = ml_mat
  ))
}

# ----------------------------------------------------------------
# 1. LOAD METADATA
# ----------------------------------------------------------------
cat("Loading HMP2 metadata...\n")
meta <- read.csv(file.path(BASE, "data/hmp2_metadata_clean.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
cat("  Metadata rows:", nrow(meta), "\n")

# ----------------------------------------------------------------
# 2. PROCESS ECS (Enzyme Commission rollups)
# ----------------------------------------------------------------
ecs_result <- preprocess_hmp2_mpx(
  tsv_path = file.path(BASE, "proteomics/HMP2_proteomics_ecs.tsv"),
  id_col   = "Gene",
  label    = "ECS",
  meta_df  = meta
)

# ----------------------------------------------------------------
# 3. PROCESS KOS (KEGG Orthology rollups)
# ----------------------------------------------------------------
kos_result <- preprocess_hmp2_mpx(
  tsv_path = file.path(BASE, "proteomics/HMP2_proteomics_kos.tsv"),
  id_col   = "KO",
  label    = "KOS",
  meta_df  = meta
)

# ----------------------------------------------------------------
# 4. SAVE SHARED COLUMN METADATA (same 450 samples for both)
# ----------------------------------------------------------------
cat("\nSaving shared sample metadata...\n")
write.csv(ecs_result$col_data,
          file.path(OUT_DIR, "HMP2_metaproteomics_col_data.csv"))

# ----------------------------------------------------------------
# 5. FINAL SUMMARY
# ----------------------------------------------------------------
cat("\n================================================================\n")
cat("=== DONE — HMP2 Metaproteomics Summary ===\n")
cat("================================================================\n")
cat("ECS:\n")
cat("  Features after QC   :", nrow(ecs_result$de_results), "\n")
cat("  Significant features:", nrow(ecs_result$sig_features), "\n")
cat("  Samples             :", nrow(ecs_result$col_data), "\n")
cat("  IBD / nonIBD        :",
    sum(ecs_result$col_data$diagnosis_binary == "IBD"), "/",
    sum(ecs_result$col_data$diagnosis_binary == "nonIBD"), "\n")
cat("\nKOS:\n")
cat("  Features after QC   :", nrow(kos_result$de_results), "\n")
cat("  Significant features:", nrow(kos_result$sig_features), "\n")
cat("  Samples             :", nrow(kos_result$col_data), "\n")
cat("  IBD / nonIBD        :",
    sum(kos_result$col_data$diagnosis_binary == "IBD"), "/",
    sum(kos_result$col_data$diagnosis_binary == "nonIBD"), "\n")
