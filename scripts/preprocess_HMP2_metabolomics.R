# ================================================================
# HMP2 Metabolomics Preprocessing — R pipeline
# Input : processed/HMP2_metabolomics/HMP2_metabolomics_matrix.csv
#         processed/HMP2_metabolomics/HMP2_metabolomics_col_data.csv
# Output: processed/HMP2_metabolomics/
#           HMP2_metabolomics_processed.rds
#           HMP2_metabolomics_ml_matrix.csv
#           HMP2_metabolomics_DE_results.csv
#           HMP2_metabolomics_sig_features.csv
# ================================================================

library(dplyr)
library(data.table)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "HMP2_metabolomics")

cat("================================================================\n")
cat("HMP2 Metabolomics Preprocessing\n")
cat("75,447 features x 546 samples (IBD: 411 | nonIBD: 135)\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD CSV (already filtered to >10% detection by Python)
# ----------------------------------------------------------------
cat("Step 1: Loading matrix from CSV...\n")
t0  <- proc.time()
mat_df <- fread(file.path(OUT_DIR, "HMP2_metabolomics_matrix.csv"),
                header = TRUE, data.table = FALSE)
# First column is feature IDs
feat_ids    <- mat_df[[1]]
mat_df      <- mat_df[, -1, drop = FALSE]
mat         <- as.matrix(mat_df)
rownames(mat) <- feat_ids
cat("  Matrix:", nrow(mat), "features x", ncol(mat), "samples\n")
cat("  Loaded in", round((proc.time()-t0)["elapsed"], 1), "sec\n")
rm(mat_df); gc()

# ----------------------------------------------------------------
# 2. LOAD SAMPLE METADATA
# ----------------------------------------------------------------
cat("\nStep 2: Loading sample metadata...\n")
col_data <- read.csv(file.path(OUT_DIR, "HMP2_metabolomics_col_data.csv"),
                     stringsAsFactors = FALSE, row.names = 1)
col_data$sample_id <- rownames(col_data)
cat("  Samples:", nrow(col_data), "\n")
cat("  Diagnosis breakdown:\n")
print(table(col_data$diagnosis_binary))

ibd_s  <- rownames(col_data)[col_data$diagnosis_binary == "IBD"]
ctrl_s <- rownames(col_data)[col_data$diagnosis_binary == "nonIBD"]

# Align matrix columns to col_data
mat <- mat[, rownames(col_data), drop = FALSE]

# ----------------------------------------------------------------
# 3. REPLACE ZEROS WITH NA
# ----------------------------------------------------------------
cat("\nStep 3: Replacing zeros with NA...\n")
mat[mat == 0] <- NA
miss_pct <- round(mean(is.na(mat)) * 100, 1)
cat("  Missing:", miss_pct, "%\n")

# ----------------------------------------------------------------
# 4. LOG2 TRANSFORM + TOTAL SUM SCALING (CPM)
# ----------------------------------------------------------------
cat("\nStep 4: TSS normalisation + log2 transform...\n")
mat_nz    <- mat;  mat_nz[is.na(mat_nz)] <- 0
col_sums  <- colSums(mat_nz)
mat_tss   <- sweep(mat_nz, 2, col_sums, "/") * 1e6   # CPM
mat_log2  <- log2(mat_tss + 1)
mat_log2[mat_log2 == 0 & is.na(mat)] <- NA
rm(mat, mat_nz, mat_tss); gc()
cat("  Range (non-NA):", round(range(mat_log2, na.rm = TRUE), 2), "\n")

# ----------------------------------------------------------------
# 5. IMPUTE NAs (global mean - 1.8 SD)
# ----------------------------------------------------------------
cat("\nStep 5: Imputing missing values...\n")
global_mean <- mean(mat_log2, na.rm = TRUE)
global_sd   <- sd(mat_log2,   na.rm = TRUE)
impute_val  <- global_mean - 1.8 * global_sd
cat("  Impute value:", round(impute_val, 3), "\n")
mat_log2[is.na(mat_log2)] <- impute_val

# ----------------------------------------------------------------
# 6. DIFFERENTIAL ANALYSIS — Wilcoxon per feature
# ----------------------------------------------------------------
cat("\nStep 6: Wilcoxon test (IBD vs nonIBD)...\n")
cat("  Testing", nrow(mat_log2), "features — this will take a few minutes...\n")
t1 <- proc.time()

mean_ibd  <- rowMeans(mat_log2[, ibd_s,  drop = FALSE])
mean_ctrl <- rowMeans(mat_log2[, ctrl_s, drop = FALSE])

pvals <- apply(mat_log2, 1, function(r) {
  iv <- r[ibd_s];  cv <- r[ctrl_s]
  # Skip features where almost all values are imputed
  if (sum(iv > impute_val + 0.01) < 5 ||
      sum(cv > impute_val + 0.01) < 5) return(NA)
  tryCatch(wilcox.test(iv, cv, exact = FALSE)$p.value,
           error = function(e) NA)
})

cat("  Done in", round((proc.time()-t1)["elapsed"]/60, 1), "min\n")

de <- data.frame(
  feature   = rownames(mat_log2),
  mean_IBD  = mean_ibd,
  mean_ctrl = mean_ctrl,
  log2FC    = mean_ibd - mean_ctrl,
  pval      = pvals,
  stringsAsFactors = FALSE
)
de$padj <- p.adjust(de$pval, method = "BH")
de       <- de[order(de$padj), ]

sig <- de[!is.na(de$padj) & de$padj < 0.05 & abs(de$log2FC) > 0.5, ]
cat("  Significant (padj<0.05, |LFC|>0.5):", nrow(sig), "\n")
cat("  Top 10:\n")
print(head(sig[, c("feature", "log2FC", "pval", "padj")], 10))

# ----------------------------------------------------------------
# 7. ML MATRIX — top 100 by padj
# ----------------------------------------------------------------
cat("\nStep 7: Building ML matrix (top 100 features)...\n")
top100 <- head(de[!is.na(de$padj), "feature"], 100)
ml_mat <- t(mat_log2[top100, ])
cat("  ML matrix:", nrow(ml_mat), "x", ncol(ml_mat), "\n")

# ----------------------------------------------------------------
# 8. SAVE
# ----------------------------------------------------------------
cat("\nStep 8: Saving outputs...\n")

saveRDS(list(
  normalized_matrix = mat_log2,
  ml_matrix         = ml_mat,
  col_data          = col_data,
  de_results        = de,
  sig_features      = sig,
  dataset_info      = list(
    name     = "HMP2_metabolomics",
    n_feat   = nrow(mat_log2),
    n_samp   = ncol(mat_log2),
    n_IBD    = length(ibd_s),
    n_nonIBD = length(ctrl_s),
    n_sig    = nrow(sig)
  )
), file.path(OUT_DIR, "HMP2_metabolomics_processed.rds"))

write.csv(as.data.frame(ml_mat),
          file.path(OUT_DIR, "HMP2_metabolomics_ml_matrix.csv"))
write.csv(de,  file.path(OUT_DIR, "HMP2_metabolomics_DE_results.csv"),
          row.names = FALSE)
write.csv(sig, file.path(OUT_DIR, "HMP2_metabolomics_sig_features.csv"),
          row.names = FALSE)

cat("  Saved: HMP2_metabolomics_processed.rds\n")
cat("  Saved: HMP2_metabolomics_ml_matrix.csv\n")
cat("  Saved: HMP2_metabolomics_DE_results.csv\n")
cat("  Saved: HMP2_metabolomics_sig_features.csv\n")

cat("\n=== DONE ===\n")
cat("Features:", nrow(mat_log2), "| Samples:", ncol(mat_log2),
    "| Significant:", nrow(sig), "\n")
