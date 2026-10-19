# ================================================================
# PREPROCESSING: MTBLS8713 — Targeted LC-MS Metabolomics
# IBD (IBDFM*) vs Control (CTRFM*)
# 245 named metabolites x 179 samples
# Input : metabolomics/m_MTBLS8713_LC-MS___metabolite_profiling_v2_maf.tsv
# Output: processed/MTBLS8713_metabolomics/
# ================================================================

library(dplyr)
library(data.table)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "MTBLS8713_metabolomics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("MTBLS8713 — Targeted LC-MS Metabolomics\n")
cat("IBD (IBDFM) vs Control (CTRFM)\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD
# ----------------------------------------------------------------
cat("Step 1: Loading TSV file...\n")
tsv_path <- file.path(BASE,
  "metabolomics/m_MTBLS8713_LC-MS___metabolite_profiling_v2_maf.tsv")

dat <- fread(tsv_path, sep = "\t", header = TRUE,
             data.table = FALSE, check.names = FALSE)

cat("  Rows:", nrow(dat), "| Columns:", ncol(dat), "\n")

# ----------------------------------------------------------------
# 2. IDENTIFY SAMPLE COLUMNS AND METABOLITE NAMES
# ----------------------------------------------------------------
cat("\nStep 2: Identifying sample columns...\n")

# Annotation columns are first 18; sample columns start after
# Column 5 (index 4) = metabolite_identification
# Sample columns: CTRFM* and IBDFM*
ctrl_cols <- grep("^CTRFM", colnames(dat), value = TRUE)
ibd_cols  <- grep("^IBDFM", colnames(dat), value = TRUE)
all_sample_cols <- c(ctrl_cols, ibd_cols)

cat("  Control samples (CTRFM):", length(ctrl_cols), "\n")
cat("  IBD samples (IBDFM):",     length(ibd_cols),  "\n")
cat("  Total samples:",           length(all_sample_cols), "\n")

# Metabolite names from column 5
met_names <- dat[["metabolite_identification"]]
met_names <- make.unique(met_names)
cat("  Metabolites:", length(met_names), "\n")
cat("  First 5:", paste(met_names[1:5], collapse=", "), "\n")

# ----------------------------------------------------------------
# 3. BUILD NUMERIC MATRIX
# ----------------------------------------------------------------
cat("\nStep 3: Building numeric matrix...\n")

mat_df <- as.data.frame(lapply(
  dat[, all_sample_cols, drop = FALSE],
  function(x) suppressWarnings(as.numeric(as.character(x)))
))
rownames(mat_df) <- met_names
mat <- as.matrix(mat_df)

cat("  Matrix:", nrow(mat), "metabolites x", ncol(mat), "samples\n")

# ----------------------------------------------------------------
# 4. BUILD SAMPLE METADATA
# ----------------------------------------------------------------
cat("\nStep 4: Building sample metadata...\n")

col_data <- data.frame(
  sample_id        = all_sample_cols,
  diagnosis        = ifelse(grepl("^IBDFM", all_sample_cols), "IBD", "Control"),
  diagnosis_binary = ifelse(grepl("^IBDFM", all_sample_cols), "IBD", "nonIBD"),
  row.names        = all_sample_cols,
  stringsAsFactors = FALSE
)

cat("  Sample breakdown:\n")
print(table(col_data$diagnosis_binary))

# ----------------------------------------------------------------
# 5. HANDLE MISSING VALUES
#    Data is pre-normalised (range ~0.07 to 3.6)
#    Zeros = missing/below detection, not true zeros
# ----------------------------------------------------------------
cat("\nStep 5: Replacing zeros with NA...\n")
mat[mat == 0 | mat == "0"] <- NA
miss_pct <- round(mean(is.na(mat)) * 100, 1)
cat("  Missing values:", miss_pct, "%\n")

# ----------------------------------------------------------------
# 6. FILTER METABOLITES — remove near-zero variance
#    Also remove any metabolite missing in >50% of samples
# ----------------------------------------------------------------
cat("\nStep 6: Filtering metabolites...\n")
cat("  Before:", nrow(mat), "\n")

# Detection rate filter
det_rate <- rowMeans(!is.na(mat))
mat      <- mat[det_rate >= 0.5, ]
cat("  After detection filter (>50%):", nrow(mat), "\n")

# Near-zero variance filter
row_vars <- apply(mat, 1, var, na.rm = TRUE)
mat      <- mat[!is.na(row_vars) & row_vars > 1e-6, ]
cat("  After variance filter:", nrow(mat), "\n")

# ----------------------------------------------------------------
# 7. LOG2 TRANSFORM
#    Data is already normalised ratios — log2 improves normality
# ----------------------------------------------------------------
cat("\nStep 7: Log2 transforming...\n")
# Values are ratios so we use log2(x + small constant) to handle near-zeros
min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
cat("  Minimum non-zero value:", round(min_nonzero, 4), "\n")
mat_log2 <- log2(mat + min_nonzero * 0.01)
cat("  Range after log2:", round(range(mat_log2, na.rm=TRUE), 2), "\n")

# ----------------------------------------------------------------
# 8. IMPUTE REMAINING NAs
# ----------------------------------------------------------------
cat("\nStep 8: Imputing remaining missing values...\n")
global_mean <- mean(mat_log2, na.rm = TRUE)
global_sd   <- sd(mat_log2,   na.rm = TRUE)
impute_val  <- global_mean - 1.8 * global_sd
cat("  Impute value:", round(impute_val, 3), "\n")
mat_log2[is.na(mat_log2)] <- impute_val

# ----------------------------------------------------------------
# 9. DIFFERENTIAL ANALYSIS — Wilcoxon test
# ----------------------------------------------------------------
cat("\nStep 9: Wilcoxon rank-sum test (IBD vs nonIBD)...\n")

ibd_s  <- col_data$sample_id[col_data$diagnosis_binary == "IBD"]
ctrl_s <- col_data$sample_id[col_data$diagnosis_binary == "nonIBD"]

mean_ibd  <- rowMeans(mat_log2[, ibd_s,  drop = FALSE])
mean_ctrl <- rowMeans(mat_log2[, ctrl_s, drop = FALSE])

pvals <- apply(mat_log2, 1, function(row) {
  ibd_v  <- as.numeric(row[ibd_s])
  ctrl_v <- as.numeric(row[ctrl_s])
  if (sum(!is.na(ibd_v)) < 3 || sum(!is.na(ctrl_v)) < 3) return(NA)
  tryCatch(wilcox.test(ibd_v, ctrl_v, exact = FALSE)$p.value,
           error = function(e) NA)
})

de_results <- data.frame(
  metabolite = rownames(mat_log2),
  mean_IBD   = mean_ibd,
  mean_ctrl  = mean_ctrl,
  log2FC     = mean_ibd - mean_ctrl,
  pval       = pvals,
  stringsAsFactors = FALSE
)
de_results$padj <- p.adjust(de_results$pval, method = "BH")
de_results       <- de_results[order(de_results$padj), ]

sig_features <- de_results[!is.na(de_results$padj) &
                              de_results$padj < 0.05 &
                              abs(de_results$log2FC) > 0.3, ]

cat("  Significant metabolites (padj<0.05, |LFC|>0.3):", nrow(sig_features), "\n")
if (nrow(sig_features) > 0) {
  cat("  Top 10:\n")
  print(head(sig_features[, c("metabolite","log2FC","pval","padj")], 10))
}

# ----------------------------------------------------------------
# 10. ML MATRIX — top 50 features
# ----------------------------------------------------------------
cat("\nStep 10: Building ML matrix (top 50 features)...\n")
top50  <- de_results[!is.na(de_results$padj), ][1:min(50, sum(!is.na(de_results$padj))), "metabolite"]
ml_mat <- t(mat_log2[top50, ])
cat("  ML matrix:", nrow(ml_mat), "samples x", ncol(ml_mat), "metabolites\n")

# ----------------------------------------------------------------
# 11. SAVE
# ----------------------------------------------------------------
cat("\nStep 11: Saving outputs...\n")

saveRDS(list(
  normalized_matrix = mat_log2,
  ml_matrix         = ml_mat,
  col_data          = col_data,
  de_results        = de_results,
  sig_features      = sig_features,
  dataset_info      = list(
    name       = "MTBLS8713",
    data_type  = "targeted_metabolomics",
    n_features = nrow(mat_log2),
    n_samples  = ncol(mat_log2),
    n_IBD      = sum(col_data$diagnosis_binary == "IBD"),
    n_nonIBD   = sum(col_data$diagnosis_binary == "nonIBD"),
    n_sig      = nrow(sig_features)
  )
), file.path(OUT_DIR, "MTBLS8713_processed.rds"))

write.csv(as.data.frame(ml_mat),
          file.path(OUT_DIR, "MTBLS8713_ml_matrix.csv"))
write.csv(de_results,
          file.path(OUT_DIR, "MTBLS8713_DE_results.csv"),
          row.names = FALSE)
write.csv(sig_features,
          file.path(OUT_DIR, "MTBLS8713_sig_features.csv"),
          row.names = FALSE)
write.csv(col_data,
          file.path(OUT_DIR, "MTBLS8713_col_data.csv"))

cat("  Saved: MTBLS8713_processed.rds\n")
cat("  Saved: MTBLS8713_ml_matrix.csv\n")
cat("  Saved: MTBLS8713_DE_results.csv\n")
cat("  Saved: MTBLS8713_sig_features.csv\n")

cat("\n=== DONE ===\n")
cat("Metabolites quantified:", nrow(mat_log2), "\n")
cat("Samples:", ncol(mat_log2),
    "(IBD:", sum(col_data$diagnosis_binary=="IBD"),
    "| nonIBD:", sum(col_data$diagnosis_binary=="nonIBD"), ")\n")
cat("Significant metabolites:", nrow(sig_features), "\n")
