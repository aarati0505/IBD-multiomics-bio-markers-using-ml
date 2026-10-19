# ============================================================
# PREPROCESSING: Proteomics — MaxQuant proteinGroups.txt
# UC patients vs healthy controls, triplicates
# Input : proteinGroups.txt
# Output: data/proteomics_processed.rds
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)

dest <- "C:/Users/Aarati/Desktop/IBD datasets"
dir.create(file.path(dest, "data"), showWarnings = FALSE)

cat("=== Proteomics Preprocessing ===\n\n")

# ----------------------------------------------------------
# 1. LOAD proteinGroups.txt
# ----------------------------------------------------------
cat("1. Loading proteinGroups.txt...\n")
pg <- fread(file.path(dest, "proteinGroups.txt"), sep = "\t",
            header = TRUE, data.table = FALSE)
cat("   Dimensions:", nrow(pg), "protein groups x", ncol(pg), "columns\n")

# ----------------------------------------------------------
# 2. IDENTIFY LFQ INTENSITY COLUMNS
# ----------------------------------------------------------
# MaxQuant names them "LFQ intensity <sample>"
# Your samples are: UC_1_1, UC_1_2... and Ctrl_1_1, Ctrl_1_2...
lfq_cols <- grep("^LFQ intensity", colnames(pg), value = TRUE)

if (length(lfq_cols) == 0) {
  # If no LFQ columns, use raw Intensity columns
  lfq_cols <- grep("^Intensity", colnames(pg), value = TRUE)
  lfq_cols <- lfq_cols[lfq_cols != "Intensity"]  # remove aggregate column
  cat("   Using Intensity columns (no LFQ found)\n")
} else {
  cat("   Found", length(lfq_cols), "LFQ intensity columns\n")
}

# ----------------------------------------------------------
# 3. EXTRACT GENE NAMES AND FILTER CONTAMINANTS
# ----------------------------------------------------------
cat("\n2. Filtering contaminants and reverse hits...\n")

# MaxQuant flags contaminants (+) and reverse decoys (+)
if ("Potential contaminant" %in% colnames(pg)) {
  pg <- pg %>% filter(`Potential contaminant` != "+")
}
if ("Reverse" %in% colnames(pg)) {
  pg <- pg %>% filter(Reverse != "+")
}
if ("Only identified by site" %in% colnames(pg)) {
  pg <- pg %>% filter(`Only identified by site` != "+")
}

cat("   After contaminant removal:", nrow(pg), "protein groups\n")

# Get gene names
gene_col <- if ("Gene names" %in% colnames(pg)) "Gene names" else "Protein IDs"
pg$gene_name <- pg[[gene_col]]

# ----------------------------------------------------------
# 4. EXTRACT INTENSITY MATRIX
# ----------------------------------------------------------
cat("\n3. Extracting intensity matrix...\n")

mat <- pg[, lfq_cols, drop = FALSE]
rownames(mat) <- make.unique(pg$gene_name)

# Replace 0 with NA (0 in MaxQuant = not detected, not truly zero)
mat[mat == 0] <- NA

cat("   Matrix:", nrow(mat), "proteins x", ncol(mat), "samples\n")

# ----------------------------------------------------------
# 5. FILTER PROTEINS BY DETECTION RATE
# ----------------------------------------------------------
cat("\n4. Filtering by detection rate (>50% samples)...\n")

detection_rate <- rowMeans(!is.na(mat))
mat_filtered <- mat[detection_rate >= 0.5, ]
cat("   After filter:", nrow(mat_filtered), "proteins\n")

# ----------------------------------------------------------
# 6. LOG2 TRANSFORM
# ----------------------------------------------------------
cat("\n5. Log2 transforming...\n")
mat_log2 <- log2(mat_filtered + 1)

# ----------------------------------------------------------
# 7. IMPUTE MISSING VALUES
# ----------------------------------------------------------
# Use minimum value imputation per protein (conservative, standard in proteomics)
cat("6. Imputing missing values (min-value imputation)...\n")

mat_imputed <- mat_log2
for (i in seq_len(nrow(mat_imputed))) {
  row_vals <- mat_imputed[i, ]
  min_val  <- min(row_vals, na.rm = TRUE)
  mat_imputed[i, is.na(row_vals)] <- min_val - 1  # below detection floor
}

# ----------------------------------------------------------
# 8. BUILD SAMPLE LABELS FROM COLUMN NAMES
# ----------------------------------------------------------
cat("\n7. Building sample labels...\n")

# Column names like "LFQ intensity UC_8_2" or "LFQ intensity Ctrl_1_1"
sample_names <- gsub("LFQ intensity |Intensity ", "", lfq_cols)

col_data <- data.frame(
  sample_id        = sample_names,
  group            = ifelse(grepl("^Ctrl", sample_names), "nonIBD", "IBD"),
  diagnosis_binary = ifelse(grepl("^Ctrl", sample_names), "nonIBD", "IBD"),
  row.names        = sample_names,
  stringsAsFactors = FALSE
)
colnames(mat_imputed) <- sample_names

cat("   Label counts:\n")
print(table(col_data$diagnosis_binary))

# ----------------------------------------------------------
# 9. DIFFERENTIAL EXPRESSION (t-test per protein)
# ----------------------------------------------------------
cat("\n8. Running differential analysis (t-test)...\n")

ibd_samples  <- col_data$sample_id[col_data$diagnosis_binary == "IBD"]
ctrl_samples <- col_data$sample_id[col_data$diagnosis_binary == "nonIBD"]

de_results <- data.frame(
  protein   = rownames(mat_imputed),
  mean_IBD  = rowMeans(mat_imputed[, ibd_samples,  drop = FALSE], na.rm = TRUE),
  mean_ctrl = rowMeans(mat_imputed[, ctrl_samples, drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)
de_results$log2FC <- de_results$mean_IBD - de_results$mean_ctrl

pvals <- apply(mat_imputed, 1, function(row) {
  ibd_vals  <- as.numeric(row[ibd_samples])
  ctrl_vals <- as.numeric(row[ctrl_samples])
  if (sum(!is.na(ibd_vals)) < 2 || sum(!is.na(ctrl_vals)) < 2) return(NA)
  tryCatch(t.test(ibd_vals, ctrl_vals)$p.value, error = function(e) NA)
})

de_results$pval <- pvals
de_results$padj <- p.adjust(pvals, method = "BH")
de_results <- de_results %>% arrange(padj)

sig_proteins <- de_results %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1)
cat("   Significant proteins:", nrow(sig_proteins), "\n")
cat("   Top 10:\n")
print(head(sig_proteins[, c("protein","log2FC","padj")], 10))

# ----------------------------------------------------------
# 10. SELECT TOP FEATURES FOR ML
# ----------------------------------------------------------
cat("\n9. Selecting top 100 features for ML...\n")
top_proteins <- de_results %>%
  filter(!is.na(padj)) %>%
  head(100) %>%
  pull(protein)

mat_ml <- t(mat_imputed[top_proteins, ])

# ----------------------------------------------------------
# 11. SAVE
# ----------------------------------------------------------
cat("\n10. Saving outputs...\n")
saveRDS(list(
  normalized_matrix = mat_imputed,
  ml_matrix         = mat_ml,
  col_data          = col_data,
  de_results        = de_results,
  sig_proteins      = sig_proteins
), file.path(dest, "data/proteomics_processed.rds"))

write.csv(de_results,   file.path(dest, "data/proteomics_DE_results.csv"),   row.names = FALSE)
write.csv(sig_proteins, file.path(dest, "data/proteomics_sig_proteins.csv"), row.names = FALSE)

cat("   Saved: data/proteomics_processed.rds\n")
cat("\n=== DONE ===\n")
