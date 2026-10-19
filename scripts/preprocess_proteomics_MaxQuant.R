# ================================================================
# PROTEOMICS PREPROCESSING — PXD001608 (UC vs Healthy Colon Biopsy)
# Bennike et al. 2015 — Neutrophil Extracellular Traps in UC
# Input  : proteomics/PXD001608_UC_colon_biopsy/proteinGroups.txt
# Outputs: processed/proteomics_MaxQuant/
#            proteomics_MaxQuant_processed.rds
#            proteomics_MaxQuant_ml_matrix.csv
#            proteomics_MaxQuant_DE_results.csv
#            proteomics_MaxQuant_sig_proteins.csv
#            proteomics_MaxQuant_col_data.csv
# ================================================================

library(dplyr)
library(data.table)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "proteomics_MaxQuant")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("Proteomics Preprocessing — PXD001608 (UC vs Healthy Colon Biopsy)\n")
cat("Bennike et al. 2015 — Neutrophil Extracellular Traps in UC\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD proteinGroups.txt
#    Row 1 = column headers
#    Row 2 = #!{Type} annotation row (Perseus format) — must skip
#    Row 3+ = data
# ----------------------------------------------------------------
cat("Step 1: Loading proteinGroups.txt...\n")

pg_file <- file.path(BASE,
  "proteomics/PXD001608_UC_colon_biopsy/proteinGroups.txt")

# Read raw lines to detect Perseus type annotation row
raw_lines <- readLines(pg_file, n = 3)
cat("  Header preview:", substr(raw_lines[1], 1, 80), "\n")
cat("  Row 2 preview:", substr(raw_lines[2], 1, 40), "\n")

has_type_row <- grepl("^#!\\{Type\\}", raw_lines[2])
cat("  Perseus type annotation row detected:", has_type_row, "\n")

# fread reads header from row 1, so to skip row 2 we filter afterwards
pg <- fread(pg_file, sep = "\t", header = TRUE, data.table = FALSE)

# Drop the #!{Type} annotation row if present
if (has_type_row) {
  pg <- pg[!grepl("^#!\\{Type\\}", pg[[1]]), ]
}

cat("  Rows:", nrow(pg), "| Columns:", ncol(pg), "\n")

# ----------------------------------------------------------------
# 2. REMOVE CONTAMINANTS AND DECOYS
#    Column names: "Contaminant", "Reverse", "Only identified by site"
#    (Note: this dataset uses "Contaminant", not "Potential contaminant")
# ----------------------------------------------------------------
cat("\nStep 2: Removing contaminants and decoys...\n")
cat("  Before:", nrow(pg), "\n")

cont_col <- grep("^[Cc]ontaminant$|^Potential contaminant$",
                 colnames(pg), value = TRUE)
rev_col  <- grep("^Reverse$",               colnames(pg), value = TRUE)
site_col <- grep("^Only identified by site", colnames(pg), value = TRUE)

cat("  Using contaminant col:", paste(cont_col, collapse=", "), "\n")
cat("  Using reverse col:",     paste(rev_col,  collapse=", "), "\n")
cat("  Using site col:",        paste(site_col, collapse=", "), "\n")

if (length(cont_col) > 0) pg <- pg[pg[[cont_col[1]]] != "+", ]
if (length(rev_col)  > 0) pg <- pg[pg[[rev_col[1]]]  != "+", ]
if (length(site_col) > 0) pg <- pg[pg[[site_col[1]]] != "+", ]

cat("  After:", nrow(pg), "\n")

# ----------------------------------------------------------------
# 3. EXTRACT GENE NAMES from Protein IDs (FASTA GN= tag)
#    Protein IDs column contains FASTA entries like:
#    >sp|P01850|TRBC1_HUMAN ... GN=TRBC1 PE=1 SV=3
#    We parse the first GN= hit; fall back to UniProt accession.
# ----------------------------------------------------------------
cat("\nStep 3: Extracting gene names from Protein IDs...\n")

protein_ids <- pg[["Protein IDs"]]

parse_gene_name <- function(fasta_str) {
  # Try GN= tag first
  gn <- regmatches(fasta_str, regexpr("GN=([A-Za-z0-9_\\-\\.]+)", fasta_str))
  if (length(gn) > 0 && nchar(gn) > 0) {
    return(sub("GN=", "", gn))
  }
  # Fall back: UniProt accession from |ACC|
  acc <- regmatches(fasta_str, regexpr("\\|([A-Z0-9]+)\\|", fasta_str))
  if (length(acc) > 0 && nchar(acc) > 0) {
    return(gsub("\\|", "", acc))
  }
  # Last resort: first 20 chars
  return(substr(fasta_str, 1, 20))
}

gene_names_raw <- sapply(protein_ids, parse_gene_name, USE.NAMES = FALSE)
gene_names     <- make.unique(gene_names_raw)

cat("  Parsed", length(gene_names), "gene names\n")
cat("  Examples:", paste(head(gene_names, 5), collapse=", "), "\n")
cat("  Unique names:", length(unique(gene_names_raw)), "\n")

# ----------------------------------------------------------------
# 4. EXTRACT iBAQ INTENSITY COLUMNS
#    Use "iBAQ <sample>" columns (exclude the bare "iBAQ" sum column)
# ----------------------------------------------------------------
cat("\nStep 4: Extracting iBAQ intensity columns...\n")

ibaq_cols <- grep("^iBAQ [A-Za-z]", colnames(pg), value = TRUE)

if (length(ibaq_cols) == 0) {
  # Fall back to per-sample Intensity columns
  ibaq_cols <- grep("^Intensity [A-Za-z]", colnames(pg), value = TRUE)
  intensity_type <- "Intensity"
  cat("  No iBAQ sample cols found — falling back to Intensity columns\n")
} else {
  intensity_type <- "iBAQ"
  cat("  Using iBAQ columns:", length(ibaq_cols), "samples\n")
}

sample_names <- gsub(paste0("^", intensity_type, " "), "", ibaq_cols)
cat("  Sample names (first 5):", paste(sample_names[1:5], collapse=", "), "\n")
cat("  Total samples:", length(sample_names), "\n")

# ----------------------------------------------------------------
# 5. BUILD NUMERIC MATRIX
# ----------------------------------------------------------------
cat("\nStep 5: Building numeric matrix...\n")

mat_df <- as.data.frame(lapply(
  pg[, ibaq_cols, drop = FALSE],
  function(x) as.numeric(as.character(x))
))
rownames(mat_df) <- gene_names
colnames(mat_df) <- sample_names
mat <- as.matrix(mat_df)

cat("  Matrix dimensions:", nrow(mat), "x", ncol(mat), "\n")

# ----------------------------------------------------------------
# 6. BUILD SAMPLE METADATA
#    UC_<patient>_<replicate>   -> IBD
#    Ctrl_<patient>_<replicate> -> nonIBD
# ----------------------------------------------------------------
cat("\nStep 6: Building sample metadata...\n")

col_data <- data.frame(
  sample_id        = sample_names,
  patient_id       = sub("_\\d+$", "", sample_names),        # UC_8 / Ctrl_1
  replicate        = sub(".*_(\\d+)$", "\\1", sample_names), # 1 / 2 / 3
  diagnosis        = ifelse(grepl("^UC_", sample_names), "UC", "Ctrl"),
  diagnosis_binary = ifelse(grepl("^UC_", sample_names), "IBD", "nonIBD"),
  row.names        = sample_names,
  stringsAsFactors = FALSE
)

cat("  Sample breakdown:\n")
print(table(col_data$diagnosis_binary))
cat("  Unique patients:", length(unique(col_data$patient_id)), "\n")

# ----------------------------------------------------------------
# 7. REPLACE ZEROS WITH NA (0 = "not detected" in MaxQuant)
# ----------------------------------------------------------------
cat("\nStep 7: Replacing zeros with NA...\n")
mat[mat == 0] <- NA
missing_pct <- round(mean(is.na(mat)) * 100, 1)
cat("  Missing (NA) after zero replacement:", missing_pct, "%\n")

# ----------------------------------------------------------------
# 8. FILTER: keep proteins detected in >50% of samples
#    Apply within each group to avoid group-specific proteins
#    being lost, then take the union.
# ----------------------------------------------------------------
cat("\nStep 8: Filtering by detection rate (>50% in at least one group)...\n")
cat("  Before:", nrow(mat), "proteins\n")

ibd_idx  <- col_data$sample_id[col_data$diagnosis_binary == "IBD"]
ctrl_idx <- col_data$sample_id[col_data$diagnosis_binary == "nonIBD"]

det_ibd  <- rowMeans(!is.na(mat[, ibd_idx,  drop = FALSE]))
det_ctrl <- rowMeans(!is.na(mat[, ctrl_idx, drop = FALSE]))

keep <- (det_ibd >= 0.5) | (det_ctrl >= 0.5)
mat  <- mat[keep, ]
cat("  After:", nrow(mat), "proteins\n")

# ----------------------------------------------------------------
# 9. LOG2 TRANSFORM
# ----------------------------------------------------------------
cat("\nStep 9: Log2 transforming iBAQ intensities...\n")
mat_log2 <- log2(mat)     # zeros already NA; no +1 needed for iBAQ
cat("  Range (non-NA):", round(range(mat_log2, na.rm=TRUE), 2), "\n")

# ----------------------------------------------------------------
# 10. IMPUTE MISSING VALUES
#     Strategy: per-protein min-value imputation (shifted -1.8 SDs)
#     Represents "below-detection-limit" proteins.
# ----------------------------------------------------------------
cat("\nStep 10: Imputing missing values (min-value - 1.8 sd)...\n")

global_mean <- mean(mat_log2, na.rm = TRUE)
global_sd   <- sd(mat_log2,   na.rm = TRUE)
impute_val  <- global_mean - 1.8 * global_sd

cat("  Global mean:", round(global_mean, 2), "| SD:", round(global_sd, 2), "\n")
cat("  Imputation value:", round(impute_val, 2), "\n")

mat_imputed <- mat_log2
mat_imputed[is.na(mat_imputed)] <- impute_val
cat("  Any remaining NA:", any(is.na(mat_imputed)), "\n")

# ----------------------------------------------------------------
# 11. DIFFERENTIAL EXPRESSION — Welch t-test per protein
# ----------------------------------------------------------------
cat("\nStep 11: Running differential analysis (Welch t-test, IBD vs nonIBD)...\n")

de_results <- data.frame(
  protein   = rownames(mat_imputed),
  mean_IBD  = rowMeans(mat_imputed[, ibd_idx,  drop = FALSE]),
  mean_ctrl = rowMeans(mat_imputed[, ctrl_idx, drop = FALSE]),
  stringsAsFactors = FALSE
)
de_results$log2FC <- de_results$mean_IBD - de_results$mean_ctrl

pvals <- apply(mat_imputed, 1, function(row) {
  ibd_v  <- as.numeric(row[ibd_idx])
  ctrl_v <- as.numeric(row[ctrl_idx])
  n_ibd  <- sum(!is.na(ibd_v))
  n_ctrl <- sum(!is.na(ctrl_v))
  if (n_ibd < 2 || n_ctrl < 2) return(NA)
  tryCatch(t.test(ibd_v, ctrl_v, var.equal = FALSE)$p.value,
           error = function(e) NA)
})

de_results$pval <- pvals
de_results$padj <- p.adjust(pvals, method = "BH")
de_results       <- de_results %>% arrange(padj)

sig_proteins <- de_results %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1)

cat("  Significant proteins (padj < 0.05, |log2FC| > 1):", nrow(sig_proteins), "\n")
if (nrow(sig_proteins) > 0) {
  cat("  Top 10 hits:\n")
  print(head(sig_proteins[, c("protein", "log2FC", "pval", "padj")], 10))
}

# ----------------------------------------------------------------
# 12. SELECT TOP 100 PROTEINS FOR ML
#     Take the top 100 by adjusted p-value (irrespective of LFC
#     threshold) so the ML matrix always has 100 features.
# ----------------------------------------------------------------
cat("\nStep 12: Selecting top 100 proteins for ML matrix...\n")

top100 <- de_results %>%
  filter(!is.na(padj)) %>%
  slice_head(n = 100) %>%
  pull(protein)

ml_mat <- t(mat_imputed[top100, ])   # samples x proteins
cat("  ML matrix:", nrow(ml_mat), "samples x", ncol(ml_mat), "proteins\n")

# ----------------------------------------------------------------
# 13. SAVE OUTPUTS
# ----------------------------------------------------------------
cat("\nStep 13: Saving outputs to", OUT_DIR, "...\n")

saveRDS(list(
  normalized_matrix = mat_imputed,
  ml_matrix         = ml_mat,
  col_data          = col_data,
  de_results        = de_results,
  sig_proteins      = sig_proteins,
  dataset_info      = list(
    name       = "PXD001608",
    tissue     = "colon_biopsy",
    condition  = "UC_vs_Ctrl",
    n_IBD      = sum(col_data$diagnosis_binary == "IBD"),
    n_nonIBD   = sum(col_data$diagnosis_binary == "nonIBD"),
    n_proteins = nrow(mat_imputed),
    n_sig      = nrow(sig_proteins)
  )
), file.path(OUT_DIR, "proteomics_MaxQuant_processed.rds"))

write.csv(as.data.frame(ml_mat),
          file.path(OUT_DIR, "proteomics_MaxQuant_ml_matrix.csv"))
write.csv(de_results,
          file.path(OUT_DIR, "proteomics_MaxQuant_DE_results.csv"),
          row.names = FALSE)
write.csv(sig_proteins,
          file.path(OUT_DIR, "proteomics_MaxQuant_sig_proteins.csv"),
          row.names = FALSE)
write.csv(col_data,
          file.path(OUT_DIR, "proteomics_MaxQuant_col_data.csv"))

cat("  Saved: proteomics_MaxQuant_processed.rds\n")
cat("  Saved: proteomics_MaxQuant_ml_matrix.csv\n")
cat("  Saved: proteomics_MaxQuant_DE_results.csv\n")
cat("  Saved: proteomics_MaxQuant_sig_proteins.csv\n")
cat("  Saved: proteomics_MaxQuant_col_data.csv\n")

cat("\n=== DONE ===\n")
cat("Dataset summary:\n")
cat("  Proteins quantified:", nrow(mat_imputed), "\n")
cat("  Samples:", ncol(mat_imputed), "(IBD:", sum(col_data$diagnosis_binary=="IBD"),
    "| nonIBD:", sum(col_data$diagnosis_binary=="nonIBD"), ")\n")
cat("  Significant DEPs:", nrow(sig_proteins), "\n")
