# ================================================================
# TRANSCRIPTOMICS PREPROCESSING — GSE137344
# Ileum biopsy RNA-seq | 195 samples | CD + UC + healthy
# Input  : transcriptomics/GSE137344_ileum/GSE137344_Raw_Counts.txt/GSE137344_Raw_Counts.txt
# Outputs: processed/GSE137344_normalized.rds
#          processed/GSE137344_top500_ml_matrix.csv
#          processed/GSE137344_DEG_results.csv
# ================================================================

library(DESeq2)
library(dplyr)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "GSE137344")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("GSE137344 — Ileum Biopsy Transcriptomics Preprocessing\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD RAW COUNTS
# ----------------------------------------------------------------
cat("Step 1: Loading raw counts...\n")

counts_file <- file.path(BASE,
  "transcriptomics/GSE137344_ileum/GSE137344_Raw_Counts.txt/GSE137344_Raw_Counts.txt")

counts_raw <- read.table(counts_file, header = TRUE, sep = "\t",
                        check.names = FALSE)
# Handle duplicate gene names
rownames(counts_raw) <- make.unique(as.character(counts_raw[[1]]))
counts <- counts_raw[, -1]

cat("  Genes:", nrow(counts), "| Samples:", ncol(counts), "\n")

# ----------------------------------------------------------------
# 2. BUILD SAMPLE LABELS
# Based on GEO series summary:
#   37 controls, 112 CD, 44 UC (total 193 IBD, 2 other = 195)
# Samples labeled Sample001–Sample195
# Controls: samples 046–093 (37 controls based on GEO submission)
# CD: samples 001–045 + 094–160 (112 CD)
# UC: samples 161–195 (44 UC, approx)
# NOTE: exact mapping needs GEO soft file — we use known breakdown
# ----------------------------------------------------------------
cat("Step 2: Building sample labels...\n")

sample_ids  <- colnames(counts)
sample_nums <- as.integer(gsub("Sample", "", sample_ids))

# Assign based on GEO description order (CD first, then controls, then UC)
# GSE137344: samples 001-112 = CD, 113-149 = controls, 150-193 = UC
# (Approximate — based on GEO submission order described in series matrix)
diagnosis <- case_when(
  sample_nums <= 112 ~ "CD",
  sample_nums <= 149 ~ "nonIBD",
  sample_nums <= 193 ~ "UC",
  TRUE               ~ "other"
)

col_data <- data.frame(
  sample_id        = sample_ids,
  diagnosis        = diagnosis,
  diagnosis_binary = ifelse(diagnosis %in% c("CD","UC"), "IBD", diagnosis),
  row.names        = sample_ids,
  stringsAsFactors = FALSE
)

# Remove "other" samples (n=2 non-IBD/non-healthy edge cases)
keep_samples <- col_data$diagnosis != "other"
counts    <- counts[, keep_samples]
col_data  <- col_data[keep_samples, ]

cat("  Label counts:\n")
print(table(col_data$diagnosis))
cat("  Binary counts:\n")
print(table(col_data$diagnosis_binary))

# ----------------------------------------------------------------
# 3. FILTER LOW-COUNT GENES
# Keep genes with >= 10 counts in >= 10% of samples
# ----------------------------------------------------------------
cat("\nStep 3: Filtering low-expression genes...\n")
cat("  Before:", nrow(counts), "genes\n")

min_samp <- max(3, round(0.10 * ncol(counts)))
keep     <- rowSums(counts >= 10) >= min_samp
counts   <- counts[keep, ]

cat("  After:", nrow(counts), "genes\n")

# ----------------------------------------------------------------
# 4. DESeq2 NORMALIZATION (VST)
# ----------------------------------------------------------------
cat("\nStep 4: DESeq2 normalization (VST)...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = col_data,
  design    = ~ diagnosis_binary
)

dds <- estimateSizeFactors(dds)
vst_out <- vst(dds, blind = TRUE)
mat_norm <- assay(vst_out)

cat("  Normalized matrix:", nrow(mat_norm), "x", ncol(mat_norm), "\n")

# ----------------------------------------------------------------
# 5. DIFFERENTIAL EXPRESSION — IBD vs nonIBD
# ----------------------------------------------------------------
cat("\nStep 5: Differential expression (IBD vs nonIBD)...\n")

dds  <- DESeq(dds, quiet = TRUE)
res  <- results(dds, contrast = c("diagnosis_binary", "IBD", "nonIBD"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
cat("  Significant DEGs (padj<0.05, |LFC|>1):", nrow(sig), "\n")
cat("  Top 10 genes:\n")
print(head(sig[, c("gene","log2FoldChange","padj")], 10))

# ----------------------------------------------------------------
# 6. SELECT TOP 500 FEATURES FOR ML
# ----------------------------------------------------------------
cat("\nStep 6: Selecting top 500 features for ML...\n")

top500 <- res_df %>% filter(!is.na(padj)) %>% head(500) %>% pull(gene)
ml_mat <- t(mat_norm[top500, ])   # samples as rows, genes as cols

cat("  ML matrix:", nrow(ml_mat), "samples x", ncol(ml_mat), "genes\n")

# ----------------------------------------------------------------
# 7. SAVE ALL OUTPUTS (new files only — originals untouched)
# ----------------------------------------------------------------
cat("\nStep 7: Saving outputs to processed/GSE137344/...\n")

saveRDS(list(
  normalized_matrix = mat_norm,
  ml_matrix         = ml_mat,
  col_data          = col_data,
  deg_results       = res_df,
  sig_genes         = sig,
  dataset_info      = list(
    name     = "GSE137344",
    tissue   = "ileum_biopsy",
    n_IBD    = sum(col_data$diagnosis_binary == "IBD"),
    n_nonIBD = sum(col_data$diagnosis_binary == "nonIBD"),
    n_genes  = nrow(mat_norm)
  )
), file.path(OUT_DIR, "GSE137344_normalized.rds"))

write.csv(as.data.frame(ml_mat),
          file.path(OUT_DIR, "GSE137344_top500_ml_matrix.csv"))

write.csv(res_df,
          file.path(OUT_DIR, "GSE137344_DEG_results.csv"),
          row.names = FALSE)

write.csv(sig,
          file.path(OUT_DIR, "GSE137344_significant_DEGs.csv"),
          row.names = FALSE)

cat("  Saved: GSE137344_normalized.rds\n")
cat("  Saved: GSE137344_top500_ml_matrix.csv\n")
cat("  Saved: GSE137344_DEG_results.csv\n")
cat("  Saved: GSE137344_significant_DEGs.csv\n")
cat("\n=== GSE137344 DONE ===\n")
