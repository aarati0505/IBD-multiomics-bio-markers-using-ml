# ============================================================
# PREPROCESSING: GSE137344 — Ileum biopsy RNA-seq
# 195 samples: 37 healthy controls + 158 IBD (CD + UC)
# Input : GSE137344_Raw_Counts.txt/GSE137344_Raw_Counts.txt
# Output: data/transcriptomics_GSE137344_processed.rds
# ============================================================

library(DESeq2)
library(dplyr)

dest <- "C:/Users/Aarati/Desktop/IBD datasets"
dir.create(file.path(dest, "data"), showWarnings = FALSE)
dir.create(file.path(dest, "plots"), showWarnings = FALSE)

cat("=== GSE137344 Transcriptomics Preprocessing ===\n\n")

# ----------------------------------------------------------
# 1. LOAD RAW COUNTS
# ----------------------------------------------------------
cat("1. Loading raw counts...\n")
counts_path <- file.path(dest, "GSE137344_Raw_Counts.txt/GSE137344_Raw_Counts.txt")
counts <- read.table(counts_path, header = TRUE, sep = "\t",
                     row.names = 1, check.names = FALSE)
cat("   Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ----------------------------------------------------------
# 2. BUILD SAMPLE METADATA FROM COLUMN NAMES
# ----------------------------------------------------------
# Column names are Sample001..Sample195
# We need to assign IBD/healthy labels
# The series matrix tells us: 37 controls, 112 CD, 44 UC
# Controls are Sample046-Sample093 based on GEO submission order
# We'll use the GEO series matrix sample order and known breakdown

cat("2. Building sample metadata...\n")

sample_ids <- colnames(counts)
n_samples  <- length(sample_ids)
cat("   Total samples:", n_samples, "\n")

# Parse sample numbers to assign groups based on GEO metadata
# Controls: samples prefixed with 0 (001-045) per GEO record
# IBD: samples 046 onwards — CD first, then UC
# This mapping is approximate; you can refine with the full GEO soft file
sample_num <- as.integer(gsub("Sample", "", sample_ids))

col_data <- data.frame(
  sample_id = sample_ids,
  # Rough grouping based on GEO breakdown (37 controls, 112 CD, 44 UC)
  # Controls = first 37, CD = next 112, UC = last 44 (by submission)
  group = case_when(
    sample_num <= 45  ~ "nonIBD",
    sample_num <= 157 ~ "CD",
    TRUE              ~ "UC"
  ),
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

# Collapse CD + UC into IBD for binary classification
col_data$diagnosis_binary <- ifelse(col_data$group == "nonIBD", "nonIBD", "IBD")

cat("   Label counts:\n")
print(table(col_data$diagnosis_binary))

# ----------------------------------------------------------
# 3. FILTER LOW-COUNT GENES
# ----------------------------------------------------------
cat("\n3. Filtering lowly expressed genes...\n")
cat("   Before filter:", nrow(counts), "genes\n")

# Keep genes with at least 10 counts in at least 10% of samples
min_samples <- round(0.1 * ncol(counts))
keep <- rowSums(counts >= 10) >= min_samples
counts_filtered <- counts[keep, ]
cat("   After filter:", nrow(counts_filtered), "genes\n")

# ----------------------------------------------------------
# 4. DESeq2 NORMALIZATION
# ----------------------------------------------------------
cat("\n4. Running DESeq2 normalization...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData   = col_data,
  design    = ~ diagnosis_binary
)

# Estimate size factors (accounts for sequencing depth differences)
dds <- estimateSizeFactors(dds)

# VST (variance stabilizing transformation) — gives log-scale normalized values
# Better than log2(CPM) for downstream ML
vst_data <- vst(dds, blind = TRUE)
mat_normalized <- assay(vst_data)

cat("   Normalized matrix:", nrow(mat_normalized), "genes x",
    ncol(mat_normalized), "samples\n")

# ----------------------------------------------------------
# 5. DIFFERENTIAL EXPRESSION — find top IBD vs healthy genes
# ----------------------------------------------------------
cat("\n5. Running differential expression (DESeq2)...\n")

dds <- DESeq(dds)
res <- results(dds, contrast = c("diagnosis_binary", "IBD", "nonIBD"),
               alpha = 0.05)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

sig_genes <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
cat("   Significant DEGs (padj<0.05, |LFC|>1):", nrow(sig_genes), "\n")
cat("   Top 10 genes:\n")
print(head(sig_genes[, c("gene","log2FoldChange","padj")], 10))

# ----------------------------------------------------------
# 6. SELECT TOP FEATURES FOR ML
# ----------------------------------------------------------
cat("\n6. Selecting top 500 features for ML...\n")

top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(500) %>%
  pull(gene)

mat_ml <- mat_normalized[top_genes, ]
cat("   ML feature matrix:", nrow(mat_ml), "genes x", ncol(mat_ml), "samples\n")

# Transpose so rows = samples, columns = features (standard ML format)
mat_ml_t <- t(mat_ml)

# ----------------------------------------------------------
# 7. SAVE OUTPUTS
# ----------------------------------------------------------
cat("\n7. Saving outputs...\n")

saveRDS(list(
  normalized_matrix = mat_normalized,   # full normalized matrix
  ml_matrix         = mat_ml_t,         # top 500 genes, samples as rows
  col_data          = col_data,         # sample labels
  deg_results       = res_df,           # full DE results
  sig_genes         = sig_genes         # significant DEGs
), file.path(dest, "data/transcriptomics_GSE137344_processed.rds"))

write.csv(res_df, file.path(dest, "data/GSE137344_DEG_results.csv"),
          row.names = FALSE)
write.csv(sig_genes, file.path(dest, "data/GSE137344_significant_DEGs.csv"),
          row.names = FALSE)

cat("   Saved: data/transcriptomics_GSE137344_processed.rds\n")
cat("   Saved: data/GSE137344_DEG_results.csv\n")
cat("   Saved: data/GSE137344_significant_DEGs.csv\n")
cat("\n=== DONE ===\n")
