# ================================================================
# TRANSCRIPTOMICS PREPROCESSING — GSE186507
# Whole blood RNA-seq | 1030 samples | CD + UC + healthy
# MSCCR cohort — Mount Sinai
# Input  : GSE186507_MSCCR_Blood_adjcounts.txt (log-normalized, batch-corrected)
# Outputs: processed/GSE186507/GSE186507_normalized.rds
#          processed/GSE186507/GSE186507_top500_ml_matrix.csv
#          processed/GSE186507/GSE186507_DEG_results.csv
# ================================================================

library(limma)
library(dplyr)
library(data.table)
library(GEOquery)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "GSE186507")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
options(timeout = 300)

cat("================================================================\n")
cat("GSE186507 — Blood RNA-seq Transcriptomics Preprocessing\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD ADJUSTED COUNTS
# ----------------------------------------------------------------
cat("Step 1: Loading adjusted count matrix...\n")

adj_file <- file.path(BASE,
  "transcriptomics/GSE186507_blood_MSCCR/GSE186507_MSCCR_Blood_adjcounts.txt")

mat <- fread(adj_file, sep = " ", header = TRUE, data.table = FALSE)
rownames(mat) <- mat[[1]]
mat <- as.matrix(mat[, -1])

cat("  Genes:", nrow(mat), "| Samples:", ncol(mat), "\n")
cat("  Value range: [", round(min(mat),2), ",", round(max(mat),2), "]\n")

# ----------------------------------------------------------------
# 2. GET DIAGNOSIS LABELS — USE CACHED GEO METADATA IF AVAILABLE
# ----------------------------------------------------------------
cat("\nStep 2: Fetching diagnosis labels from GEO...\n")

cached_pdata <- file.path(OUT_DIR, "GSE186507_GEO_pdata.csv")

if (file.exists(cached_pdata)) {
  cat("  Using cached GEO metadata from previous run\n")
  pdata <- read.csv(cached_pdata, row.names = 1, stringsAsFactors = FALSE)
} else {
  library(GEOquery)
  gse   <- getGEO("GSE186507", GSEMatrix = TRUE, getGPL = FALSE)
  pdata <- pData(gse[[1]])
  write.csv(pdata, cached_pdata)
}

cat("  Rows in GEO metadata:", nrow(pdata), "\n")

# The diagnosis column confirmed from run output: "ibd_disease:ch1"
diag_col <- "ibd_disease:ch1"
if (!diag_col %in% colnames(pdata)) {
  # fallback — try diseasebi:ch1
  diag_col <- "diseasebi:ch1"
}
if (!diag_col %in% colnames(pdata)) {
  # last resort — grep
  diag_col <- grep("disease|diagno|ibd", colnames(pdata), ignore.case=TRUE, value=TRUE)[1]
}

cat("  Using diagnosis column:", diag_col, "\n")
cat("  Unique values:\n")
print(table(pdata[[diag_col]]))

# Save GEO metadata
write.csv(pdata, file.path(OUT_DIR, "GSE186507_GEO_pdata.csv"))

# ----------------------------------------------------------------
# 3. BUILD COL_DATA — match GEO rows to matrix columns
#    GEO pdata rownames are GSM accessions
#    Matrix colnames are MSCCR_reGRID_N
#    Linking via geo_accession or title field
# ----------------------------------------------------------------
cat("\nStep 3: Building sample metadata...\n")

# GEO title field should contain MSCCR_reGRID_N
if ("title" %in% colnames(pdata)) {
  pdata$msccr_id <- gsub("\\s.*", "", pdata$title)   # take first word
} else {
  pdata$msccr_id <- rownames(pdata)
}

# Standardize diagnosis labels
pdata$diagnosis <- pdata[[diag_col]]
pdata$diagnosis_binary <- case_when(
  pdata$diagnosis %in% c("UC","CD","IBD","Crohn's Disease",
                          "Ulcerative Colitis","crohn") ~ "IBD",
  pdata$diagnosis %in% c("Control","nonIBD","healthy","HC",
                          "non-IBD","control") ~ "nonIBD",
  TRUE ~ "UNKNOWN"
)

# Match to matrix column order
mat_samples <- colnames(mat)
idx         <- match(mat_samples, pdata$msccr_id)

col_data <- data.frame(
  sample_id        = mat_samples,
  diagnosis        = ifelse(is.na(idx), "UNKNOWN", pdata$diagnosis[idx]),
  diagnosis_binary = ifelse(is.na(idx), "UNKNOWN", pdata$diagnosis_binary[idx]),
  row.names        = mat_samples,
  stringsAsFactors = FALSE
)

cat("  Label summary:\n")
print(table(col_data$diagnosis_binary))

# If most are still UNKNOWN, GEO title didn't match — use positional order
if (sum(col_data$diagnosis_binary == "UNKNOWN") > 0.5 * nrow(col_data)) {
  cat("  Title matching failed — using positional order from GEO pdata\n")
  if (nrow(pdata) == ncol(mat)) {
    col_data$diagnosis        <- pdata$diagnosis
    col_data$diagnosis_binary <- pdata$diagnosis_binary
    cat("  Positional match label summary:\n")
    print(table(col_data$diagnosis_binary))
  }
}

# ----------------------------------------------------------------
# 4. FILTER LOW-VARIANCE GENES
# ----------------------------------------------------------------
cat("\nStep 4: Filtering low-variance genes...\n")
cat("  Before:", nrow(mat), "genes\n")

gene_vars    <- apply(mat, 1, var)
mat_filtered <- mat[gene_vars > quantile(gene_vars, 0.10), ]
cat("  After:", nrow(mat_filtered), "genes\n")

# ----------------------------------------------------------------
# 5. LIMMA DIFFERENTIAL EXPRESSION
# ----------------------------------------------------------------
cat("\nStep 5: Running limma differential expression...\n")

keep_samp <- col_data$diagnosis_binary %in% c("IBD", "nonIBD")
cat("  Samples used for DE:", sum(keep_samp),
    "(IBD:", sum(col_data$diagnosis_binary=="IBD" & keep_samp),
    "| nonIBD:", sum(col_data$diagnosis_binary=="nonIBD" & keep_samp), ")\n")

mat_de <- mat_filtered[, keep_samp]
groups <- factor(col_data$diagnosis_binary[keep_samp],
                 levels = c("nonIBD", "IBD"))

design <- model.matrix(~ groups)
fit    <- lmFit(mat_de, design)
fit    <- eBayes(fit)
de_res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_res <- tibble::rownames_to_column(de_res, "gene")

sig <- de_res %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
cat("  Significant DEGs (adj.P<0.05, |LFC|>1):", nrow(sig), "\n")
cat("  Top 10:\n")
print(head(sig[, c("gene","logFC","adj.P.Val")], 10))

# ----------------------------------------------------------------
# 6. SELECT TOP 500 FOR ML
# ----------------------------------------------------------------
cat("\nStep 6: Selecting top 500 features for ML...\n")

top500  <- head(de_res$gene, 500)
ml_mat  <- t(mat_filtered[top500, ])
cat("  ML matrix:", nrow(ml_mat), "samples x", ncol(ml_mat), "genes\n")

# ----------------------------------------------------------------
# 7. SAVE
# ----------------------------------------------------------------
cat("\nStep 7: Saving outputs...\n")

saveRDS(list(
  normalized_matrix = mat_filtered,
  ml_matrix         = ml_mat,
  col_data          = col_data,
  deg_results       = de_res,
  sig_genes         = sig,
  dataset_info      = list(
    name      = "GSE186507",
    tissue    = "whole_blood",
    n_genes   = nrow(mat_filtered),
    n_samples = ncol(mat_filtered)
  )
), file.path(OUT_DIR, "GSE186507_normalized.rds"))

write.csv(col_data,
          file.path(OUT_DIR, "GSE186507_col_data.csv"))
write.csv(de_res,
          file.path(OUT_DIR, "GSE186507_DEG_results.csv"), row.names = FALSE)
write.csv(sig,
          file.path(OUT_DIR, "GSE186507_significant_DEGs.csv"), row.names = FALSE)
write.csv(as.data.frame(ml_mat),
          file.path(OUT_DIR, "GSE186507_top500_ml_matrix.csv"))

cat("  Saved: GSE186507_normalized.rds\n")
cat("  Saved: GSE186507_col_data.csv\n")
cat("  Saved: GSE186507_DEG_results.csv\n")
cat("  Saved: GSE186507_significant_DEGs.csv\n")
cat("  Saved: GSE186507_top500_ml_matrix.csv\n")
cat("\n=== GSE186507 DONE ===\n")
