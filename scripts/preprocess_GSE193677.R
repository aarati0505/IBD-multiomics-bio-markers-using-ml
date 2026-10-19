# ================================================================
# TRANSCRIPTOMICS PREPROCESSING — GSE193677
# Intestinal biopsy RNA-seq | 2490 samples | CD + UC + healthy
# MSCCR cohort — Mount Sinai — PRIMARY TRAINING DATASET
# Input  : GSE193677_MSCCR_Biopsy_adjcounts.txt (log-normalized, batch-corrected)
# Outputs: processed/GSE193677/GSE193677_normalized.rds
#          processed/GSE193677/GSE193677_top500_ml_matrix.csv
#          processed/GSE193677/GSE193677_DEG_results.csv
# ================================================================

library(limma)
library(dplyr)
library(data.table)
library(GEOquery)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "processed", "GSE193677")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
options(timeout = 300)

cat("================================================================\n")
cat("GSE193677 — Intestinal Biopsy Transcriptomics Preprocessing\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD ADJUSTED COUNTS (~640MB — use fread)
# ----------------------------------------------------------------
cat("Step 1: Loading adjusted count matrix (~640MB)...\n")
cat("  Please wait 2-3 minutes...\n")

adj_file <- file.path(BASE,
  "transcriptomics/GSE193677_biopsy_MSCCR/GSE193677_MSCCR_Biopsy_adjcounts.txt/GSE193677_MSCCR_Biopsy_adjcounts.txt")

mat <- fread(adj_file, sep = " ", header = TRUE, data.table = FALSE)
rownames(mat) <- mat[[1]]
mat <- as.matrix(mat[, -1])

cat("  Genes:", nrow(mat), "| Samples:", ncol(mat), "\n")
cat("  Value range: [", round(min(mat),2), ",", round(max(mat),2), "]\n")

# ----------------------------------------------------------------
# 2. GET DIAGNOSIS LABELS FROM GEO
# ----------------------------------------------------------------
cat("\nStep 2: Fetching diagnosis labels from GEO...\n")

gse   <- getGEO("GSE193677", GSEMatrix = TRUE, getGPL = FALSE)
pdata <- pData(gse[[1]])

cat("  Rows in GEO metadata:", nrow(pdata), "\n")

# Find diagnosis column
diag_col <- grep("ibd_disease|disease|diagno|condition|phenotype",
                 colnames(pdata), ignore.case = TRUE, value = TRUE)
diag_col <- diag_col[!grepl("^status$|submission|update|platform",
                              diag_col, ignore.case = TRUE)]
cat("  Diagnosis candidate columns:", paste(diag_col, collapse=", "), "\n")

if (length(diag_col) > 0) {
  best_col <- diag_col[1]
  cat("  Using:", best_col, "\n")
  cat("  Unique values:\n")
  print(table(pdata[[best_col]]))
} else {
  best_col <- NULL
}

write.csv(pdata, file.path(OUT_DIR, "GSE193677_GEO_pdata.csv"))
cat("  Saved GEO metadata to GSE193677_GEO_pdata.csv\n")

# ----------------------------------------------------------------
# 3. BUILD COL_DATA
# ----------------------------------------------------------------
cat("\nStep 3: Building sample metadata...\n")

sample_ids <- colnames(mat)
n          <- length(sample_ids)

if (!is.null(best_col)) {
  # Match via GEO title
  if ("title" %in% colnames(pdata)) {
    pdata$title_clean <- trimws(pdata$title)
    idx <- match(sample_ids, pdata$title_clean)
  } else {
    idx <- seq_len(min(n, nrow(pdata)))
  }

  raw_diag <- ifelse(is.na(idx), "UNKNOWN", pdata[[best_col]][idx])

  # If title match failed, try positional (GEO rows = matrix columns order)
  if (sum(raw_diag == "UNKNOWN") > 0.5 * n && nrow(pdata) == n) {
    cat("  Title match failed — using positional order\n")
    raw_diag <- pdata[[best_col]]
  }
} else {
  raw_diag <- rep("UNKNOWN", n)
}

col_data <- data.frame(
  sample_id = sample_ids,
  diagnosis = raw_diag,
  diagnosis_binary = case_when(
    raw_diag %in% c("UC","CD","IBD","Crohn's Disease",
                    "Ulcerative Colitis","crohn","uc","cd") ~ "IBD",
    raw_diag %in% c("Control","nonIBD","healthy","HC",
                    "non-IBD","control","Control") ~ "nonIBD",
    TRUE ~ "UNKNOWN"
  ),
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

cat("  Label summary:\n")
print(table(col_data$diagnosis_binary))

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
    "(IBD:", sum(col_data$diagnosis_binary == "IBD" & keep_samp),
    "| nonIBD:", sum(col_data$diagnosis_binary == "nonIBD" & keep_samp), ")\n")

if (sum(keep_samp) < 10) {
  cat("  Not enough labeled samples for DE — skipping\n")
  de_res <- data.frame()
  sig    <- data.frame()
  ml_mat <- t(mat_filtered)
} else {
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

  top500 <- head(de_res$gene, 500)
  ml_mat <- t(mat_filtered[top500, ])
  cat("  ML matrix:", nrow(ml_mat), "samples x", ncol(ml_mat), "genes\n")

  write.csv(de_res, file.path(OUT_DIR, "GSE193677_DEG_results.csv"),     row.names = FALSE)
  write.csv(sig,    file.path(OUT_DIR, "GSE193677_significant_DEGs.csv"), row.names = FALSE)
  write.csv(as.data.frame(ml_mat),
            file.path(OUT_DIR, "GSE193677_top500_ml_matrix.csv"))
}

# ----------------------------------------------------------------
# 6. SAVE
# ----------------------------------------------------------------
cat("\nStep 6: Saving outputs...\n")

saveRDS(list(
  normalized_matrix = mat_filtered,
  ml_matrix         = ml_mat,
  col_data          = col_data,
  deg_results       = de_res,
  sig_genes         = sig,
  dataset_info      = list(
    name      = "GSE193677",
    tissue    = "intestinal_biopsy",
    n_genes   = nrow(mat_filtered),
    n_samples = ncol(mat_filtered)
  )
), file.path(OUT_DIR, "GSE193677_normalized.rds"))

write.csv(col_data, file.path(OUT_DIR, "GSE193677_col_data.csv"))

cat("  Saved: GSE193677_normalized.rds\n")
cat("  Saved: GSE193677_col_data.csv\n")
cat("\n=== GSE193677 DONE ===\n")
