# ================================================================
# HMP2 Metabolomics BIOM — Fast Parser (byte-position approach)
# File structure confirmed:
#   Bytes 0-247   : JSON header with shape/format info
#   Bytes ~248    : "data":[[r,c,v],...] starts
#   Byte 525307144: "rows":[...] starts
#   End           : "columns":[...] ends
# ================================================================

library(data.table)
library(Matrix)
library(dplyr)

BASE      <- "C:/Users/Aarati/Desktop/IBD datasets"
BIOM_PATH <- file.path(BASE, "metabolomics/HMP2_metabolomics.biom.txt")
META_PATH <- file.path(BASE, "data/hmp2_metadata_clean.csv")
OUT_DIR   <- file.path(BASE, "processed", "HMP2_metabolomics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ROWS_BYTE <- 525307144L   # confirmed byte position of "rows":[
COLS_BYTE <- 528597160L   # confirmed byte position of "columns":[

cat("================================================================\n")
cat("HMP2 Metabolomics BIOM — Byte-Position Parser\n")
cat("File size:", round(file.info(BIOM_PATH)$size / 1e6, 1), "MB\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. EXTRACT ROW IDs — from bytes after ROWS_BYTE
# ----------------------------------------------------------------
cat("Step 1: Extracting feature IDs (rows section)...\n")
t0 <- proc.time()
con <- file(BIOM_PATH, "rb")
seek(con, ROWS_BYTE)
rows_size  <- COLS_BYTE - ROWS_BYTE     # exact rows section size
rows_chunk <- readChar(con, rows_size)
close(con)

row_ids <- regmatches(rows_chunk,
  gregexpr('"id":\\s*"([^"]+)"', rows_chunk, perl = TRUE))[[1]]
row_ids <- gsub('"id":\\s*"', '', row_ids)
row_ids <- gsub('"$', '', row_ids)
cat("  Features:", length(row_ids), "\n")
cat("  First 5:", paste(row_ids[1:5], collapse=", "), "\n")
rm(rows_chunk); gc()

# ----------------------------------------------------------------
# 2. EXTRACT COLUMN IDs — from COLS_BYTE to end
# ----------------------------------------------------------------
cat("\nStep 2: Extracting sample IDs (columns section)...\n")
file_size <- file.info(BIOM_PATH)$size
cols_size <- file_size - COLS_BYTE
con2 <- file(BIOM_PATH, "rb")
seek(con2, COLS_BYTE)
cols_chunk <- readChar(con2, cols_size)
close(con2)

col_ids <- regmatches(cols_chunk,
  gregexpr('"id":\\s*"([^"]+)"', cols_chunk, perl = TRUE))[[1]]
col_ids <- gsub('"id":\\s*"', '', col_ids)
col_ids <- gsub('"$', '', col_ids)
cat("  Samples:", length(col_ids), "\n")
cat("  First 5:", paste(col_ids[1:5], collapse=", "), "\n")
rm(cols_chunk); gc()

# ----------------------------------------------------------------
# 3. FIND DATA START POSITION
# ----------------------------------------------------------------
cat("\nStep 3: Finding data section start...\n")
con3   <- file(BIOM_PATH, "rb")
header <- readChar(con3, 500)
close(con3)

# Pattern handles optional space: "data": [[ or "data":[[
data_match  <- regexpr('"data":\\s*\\[\\[', header, perl = TRUE)
data_start  <- data_match + attr(data_match, "match.length")  # byte after [[
data_end    <- ROWS_BYTE - 3L   # just before ],"rows":
data_length <- data_end - data_start
cat("  Data starts at byte:", data_start, "\n")
cat("  Data ends at byte:", data_end, "\n")
cat("  Data section size:", round(data_length/1e6, 1), "MB\n")

# ----------------------------------------------------------------
# 4. READ DATA SECTION IN CHUNKS AND WRITE TO TEMP CSV
# ----------------------------------------------------------------
cat("\nStep 4: Streaming data section to temp CSV...\n")
t1 <- proc.time()

tmp_csv <- file.path(OUT_DIR, "biom_data_temp.csv")
con_in  <- file(BIOM_PATH, "rb")
seek(con_in, data_start)

con_out <- file(tmp_csv, "w")
cat("row,col,val\n", file = con_out)

chunk_size <- 20000000L   # 20MB per chunk
bytes_read <- 0L
leftover   <- ""

repeat {
  to_read <- min(chunk_size, data_length - bytes_read)
  if (to_read <= 0) break
  chunk      <- readChar(con_in, to_read)
  bytes_read <- bytes_read + nchar(chunk, type = "bytes")
  chunk      <- paste0(leftover, chunk)

  # Find last complete triplet ]]
  last_bracket <- max(gregexpr("\\]", chunk, fixed = TRUE)[[1]])
  if (last_bracket < 1) { leftover <- chunk; next }

  process_part <- substring(chunk, 1, last_bracket)
  leftover     <- substring(chunk, last_bracket + 1)

  # Convert [r,c,v] to r,c,v lines
  # Order matters: replace ],[ first, then remaining [ and ]
  process_part <- gsub("],[",  "\n", process_part, fixed = TRUE)
  process_part <- gsub("[",    "",   process_part, fixed = TRUE)
  process_part <- gsub("]",    "",   process_part, fixed = TRUE)
  # Remove any blank lines or trailing commas
  process_part <- gsub(",\n",  "\n", process_part, fixed = TRUE)
  process_part <- gsub("^\n",  "",   process_part)

  cat(process_part, file = con_out)

  pct <- round(bytes_read / data_length * 100)
  if (pct %% 10 == 0) cat("  ...", pct, "%\n")
}
close(con_in)
close(con_out)
cat("  Stream complete in", round((proc.time()-t1)["elapsed"], 1), "sec\n")

# Verify temp CSV looks correct before running fread
tmp_lines <- readLines(tmp_csv, n = 5)
cat("  Temp CSV first 5 lines:\n")
cat(paste(tmp_lines, collapse="\n"), "\n")
tmp_size <- round(file.info(tmp_csv)$size / 1e6, 1)
cat("  Temp CSV size:", tmp_size, "MB\n")

# ----------------------------------------------------------------
# 5. READ CSV WITH fread
# ----------------------------------------------------------------
cat("\nStep 5: Reading CSV with fread...\n")
t2  <- proc.time()
dt  <- fread(tmp_csv, header = TRUE, sep = ",",
             colClasses = c("integer","integer","numeric"),
             showProgress = FALSE)
# Remove any NA rows (blank lines)
dt  <- dt[!is.na(row) & !is.na(col)]
dt[, row := row + 1L]
dt[, col := col + 1L]
file.remove(tmp_csv)
cat("  fread done in", round((proc.time()-t2)["elapsed"], 1), "sec\n")
cat("  Non-zero entries:", nrow(dt), "\n")
cat("  Sparsity:",
    round((1 - nrow(dt)/(length(row_ids)*length(col_ids)))*100, 1), "%\n")

# ----------------------------------------------------------------
# 6. BUILD SPARSE MATRIX
# ----------------------------------------------------------------
cat("\nStep 6: Building sparse matrix...\n")
sp_mat <- sparseMatrix(
  i    = dt$row,
  j    = dt$col,
  x    = dt$val,
  dims = c(length(row_ids), length(col_ids)),
  dimnames = list(row_ids, col_ids)
)
rm(dt); gc()
cat("  Matrix:", nrow(sp_mat), "x", ncol(sp_mat), "\n")

# ----------------------------------------------------------------
# 7. ATTACH METADATA
# ----------------------------------------------------------------
cat("\nStep 7: Attaching metadata...\n")
meta     <- read.csv(META_PATH, stringsAsFactors = FALSE, check.names = FALSE)
mbx_meta <- meta %>%
  filter(data_type == "metabolomics") %>%
  distinct(`External ID`, .keep_all = TRUE)

matched  <- mbx_meta[match(col_ids, mbx_meta[["External ID"]]), ]
has_meta <- !is.na(matched[["External ID"]])
cat("  Samples with metadata:", sum(has_meta), "/", length(col_ids), "\n")

col_data <- data.frame(
  sample_id        = col_ids,
  participant_id   = matched[["Participant ID"]],
  diagnosis        = matched[["diagnosis"]],
  diagnosis_binary = matched[["diagnosis_binary"]],
  week_num         = matched[["week_num"]],
  sex              = matched[["sex"]],
  site_name        = matched[["site_name"]],
  row.names        = col_ids,
  stringsAsFactors = FALSE
)
valid_samples <- col_ids[has_meta]
col_data      <- col_data[has_meta, ]
sp_mat        <- sp_mat[, valid_samples, drop = FALSE]
cat("  Retained:", ncol(sp_mat), "samples\n")
cat("  Diagnosis breakdown:\n")
print(table(col_data$diagnosis_binary))

# ----------------------------------------------------------------
# 8. FILTER BY DETECTION RATE (>10% in at least one group)
# ----------------------------------------------------------------
cat("\nStep 8: Filtering features...\n")
cat("  Before:", nrow(sp_mat), "\n")
ibd_s  <- col_data$sample_id[col_data$diagnosis_binary == "IBD"]
ctrl_s <- col_data$sample_id[col_data$diagnosis_binary == "nonIBD"]

det_ibd  <- rowMeans(sp_mat[, ibd_s,  drop = FALSE] != 0)
det_ctrl <- rowMeans(sp_mat[, ctrl_s, drop = FALSE] != 0)
sp_mat   <- sp_mat[(det_ibd >= 0.10) | (det_ctrl >= 0.10), , drop = FALSE]
cat("  After:", nrow(sp_mat), "\n")

cat("\nConverting to dense matrix...\n")
mat <- as.matrix(sp_mat)
rm(sp_mat); gc()

# ----------------------------------------------------------------
# 9. NORMALISE: zero->NA, TSS, log2, impute
# ----------------------------------------------------------------
cat("\nStep 9: Normalising...\n")
mat[mat == 0] <- NA
cat("  Missing:", round(mean(is.na(mat))*100,1), "%\n")

mat_nz   <- mat; mat_nz[is.na(mat_nz)] <- 0
mat_tss  <- sweep(mat_nz, 2, colSums(mat_nz), "/") * 1e6
mat_log2 <- log2(mat_tss + 1)
mat_log2[mat_log2 == 0 & is.na(mat)] <- NA
rm(mat, mat_nz, mat_tss); gc()
cat("  Range (non-NA):", round(range(mat_log2, na.rm=TRUE), 2), "\n")

global_mean <- mean(mat_log2, na.rm=TRUE)
global_sd   <- sd(mat_log2,   na.rm=TRUE)
impute_val  <- global_mean - 1.8 * global_sd
cat("  Imputing NAs:", round(impute_val, 3), "\n")
mat_log2[is.na(mat_log2)] <- impute_val

# ----------------------------------------------------------------
# 10. DIFFERENTIAL — Wilcoxon
# ----------------------------------------------------------------
cat("\nStep 10: Wilcoxon test (IBD vs nonIBD)...\n")
mean_ibd  <- rowMeans(mat_log2[, ibd_s,  drop = FALSE])
mean_ctrl <- rowMeans(mat_log2[, ctrl_s, drop = FALSE])

pvals <- apply(mat_log2, 1, function(r) {
  iv <- r[ibd_s]; cv <- r[ctrl_s]
  if (sum(iv != impute_val) < 5 || sum(cv != impute_val) < 5) return(NA)
  tryCatch(wilcox.test(iv, cv, exact=FALSE)$p.value, error=function(e) NA)
})

de <- data.frame(feature=rownames(mat_log2), mean_IBD=mean_ibd,
                 mean_ctrl=mean_ctrl, log2FC=mean_ibd-mean_ctrl,
                 pval=pvals, stringsAsFactors=FALSE)
de$padj <- p.adjust(de$pval, method="BH")
de       <- de[order(de$padj), ]

sig <- de[!is.na(de$padj) & de$padj < 0.05 & abs(de$log2FC) > 0.5, ]
cat("  Significant (padj<0.05, |LFC|>0.5):", nrow(sig), "\n")
cat("  Top 10:\n")
print(head(sig[, c("feature","log2FC","pval","padj")], 10))

# ----------------------------------------------------------------
# 11. ML MATRIX + SAVE
# ----------------------------------------------------------------
cat("\nStep 11: ML matrix (top 100)...\n")
top100 <- head(de[!is.na(de$padj), "feature"], 100)
ml_mat <- t(mat_log2[top100, ])
cat("  ML matrix:", nrow(ml_mat), "x", ncol(ml_mat), "\n")

cat("\nStep 12: Saving...\n")
saveRDS(list(
  normalized_matrix = mat_log2, ml_matrix = ml_mat,
  col_data = col_data, de_results = de, sig_features = sig,
  dataset_info = list(
    name="HMP2_metabolomics", n_feat=nrow(mat_log2),
    n_samp=ncol(mat_log2),
    n_IBD=sum(col_data$diagnosis_binary=="IBD"),
    n_nonIBD=sum(col_data$diagnosis_binary=="nonIBD"),
    n_sig=nrow(sig)
  )
), file.path(OUT_DIR, "HMP2_metabolomics_processed.rds"))

write.csv(as.data.frame(ml_mat),
          file.path(OUT_DIR, "HMP2_metabolomics_ml_matrix.csv"))
write.csv(de,  file.path(OUT_DIR, "HMP2_metabolomics_DE_results.csv"),  row.names=FALSE)
write.csv(sig, file.path(OUT_DIR, "HMP2_metabolomics_sig_features.csv"), row.names=FALSE)
write.csv(col_data, file.path(OUT_DIR, "HMP2_metabolomics_col_data.csv"))

cat("  Saved all outputs to:", OUT_DIR, "\n")
cat("\n=== DONE ===\n")
cat("Features:", nrow(mat_log2), "| Samples:", ncol(mat_log2),
    "| Sig:", nrow(sig), "\n")
