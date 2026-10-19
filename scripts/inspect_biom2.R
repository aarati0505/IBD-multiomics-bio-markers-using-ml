# Try faster approach: parse BIOM JSON header only to understand structure
# then extract data using direct matrix access

library(jsonlite)

biom_path <- "C:/Users/Aarati/Desktop/IBD datasets/metabolomics/HMP2_metabolomics.biom.txt"
cat("Reading file header...\n")

# Read first 50000 chars to get metadata sections
con  <- file(biom_path, "r")
chunk <- readChar(con, 50000)
close(con)

# Extract shape
shape_match <- regmatches(chunk, regexpr('"shape":\\s*\\[\\d+,\\s*\\d+\\]', chunk))
cat("Shape:", shape_match, "\n")

# Extract first few row IDs
row_ids <- regmatches(chunk, gregexpr('"id":\\s*"C18n_[^"]+"|"id":\\s*"[A-Z][0-9][^"]{3,}"', chunk))[[1]]
cat("First 5 feature IDs:", paste(row_ids[1:min(5,length(row_ids))], collapse=", "), "\n")

# Extract first few column (sample) IDs
col_match <- regmatches(chunk, gregexpr('"columns".*', chunk))
cat("\nLooking for columns section in first 50k chars:", length(col_match[[1]]) > 0, "\n")

# The columns section is at the END of the file - read last 10000 chars
con2 <- file(biom_path, "rb")
file_size <- file.info(biom_path)$size
seek(con2, file_size - 10000)
tail_chunk <- readChar(con2, 10000)
close(con2)

# Extract last few sample IDs from tail
last_ids <- regmatches(tail_chunk, gregexpr('"id":\\s*"[A-Z0-9]+"', tail_chunk))[[1]]
cat("Last 5 sample IDs:", paste(tail(last_ids, 5), collapse=", "), "\n")

# Now try biomformat with timeout
cat("\nNow trying biomformat on smaller w_metadata file...\n")
library(biomformat)
t0 <- proc.time()
b2 <- read_biom("C:/Users/Aarati/Desktop/IBD datasets/metabolomics/HMP2_metabolomics_w_metadata.biom/HMP2_metabolomics_w_metadata.biom")
cat("Load time:", round((proc.time()-t0)["elapsed"], 1), "sec\n")
cat("Shape:", nrow(b2), "x", ncol(b2), "\n")

samp <- sample_metadata(b2)
cat("Sample metadata:", ifelse(is.null(samp), "NULL", "EXISTS"), "\n")
if (!is.null(samp) && is.data.frame(samp)) {
  cat("Columns:", paste(names(samp), collapse=", "), "\n")
}
