library(biomformat)
library(dplyr)

cat("Loading BIOM file...\n")
biom_file <- "C:/Users/Aarati/Desktop/IBD datasets/metabolomics/HMP2_metabolomics.biom.txt"
b <- read_biom(biom_file)

cat("Shape:", nrow(b), "x", ncol(b), "\n")
cat("Sample IDs (first 5):", paste(colnames(b)[1:5], collapse=", "), "\n")
cat("Feature IDs (first 5):", paste(rownames(b)[1:5], collapse=", "), "\n")

# Observation metadata
obs_meta <- observation_metadata(b)
cat("\nObservation metadata class:", class(obs_meta), "\n")
if (!is.null(obs_meta)) {
  if (is.data.frame(obs_meta)) {
    cat("Metadata columns:", paste(names(obs_meta), collapse=", "), "\n")
    print(head(obs_meta, 3))
  } else {
    cat("Type:", typeof(obs_meta), "\n")
    cat("First entry:", paste(obs_meta[[1]], collapse=", "), "\n")
  }
}

# Sample metadata
samp_meta <- sample_metadata(b)
cat("\nSample metadata class:", class(samp_meta), "\n")
if (!is.null(samp_meta)) {
  if (is.data.frame(samp_meta)) {
    cat("Columns:", paste(names(samp_meta)[1:min(10,ncol(samp_meta))], collapse=", "), "\n")
    print(head(samp_meta[, 1:min(5,ncol(samp_meta))], 3))
  }
}

# Data sparsity
cat("\nChecking sparsity on small subset...\n")
mat_small <- as.matrix(biom_data(b))[1:200, 1:20]
zero_pct <- round(mean(mat_small == 0) * 100, 1)
cat("Zero rate:", zero_pct, "%\n")
cat("Value range:", range(mat_small), "\n")
cat("Non-zero values sample:", head(mat_small[mat_small > 0], 10), "\n")
