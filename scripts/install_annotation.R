# Install org.Hs.eg.db — the standard human gene annotation package
# This is a local database — no internet needed after install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
cat("Done.\n")
