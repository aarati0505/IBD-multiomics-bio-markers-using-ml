# Try installing with SSL verification disabled
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Method 1: disable SSL for BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Try with httr config to bypass SSL
old_opt <- getOption("download.file.extra")
options(download.file.extra = "--insecure")

tryCatch({
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE,
                       checkBuilt = FALSE)
  cat("org.Hs.eg.db installed successfully\n")
}, error = function(e) {
  cat("BiocManager failed:", conditionMessage(e), "\n")
  cat("Trying mygene package instead...\n")
  install.packages("mygene", repos = "https://cloud.r-project.org")
})

options(download.file.extra = old_opt)
cat("Done.\n")
