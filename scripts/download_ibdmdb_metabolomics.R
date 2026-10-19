library(curl)
options(timeout = 60)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "metabolomics")
dir.create(OUT_DIR, showWarnings = FALSE)

h <- new_handle()
handle_setopt(h, ssl_verifypeer = 0, ssl_verifyhost = 0, nobody = 1, timeout = 20)

# Try every known URL pattern for HMP2 metabolomics
urls <- c(
  # Old product paths
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/1723/HMP2_metabolomics.csv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/HMP2_metabolomics.csv.gz",
  # Without version number
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/HMP2_metabolomics.csv",
  # From get_data.sh script we found earlier
  "https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz",
  # Alternate naming
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/1723/iHMP_metabolomics.pcl.csv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/iHMP_metabolomics.pcl.csv.gz",
  # Try downloads instead of tunnel
  "https://ibdmdb.org/downloads/products/HMP2/MBX/1723/HMP2_metabolomics.csv.gz",
  # Try raw path
  "https://ibdmdb.org/data/HMP2/MBX/HMP2_metabolomics.csv.gz"
)

cat("Testing ibdmdb.org metabolomics URLs...\n\n")

found_url <- NULL
for (u in urls) {
  cat("Testing:", u, "... ")
  tryCatch({
    h2 <- new_handle()
    handle_setopt(h2, ssl_verifypeer = 0, ssl_verifyhost = 0,
                  nobody = 1, timeout = 15)
    r <- curl_fetch_memory(u, handle = h2)
    cat("HTTP", r$status_code, "\n")
    if (r$status_code == 200) {
      found_url <- u
      cat("  *** WORKING URL FOUND ***\n")
      break
    }
  }, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
}

if (!is.null(found_url)) {
  cat("\nDownloading from:", found_url, "\n")
  outfile <- file.path(OUT_DIR, "HMP2_metabolomics.csv.gz")
  h3 <- new_handle()
  handle_setopt(h3, ssl_verifypeer = 0, ssl_verifyhost = 0)
  options(timeout = 600)
  tryCatch({
    curl_download(found_url, destfile = outfile, handle = h3)
    cat("Downloaded:", round(file.size(outfile)/1e6, 1), "MB\n")
  }, error = function(e) cat("Download failed:", conditionMessage(e), "\n"))
} else {
  cat("\nNo working URL found on ibdmdb.org\n")
  cat("Try opening this URL directly in your browser:\n")
  cat("https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz\n")
}
