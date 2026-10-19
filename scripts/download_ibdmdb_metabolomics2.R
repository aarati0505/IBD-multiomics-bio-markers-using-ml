library(curl)
options(timeout = 600)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "metabolomics")
dir.create(OUT_DIR, showWarnings = FALSE)

h <- new_handle()
handle_setopt(h, ssl_verifypeer = 0, ssl_verifyhost = 0)

# URLs visible from the browser address bar on that page
files <- list(
  list(
    name = "HMP2_metabolomics.biom",
    url  = "https://g-22ca.190ebd.5bc.data.globus.org/ibdmdb/products/HMP2/MBX/HMP2_metabolomics.biom"
  ),
  list(
    name = "HMP2_metabolomics_w_metadata.biom.gz",
    url  = "https://g-22ca.190ebd.5bc.data.globus.org/ibdmdb/products/HMP2/MBX/HMP2_metabolomics_w_metadata.biom.gz"
  ),
  list(
    name = "C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx",
    url  = "https://g-22ca.190ebd.5bc.data.globus.org/ibdmdb/products/HMP2/MBX/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"
  )
)

for (f in files) {
  outfile <- file.path(OUT_DIR, f$name)
  if (file.exists(outfile) && file.size(outfile) > 1000) {
    cat("Already exists:", f$name, "\n")
    next
  }
  cat("Downloading:", f$name, "...\n")
  tryCatch({
    curl_download(f$url, destfile = outfile, handle = h)
    cat("  OK —", round(file.size(outfile)/1e6, 2), "MB\n")
  }, error = function(e) cat("  FAILED:", conditionMessage(e), "\n"))
}

cat("\nDone.\n")
