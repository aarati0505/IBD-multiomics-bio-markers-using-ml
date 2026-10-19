library(curl)
options(timeout = 600)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "proteomics", "PXD001608_colon_biopsy")
dir.create(OUT_DIR, showWarnings = FALSE)

h <- new_handle()
handle_setopt(h, ssl_verifypeer = 0, ssl_verifyhost = 0)

# PXD001608 — colon mucosa IBD vs healthy, label-free LC-MS/MS
# Files on PRIDE FTP: ftp.pride.ebi.ac.uk/pride/data/archive/2016/09/PXD001608/
# The processed result file is the MaxQuant output submitted by the authors

files <- list(
  list(
    name = "proteinGroups_PXD001608.txt",
    url  = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/09/PXD001608/proteinGroups.txt"
  ),
  list(
    name = "peptides_PXD001608.txt",
    url  = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/09/PXD001608/peptides.txt"
  ),
  list(
    name = "parameters_PXD001608.txt",
    url  = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/09/PXD001608/parameters.txt"
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
    cat("  OK —", round(file.size(outfile)/1e6, 1), "MB\n")
  }, error = function(e) cat("  FAILED:", conditionMessage(e), "\n"))
}

cat("\nDone.\n")
