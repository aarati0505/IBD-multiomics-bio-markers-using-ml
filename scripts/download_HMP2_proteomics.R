library(curl)
options(timeout = 600)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "proteomics", "HMP2_stool_metaproteomics")
dir.create(OUT_DIR, showWarnings = FALSE)

h <- new_handle()
handle_setopt(h, ssl_verifypeer = 0, ssl_verifyhost = 0)

# HMP2 MPX (metaproteomics) — stool protein counts
# File naming from hmp2_analysis GitHub repo: iHMP_all_proteomics_*pep-pro_1p_FDR_*ppm.pcl.tsv.gz
# Try known ibdmdb.org tunnel paths for MPX data

urls_to_try <- c(
  "https://ibdmdb.org/tunnel/products/HMP2/MPX/1543/hmp2_MPX_not_normalized_ecs.names.pcl.tsv",
  "https://ibdmdb.org/tunnel/products/HMP2/MPX/1543/hmp2_MPX_not_normalized_kos.names.pcl.tsv",
  "https://ibdmdb.org/tunnel/products/HMP2/MPX/1543/iHMP_proteomics.pcl.tsv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/Proteomics/1543/HMP2_proteomics_PCAs_table.csv"
)

cat("Testing ibdmdb.org MPX URLs...\n\n")
for (u in urls_to_try) {
  cat("Testing:", basename(u), "...")
  tryCatch({
    h2 <- new_handle()
    handle_setopt(h2, ssl_verifypeer = 0, ssl_verifyhost = 0, nobody = 1, timeout = 15)
    r <- curl_fetch_memory(u, handle = h2)
    cat(" HTTP", r$status_code, "\n")
    if (r$status_code == 200) {
      cat("  *** FOUND — downloading now ***\n")
      outfile <- file.path(OUT_DIR, basename(u))
      h3 <- new_handle()
      handle_setopt(h3, ssl_verifypeer = 0, ssl_verifyhost = 0)
      curl_download(u, destfile = outfile, handle = h3)
      cat("  Saved:", round(file.size(outfile)/1e6, 1), "MB\n")
    }
  }, error = function(e) cat(" FAILED:", conditionMessage(e), "\n"))
}

cat("\nDone.\n")
