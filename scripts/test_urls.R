library(curl)

urls <- c(
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/1723/HMP2_metabolomics.csv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/MPX/1543/HMP2_proteomics_PCAs_table.csv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/HTX/1750/host_tx_counts.tsv.gz",
  "https://ibdmdb.org/tunnel/products/HMP2/MBX/1723/",
  "https://ibdmdb.org/tunnel/products/HMP2/"
)

for (u in urls) {
  cat("Testing:", u, "\n")
  tryCatch({
    h <- new_handle()
    handle_setopt(h, ssl_verifypeer = 0, ssl_verifyhost = 0, nobody = 1)
    r <- curl_fetch_memory(u, handle = h)
    cat("  HTTP", r$status_code, "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
