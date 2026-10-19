options(timeout = 600)
library(curl)

dest     <- "C:/Users/Aarati/Desktop/IBD datasets"
data_dir <- file.path(dest, "data")
dir.create(data_dir, showWarnings = FALSE)

# curl handle with SSL verification disabled (fixes expired cert issues)
h <- new_handle()
handle_setopt(h, ssl_verifypeer = 0, ssl_verifyhost = 0)

github_base <- "https://raw.githubusercontent.com/biobakery/omnibus-and-maaslin2-rscripts-and-hmp2-data/master"

files <- list(
  # Already downloaded — will be skipped
  list(
    name = "HMP2_taxonomic_profiles_week0.csv",
    url  = file.path(github_base, "taxonomic_profiles_pcl_week0.csv"),
    desc = "Metagenomics — species relative abundances (MetaPhlAn), baseline"
  ),
  list(
    name = "HMP2_dna_pathabundance_week0.csv",
    url  = file.path(github_base, "dna_pathabundance_relab_pcl_week0.csv"),
    desc = "Metagenomics — pathway abundances from DNA (HUMAnN), baseline"
  ),
  # Metabolomics — correct path found from original HMP2 get_data.sh script
  list(
    name = "HMP2_metabolomics.csv.gz",
    url  = "https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz",
    desc = "Metabolomics — fecal metabolite abundances (LC-MS)"
  ),
  # Metadata
  list(
    name = "hmp2_metadata.csv",
    url  = "https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv",
    desc = "HMP2 master metadata — links all samples to diagnosis"
  )
)

for (f in files) {
  outfile <- file.path(data_dir, f$name)

  # Skip if already downloaded successfully
  if (file.exists(outfile) && file.size(outfile) > 1000) {
    cat("Already exists:", f$name, "-", round(file.size(outfile)/1e6, 2), "MB — skipping\n\n")
    next
  }

  cat("Downloading:", f$name, "\n")
  cat(" ", f$desc, "\n")
  tryCatch({
    curl_download(f$url, destfile = outfile, handle = h)
    cat("  OK —", round(file.size(outfile) / 1e6, 2), "MB\n")
  }, error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n")
  })
  cat("\n")
}

cat("Done.\n")
