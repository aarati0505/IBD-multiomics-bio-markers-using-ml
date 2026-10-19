# ================================================================
# Export full annotated list of 294 concordant gene-protein pairs
# Adds full protein names and descriptions via biomaRt
# Output: figures/concordance/concordant_294_annotated.csv
# ================================================================

library(dplyr)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "figures", "concordance")

# ----------------------------------------------------------------
# 1. Load concordance table and filter to concordant only
# ----------------------------------------------------------------
conc <- read.csv(file.path(OUT_DIR, "concordance_table_both_sig.csv"),
                 stringsAsFactors = FALSE)

concordant <- conc %>%
  filter(concordant == TRUE) %>%
  arrange(rna_padj) %>%
  mutate(
    rank          = row_number(),
    rna_lfc       = round(rna_lfc,  3),
    prot_lfc      = round(prot_lfc, 3),
    rna_padj      = signif(rna_padj,  3),
    prot_padj     = signif(prot_padj, 3)
  )

cat("Concordant pairs:", nrow(concordant), "\n")

# ----------------------------------------------------------------
# 2. Try to get full protein names via biomaRt
# ----------------------------------------------------------------
protein_names <- NULL

tryCatch({
  library(biomaRt)
  cat("Connecting to Ensembl biomaRt...\n")
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  anno <- getBM(
    attributes = c("hgnc_symbol",
                   "description",
                   "gene_biotype",
                   "chromosome_name"),
    filters    = "hgnc_symbol",
    values     = concordant$gene,
    mart       = mart
  )
  # Clean description (remove source tag at end)
  anno$description <- sub("\\s*\\[.*\\]$", "", anno$description)
  anno <- anno %>%
    filter(!duplicated(hgnc_symbol)) %>%
    rename(gene = hgnc_symbol,
           protein_full_name = description,
           gene_type = gene_biotype,
           chromosome = chromosome_name)

  protein_names <- anno
  cat("  Retrieved annotations for", nrow(anno), "genes\n")

}, error = function(e) {
  cat("  biomaRt not available or connection failed:", conditionMessage(e), "\n")
  cat("  Will use gene symbols only.\n")
})

# ----------------------------------------------------------------
# 3. Build final table
# ----------------------------------------------------------------
if (!is.null(protein_names)) {
  final <- concordant %>%
    left_join(protein_names, by = "gene") %>%
    select(
      rank,
      gene_symbol      = gene,
      protein_full_name,
      direction,
      mRNA_log2FC      = rna_lfc,
      protein_log2FC   = prot_lfc,
      mRNA_padj        = rna_padj,
      protein_padj     = prot_padj,
      chromosome,
      gene_type
    )
} else {
  # Fallback: manual protein name lookup for top hits
  known_names <- data.frame(
    gene = c("UQCRB","IST1","LZIC","FERMT1","HSPD1","EDEM2","ROCK2",
             "HDAC2","PRPF40A","GSPT1","NFIA","DRAP1","ROCK1","HMGCS2",
             "LDHD","COL7A1","HSPB1","VIM","ACTB","GAPDH","ENO1",
             "PKM","ALDOA","TPI1","PGK1","GPI","LDHA","MDH2","IDH2",
             "CS","ACO2","SDHA","SDHB","FH","MDH1","PC","PCCA","PCCB"),
    protein_full_name = c(
      "Cytochrome b-c1 complex subunit 7",
      "IST1 homolog — ESCRT-III component",
      "Leucine zipper inhibitor",
      "Fermitin family member 1 (Kindlin-1) — cell-matrix adhesion",
      "Heat shock protein 60 — mitochondrial chaperone",
      "ER degradation-enhancing alpha-mannosidase-like protein 2",
      "Rho-associated protein kinase 2 — cytoskeleton regulation",
      "Histone deacetylase 2 — chromatin remodelling",
      "Pre-mRNA-processing factor 40 homolog A",
      "Eukaryotic peptide chain release factor GTP-binding subunit",
      "Nuclear factor I/A — transcription factor",
      "Dr1-associated protein — transcription repressor",
      "Rho-associated protein kinase 1",
      "Hydroxymethylglutaryl-CoA synthase — cholesterol synthesis",
      "Lactate dehydrogenase D — D-lactate metabolism",
      "Collagen type VII alpha 1 chain — basement membrane",
      "Heat shock protein beta-1",
      "Vimentin — intermediate filament structural protein",
      "Actin beta — cytoskeletal protein",
      "Glyceraldehyde-3-phosphate dehydrogenase — glycolysis",
      "Alpha-enolase — glycolysis",
      "Pyruvate kinase — glycolysis",
      "Fructose-bisphosphate aldolase A — glycolysis",
      "Triosephosphate isomerase — glycolysis",
      "Phosphoglycerate kinase 1 — glycolysis",
      "Glucose-6-phosphate isomerase — glycolysis",
      "Lactate dehydrogenase A — glycolysis",
      "Malate dehydrogenase 2 — TCA cycle",
      "Isocitrate dehydrogenase 2 — TCA cycle",
      "Citrate synthase — TCA cycle",
      "Aconitase 2 — TCA cycle",
      "Succinate dehydrogenase complex subunit A — TCA cycle",
      "Succinate dehydrogenase complex subunit B — TCA cycle",
      "Fumarate hydratase — TCA cycle",
      "Malate dehydrogenase 1 — TCA cycle",
      "Pyruvate carboxylase — gluconeogenesis",
      "Propionyl-CoA carboxylase alpha — amino acid catabolism",
      "Propionyl-CoA carboxylase beta — amino acid catabolism"
    ),
    stringsAsFactors = FALSE
  )

  final <- concordant %>%
    left_join(known_names, by = "gene") %>%
    mutate(protein_full_name = ifelse(
      is.na(protein_full_name),
      paste0(gene, " protein (see UniProt: ", gene, ")"),
      protein_full_name
    )) %>%
    select(
      rank,
      gene_symbol    = gene,
      protein_full_name,
      direction,
      mRNA_log2FC    = rna_lfc,
      protein_log2FC = prot_lfc,
      mRNA_padj      = rna_padj,
      protein_padj   = prot_padj
    )
}

# ----------------------------------------------------------------
# 4. Save
# ----------------------------------------------------------------
write.csv(final,
          file.path(OUT_DIR, "concordant_294_annotated.csv"),
          row.names = FALSE)

cat("\nSaved: concordant_294_annotated.csv\n")
cat("Total genes:", nrow(final), "\n\n")

# Print full list
cat("=== ALL", nrow(final), "CONCORDANT GENE-PROTEIN PAIRS ===\n\n")
for (i in seq_len(nrow(final))) {
  r <- final[i, ]
  cat(sprintf("%3d. %-12s | %-55s | %s | mRNA LFC: %+.3f | Prot LFC: %+.3f\n",
              r$rank,
              r$gene_symbol,
              substr(r$protein_full_name, 1, 55),
              r$direction,
              r$mRNA_log2FC,
              r$protein_log2FC))
}
