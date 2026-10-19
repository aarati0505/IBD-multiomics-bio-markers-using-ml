library(dplyr)

BASE <- "C:/Users/Aarati/Desktop/IBD datasets"

deg137 <- read.csv(file.path(BASE, "processed/GSE137344/GSE137344_DEG_results.csv"),
                   stringsAsFactors = FALSE)
allpxd <- read.csv(file.path(BASE, "processed/proteomics_MaxQuant/proteomics_MaxQuant_DE_results.csv"),
                   stringsAsFactors = FALSE)

cat("GSE137344 genes:", nrow(deg137), "\n")
cat("PXD001608 proteins:", nrow(allpxd), "\n")

# Match by gene symbol
matched <- inner_join(
  deg137  %>% select(gene, rna_lfc = log2FoldChange, rna_padj = padj),
  allpxd  %>% select(protein, prot_lfc = log2FC, prot_padj = padj),
  by = c("gene" = "protein")
)
cat("Matched gene-protein pairs:", nrow(matched), "\n")

# Concordance
matched <- matched %>%
  filter(!is.na(rna_lfc), !is.na(prot_lfc)) %>%
  mutate(concordant = (rna_lfc > 0) == (prot_lfc > 0))

cat("Concordant (same direction):", sum(matched$concordant), "\n")
cat("Discordant (opposite direction):", sum(!matched$concordant), "\n")
cat("Concordance rate:", round(mean(matched$concordant) * 100, 1), "%\n")

# Both significant
both_sig <- matched %>%
  filter(!is.na(rna_padj), !is.na(prot_padj), rna_padj < 0.05, prot_padj < 0.05)
cat("\nBoth sig (padj < 0.05):", nrow(both_sig), "\n")
if (nrow(both_sig) > 0) {
  cat("Both sig concordant:", sum(both_sig$concordant), "\n")
  cat("Both sig discordant:", sum(!both_sig$concordant), "\n")
  cat("Top pairs:\n")
  print(both_sig %>% arrange(rna_padj) %>%
        select(gene, rna_lfc, prot_lfc, rna_padj, prot_padj, concordant) %>%
        head(20))
}
