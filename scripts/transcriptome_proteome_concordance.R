# ================================================================
# TRANSCRIPT-PROTEIN CONCORDANCE ANALYSIS
# Do the genes upregulated in transcriptomics code for the
# proteins found upregulated in proteomics?
#
# Approach:
#   - GSE137344 (ileum, gene symbols) vs PXD001608 (tissue proteomics)
#   - Direct symbol match possible — same namespace
#   - GSE186507 / GSE193677 use Ensembl IDs → convert via biomaRt
#     then cross-reference
#
# Outputs:
#   figures/concordance/concordance_heatmap.pdf
#   figures/concordance/concordance_scatter.pdf
#   figures/concordance/concordance_table.csv
# ================================================================

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "figures", "concordance")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("Transcript-Protein Concordance Analysis\n")
cat("GSE137344 (ileum mRNA) vs PXD001608 (colon tissue proteomics)\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# 1. LOAD DATA
# ----------------------------------------------------------------
cat("Step 1: Loading data...\n")

deg137 <- read.csv(file.path(BASE,
           "processed/GSE137344/GSE137344_DEG_results.csv"),
           stringsAsFactors = FALSE)
sig137 <- read.csv(file.path(BASE,
           "processed/GSE137344/GSE137344_significant_DEGs.csv"),
           stringsAsFactors = FALSE)

allpxd <- read.csv(file.path(BASE,
           "processed/proteomics_MaxQuant/proteomics_MaxQuant_DE_results.csv"),
           stringsAsFactors = FALSE)
sigpxd <- read.csv(file.path(BASE,
           "processed/proteomics_MaxQuant/proteomics_MaxQuant_sig_proteins.csv"),
           stringsAsFactors = FALSE)

rds137 <- readRDS(file.path(BASE,
           "processed/GSE137344/GSE137344_normalized.rds"))
pxdrds <- readRDS(file.path(BASE,
           "processed/proteomics_MaxQuant/proteomics_MaxQuant_processed.rds"))

cat("  GSE137344 genes:", nrow(deg137), "| Sig DEGs:", nrow(sig137), "\n")
cat("  PXD001608 proteins:", nrow(allpxd), "| Sig proteins:", nrow(sigpxd), "\n")

# ----------------------------------------------------------------
# 2. MATCH BY GENE SYMBOL
# ----------------------------------------------------------------
cat("\nStep 2: Matching genes to proteins by gene symbol...\n")

all_matched <- inner_join(
  deg137 %>% select(gene,
                    rna_lfc  = log2FoldChange,
                    rna_padj = padj),
  allpxd %>% select(protein,
                    prot_lfc  = log2FC,
                    prot_padj = padj),
  by = c("gene" = "protein")
) %>% filter(!is.na(rna_lfc), !is.na(prot_lfc))

cat("  Total matched gene-protein pairs:", nrow(all_matched), "\n")

# ----------------------------------------------------------------
# 3. CONCORDANCE STATS
# ----------------------------------------------------------------
cat("\nStep 3: Calculating concordance...\n")

all_matched <- all_matched %>%
  mutate(
    concordant   = (rna_lfc > 0) == (prot_lfc > 0),
    direction    = case_when(
      rna_lfc  > 0 & prot_lfc  > 0 ~ "Both UP",
      rna_lfc  < 0 & prot_lfc  < 0 ~ "Both DOWN",
      rna_lfc  > 0 & prot_lfc  < 0 ~ "RNA up / Protein down",
      rna_lfc  < 0 & prot_lfc  > 0 ~ "RNA down / Protein up"
    ),
    both_sig = !is.na(rna_padj) & !is.na(prot_padj) &
               rna_padj < 0.05  & prot_padj < 0.05
  )

cat("  Overall concordance:",
    round(mean(all_matched$concordant) * 100, 1), "%\n")
cat("  Direction breakdown:\n")
print(table(all_matched$direction))

# Subset: both significant
both_sig <- all_matched %>% filter(both_sig)
cat("\n  Both significant (padj < 0.05 in both):", nrow(both_sig), "\n")
cat("  Concordant among both-sig:",
    sum(both_sig$concordant), "/", nrow(both_sig),
    paste0("(", round(mean(both_sig$concordant) * 100, 1), "%)"), "\n")
cat("  Direction breakdown (both sig):\n")
print(table(both_sig$direction))

# ----------------------------------------------------------------
# 4. SCATTER PLOT — mRNA log2FC vs Protein log2FC
# ----------------------------------------------------------------
cat("\nStep 4: Generating scatter plot...\n")

# Colour points by significance status
all_matched <- all_matched %>%
  mutate(
    sig_status = case_when(
      both_sig & concordant   ~ "Both sig, concordant",
      both_sig & !concordant  ~ "Both sig, discordant",
      rna_padj  < 0.05 & !is.na(rna_padj)  ~ "RNA sig only",
      prot_padj < 0.05 & !is.na(prot_padj) ~ "Protein sig only",
      TRUE ~ "Not significant"
    )
  )

sig_colors <- c(
  "Both sig, concordant"  = "#D73027",
  "Both sig, discordant"  = "#FFA500",
  "RNA sig only"          = "#4575B4",
  "Protein sig only"      = "#91BFDB",
  "Not significant"       = "grey80"
)

# Label top 20 concordant both-sig pairs
label_df <- both_sig %>%
  filter(concordant) %>%
  arrange(rna_padj) %>%
  slice_head(n = 20)

p_scatter <- ggplot(all_matched,
                    aes(x = rna_lfc, y = prot_lfc, colour = sig_status)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_smooth(data = both_sig, aes(x = rna_lfc, y = prot_lfc),
              method = "lm", se = TRUE, colour = "black",
              linewidth = 0.8, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = rna_lfc, y = prot_lfc, label = gene),
    size = 2.8, colour = "black", max.overlaps = 20,
    box.padding = 0.3
  ) +
  scale_colour_manual(values = sig_colors, name = "Significance") +
  labs(
    title    = "Transcript-Protein Concordance",
    subtitle = "GSE137344 (ileum mRNA) vs PXD001608 (colon tissue proteomics)",
    x        = "mRNA log2 Fold Change (IBD vs nonIBD)",
    y        = "Protein log2 Fold Change (IBD vs nonIBD)",
    caption  = paste0("n = ", nrow(all_matched), " matched gene-protein pairs | ",
                      nrow(both_sig), " both significant (padj < 0.05)")
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"))

# Check if ggrepel is available; if not, skip labels
tryCatch({
  ggsave(file.path(OUT_DIR, "concordance_scatter.pdf"),
         p_scatter, width = 10, height = 8)
  cat("  Saved: concordance_scatter.pdf\n")
}, error = function(e) {
  # Fall back without labels
  p_scatter2 <- p_scatter +
    ggplot2::geom_text(
      data = label_df %>% slice_head(n = 10),
      aes(x = rna_lfc, y = prot_lfc, label = gene),
      size = 2.5, colour = "black", vjust = -0.5, inherit.aes = FALSE
    )
  ggsave(file.path(OUT_DIR, "concordance_scatter.pdf"),
         p_scatter2, width = 10, height = 8)
  cat("  Saved: concordance_scatter.pdf (without ggrepel labels)\n")
})

# ----------------------------------------------------------------
# 5. HEATMAP OF CONCORDANT BOTH-SIG PAIRS
# ----------------------------------------------------------------
cat("\nStep 5: Generating concordance heatmap...\n")

# Take top 30 concordant both-sig genes (by combined significance)
top_concordant <- both_sig %>%
  filter(concordant) %>%
  mutate(combined_score = -log10(rna_padj) + (-log10(prot_padj))) %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = 30)

cat("  Top concordant both-sig pairs:", nrow(top_concordant), "\n")
cat("  Direction breakdown:\n")
print(table(top_concordant$direction))

if (nrow(top_concordant) >= 5) {

  # Get normalised expression for these genes from GSE137344
  norm137  <- rds137$normalized_matrix
  cd137    <- rds137$col_data
  cd137$Diagnosis <- ifelse(cd137$diagnosis_binary == "IBD", "IBD", "nonIBD")
  rownames(cd137)  <- cd137$sample_id

  # Get protein matrix from PXD001608
  norm_pxd <- pxdrds$normalized_matrix
  cd_pxd   <- pxdrds$col_data
  cd_pxd$Diagnosis <- ifelse(cd_pxd$diagnosis_binary == "IBD", "IBD", "nonIBD")
  rownames(cd_pxd)  <- cd_pxd$sample_id

  genes_to_show <- top_concordant$gene
  avail_rna  <- genes_to_show[genes_to_show %in% rownames(norm137)]
  avail_prot <- genes_to_show[genes_to_show %in% rownames(norm_pxd)]

  # Z-score scale each matrix
  zscore <- function(m) {
    t(scale(t(m)))
  }

  zrna  <- zscore(norm137[avail_rna,  , drop = FALSE])
  zprot <- zscore(norm_pxd[avail_prot, , drop = FALSE])
  zrna[zrna   >  3] <-  3; zrna[zrna   < -3] <- -3
  zprot[zprot >  3] <-  3; zprot[zprot < -3] <- -3

  heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  anno_colors    <- list(Diagnosis = c(IBD = "#D73027", nonIBD = "#4575B4"))

  # mRNA heatmap
  anno_rna <- data.frame(Diagnosis = cd137$Diagnosis, row.names = rownames(cd137))
  anno_rna <- anno_rna[order(anno_rna$Diagnosis), , drop = FALSE]
  zrna_ord <- zrna[, rownames(anno_rna), drop = FALSE]

  pdf(file.path(OUT_DIR, "concordance_heatmap_mRNA.pdf"), width = 14, height = 9)
  pheatmap(zrna_ord,
           color            = heatmap_colors,
           annotation_col   = anno_rna,
           annotation_colors= anno_colors,
           cluster_rows     = TRUE,
           cluster_cols     = FALSE,
           show_colnames    = FALSE,
           show_rownames    = TRUE,
           fontsize_row     = 8,
           fontsize         = 9,
           main             = "mRNA Expression — Concordant Gene-Protein Pairs\nGSE137344 Ileum Transcriptomics",
           border_color     = NA)
  dev.off()
  cat("  Saved: concordance_heatmap_mRNA.pdf\n")

  # Protein heatmap
  anno_pxd <- data.frame(Diagnosis = cd_pxd$Diagnosis, row.names = rownames(cd_pxd))
  anno_pxd <- anno_pxd[order(anno_pxd$Diagnosis), , drop = FALSE]
  zprot_ord <- zprot[, rownames(anno_pxd), drop = FALSE]

  pdf(file.path(OUT_DIR, "concordance_heatmap_protein.pdf"), width = 10, height = 9)
  pheatmap(zprot_ord,
           color            = heatmap_colors,
           annotation_col   = anno_pxd,
           annotation_colors= anno_colors,
           cluster_rows     = TRUE,
           cluster_cols     = FALSE,
           show_colnames    = FALSE,
           show_rownames    = TRUE,
           fontsize_row     = 8,
           fontsize         = 9,
           main             = "Protein Abundance — Concordant Gene-Protein Pairs\nPXD001608 Colon Tissue Proteomics",
           border_color     = NA)
  dev.off()
  cat("  Saved: concordance_heatmap_protein.pdf\n")
}

# ----------------------------------------------------------------
# 6. SAVE CONCORDANCE TABLE
# ----------------------------------------------------------------
cat("\nStep 6: Saving concordance table...\n")

write.csv(
  all_matched %>% arrange(desc(both_sig), rna_padj),
  file.path(OUT_DIR, "concordance_table_all.csv"),
  row.names = FALSE
)
write.csv(
  both_sig %>% arrange(concordant, rna_padj) %>%
    select(gene, rna_lfc, prot_lfc, rna_padj, prot_padj,
           direction, concordant),
  file.path(OUT_DIR, "concordance_table_both_sig.csv"),
  row.names = FALSE
)
cat("  Saved: concordance_table_all.csv\n")
cat("  Saved: concordance_table_both_sig.csv\n")

# ----------------------------------------------------------------
# 7. SUMMARY
# ----------------------------------------------------------------
cat("\n================================================================\n")
cat("=== CONCORDANCE SUMMARY ===\n")
cat("================================================================\n")
cat("Total matched gene-protein pairs:", nrow(all_matched), "\n")
cat("Overall directional concordance:  ",
    round(mean(all_matched$concordant)*100,1), "%\n\n")
cat("Both significant (padj < 0.05 in both datasets):", nrow(both_sig), "\n")
cat("  Concordant:", sum(both_sig$concordant),
    paste0("(", round(mean(both_sig$concordant)*100,1), "%)"), "\n")
cat("  Discordant:", sum(!both_sig$concordant),
    paste0("(", round(mean(!both_sig$concordant)*100,1), "%)"), "\n\n")
cat("Direction breakdown (both sig):\n")
print(table(both_sig$direction))
cat("\nTop 10 concordant both-sig pairs:\n")
print(
  both_sig %>%
    filter(concordant) %>%
    arrange(rna_padj) %>%
    select(gene, rna_lfc, prot_lfc, rna_padj, prot_padj) %>%
    mutate(across(where(is.numeric), ~round(., 4))) %>%
    head(10)
)
