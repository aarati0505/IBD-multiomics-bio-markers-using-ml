# ================================================================
# Check known IBD biomarkers across all significant results
# ================================================================
library(dplyr)

BASE <- "C:/Users/Aarati/Desktop/IBD datasets"

# Load all sig results
pxd_sig   <- read.csv(file.path(BASE,"processed/proteomics_MaxQuant/proteomics_MaxQuant_sig_proteins.csv"), stringsAsFactors=FALSE)
pxd_de    <- read.csv(file.path(BASE,"processed/proteomics_MaxQuant/proteomics_MaxQuant_DE_results.csv"),   stringsAsFactors=FALSE)
deg137    <- read.csv(file.path(BASE,"processed/GSE137344/GSE137344_significant_DEGs.csv"),               stringsAsFactors=FALSE)
deg137_all<- read.csv(file.path(BASE,"processed/GSE137344/GSE137344_DEG_results.csv"),                    stringsAsFactors=FALSE)
deg186    <- read.csv(file.path(BASE,"processed/GSE186507/GSE186507_significant_DEGs.csv"),               stringsAsFactors=FALSE)
deg186_all<- read.csv(file.path(BASE,"processed/GSE186507/GSE186507_DEG_results.csv"),                    stringsAsFactors=FALSE)
ecs_sig   <- read.csv(file.path(BASE,"processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_sig_features.csv"), stringsAsFactors=FALSE)
kos_sig   <- read.csv(file.path(BASE,"processed/HMP2_metaproteomics/HMP2_metaproteomics_kos_sig_features.csv"), stringsAsFactors=FALSE)
mtbls_sig <- read.csv(file.path(BASE,"processed/MTBLS8713_metabolomics/MTBLS8713_sig_features.csv"),      stringsAsFactors=FALSE)
conc_294  <- read.csv(file.path(BASE,"figures/concordance/concordant_294_annotated.csv"),                 stringsAsFactors=FALSE)
conc_all  <- read.csv(file.path(BASE,"figures/concordance/concordance_table_all.csv"),                    stringsAsFactors=FALSE)

cat("Loaded all files.\n")
cat("PXD sig proteins:", nrow(pxd_sig), "\n")
cat("GSE137344 sig DEGs:", nrow(deg137), "\n")
cat("GSE186507 sig DEGs:", nrow(deg186), "\n")
cat("Concordant 294:", nrow(conc_294), "\n")

# Known IBD biomarkers
biomarkers <- list(
  # VALIDATED CLINICAL TESTS
  LTF      = "Lactoferrin — validated stool IBD test, neutrophil protein",
  S100A8   = "Calprotectin subunit — FDA-approved stool IBD test",
  S100A9   = "Calprotectin subunit — FDA-approved stool IBD test",
  S100A12  = "Calgranulin-C — stool inflammatory marker",
  LCN2     = "Lipocalin-2 — validated stool IBD test (neutrophil)",
  OLFM4    = "Olfactomedin-4 — neutrophil/epithelial IBD marker",
  CHI3L1   = "YKL-40 — serum IBD disease activity marker",
  SERPINA1 = "Alpha-1 antitrypsin — serum IBD marker",
  # INFLAMMATORY CYTOKINES
  CXCL8    = "IL-8 — key neutrophil chemokine in IBD",
  TNF      = "TNF-alpha — major IBD therapeutic target",
  IL1B     = "IL-1 beta — inflammasome in IBD",
  IL6      = "Interleukin-6 — systemic inflammation",
  IL17A    = "IL-17 — Th17 pathway in IBD",
  IFNG     = "Interferon-gamma — Th1 elevated in CD",
  IL10     = "IL-10 — anti-inflammatory, deficient in IBD",
  IL23A    = "IL-23 — IBD genetic risk and drug target",
  # JAK-STAT / IBD DRUG TARGETS
  JAK1     = "JAK1 — IBD target (tofacitinib/upadacitinib)",
  JAK2     = "JAK2 — JAK-STAT pathway",
  STAT3    = "STAT3 — JAK-STAT signalling in IBD",
  # IBD GENETIC RISK GENES
  NOD2     = "NOD2 — most important CD genetic risk gene",
  ATG16L1  = "ATG16L1 — autophagy, IBD GWAS top hit",
  ATG16L2  = "ATG16L2 — autophagy, IBD associated",
  LRRK2    = "LRRK2 — IBD and Parkinson risk gene",
  PTPN2    = "PTPN2 — IBD genetic risk, T-cell regulation",
  CARD9    = "CARD9 — innate immunity, IBD risk",
  IRF4     = "IRF4 — immune regulation in IBD",
  # EPITHELIAL BARRIER
  MUC2     = "Mucin-2 — major intestinal mucus, depleted in UC",
  MUC1     = "Mucin-1 — epithelial barrier",
  TFF3     = "Trefoil factor 3 — mucosal healing marker",
  CDH1     = "E-cadherin — epithelial junction",
  CLDN1    = "Claudin-1 — tight junction barrier",
  VIL1     = "Villin-1 — intestinal epithelial marker",
  WFDC2    = "WFDC2 — epithelial protein, depleted in IBD",
  FCGBP    = "IgG Fc binding protein — depleted in UC",
  # MATRIX REMODELLING
  MMP3     = "MMP-3 — matrix metalloproteinase in IBD",
  MMP9     = "MMP-9 — neutrophil gelatinase in IBD",
  ELN      = "Elastin — ECM protein lost in UC damage",
  COL1A1   = "Collagen type I — ECM remodelling",
  COL1A2   = "Collagen type I alpha-2 — ECM",
  FN1      = "Fibronectin — ECM and tissue repair",
  # METABOLISM / OXIDATIVE STRESS
  HMGCS2   = "HMG-CoA synthase 2 — energy metabolism, depleted in IBD",
  PPARG    = "PPAR-gamma — lipid metabolism, anti-inflammatory in colon",
  HIF1A    = "HIF-1alpha — hypoxia response in IBD",
  ACOX1    = "Acyl-CoA oxidase — fatty acid oxidation",
  # ANTIMICROBIAL / MICROBIOME INTERFACE
  DEFB1    = "Defensin beta-1 — antimicrobial peptide, reduced in IBD",
  DEFA5    = "Defensin alpha-5 — Paneth cell, dysregulated in CD",
  REG3A    = "REG3A — antimicrobial lectin, IBD marker",
  # KNOWN IBD lncRNAs / ncRNAs
  NEAT1    = "NEAT1 — lncRNA, nuclear paraspeckle, IBD inflammation",
  # OTHER WELL-KNOWN
  VEGFA    = "VEGF-A — angiogenesis in IBD",
  MPO      = "Myeloperoxidase — neutrophil activation",
  HP       = "Haptoglobin — acute phase protein elevated in IBD",
  ORM1     = "Orosomucoid/Alpha-1-acid glycoprotein — IBD acute phase",
  CRP      = "C-reactive protein — classic IBD inflammation marker",
  APOA1    = "Apolipoprotein A1 — reduced in IBD inflammation",
  ALB      = "Albumin — reduced in severe IBD"
)

cat("\n================================================================\n")
cat("KNOWN IBD BIOMARKERS FOUND IN SIGNIFICANT RESULTS\n")
cat("================================================================\n\n")

results <- data.frame(
  gene = character(), description = character(),
  in_pxd_sig = logical(), pxd_lfc = numeric(), pxd_padj = numeric(),
  in_gse137_sig = logical(), rna137_lfc = numeric(), rna137_padj = numeric(),
  in_gse186_sig = logical(), rna186_lfc = numeric(), rna186_padj = numeric(),
  in_pxd_any = logical(), in_rna137_any = logical(),
  in_concordance_294 = logical(), concordance_direction = character(),
  stringsAsFactors = FALSE
)

for (gene in names(biomarkers)) {
  desc <- biomarkers[[gene]]

  # Check proteomics sig
  pxd_row  <- pxd_sig[pxd_sig$protein == gene, ]
  pxd_any  <- pxd_de[pxd_de$protein == gene, ]
  in_pxd   <- nrow(pxd_row) > 0

  # Check RNA 137
  r137_row <- deg137[deg137$gene == gene, ]
  r137_any <- deg137_all[deg137_all$gene == gene, ]
  in_r137  <- nrow(r137_row) > 0

  # Check RNA 186
  r186_row <- deg186[deg186$gene == gene, ]
  in_r186  <- nrow(r186_row) > 0

  # Check concordance
  conc_row <- conc_294[conc_294$gene_symbol == gene, ]
  in_conc  <- nrow(conc_row) > 0
  conc_dir <- if (in_conc) conc_row$direction[1] else ""

  # Only report if found in at least one sig list
  if (in_pxd || in_r137 || in_r186 || in_conc) {
    row <- data.frame(
      gene = gene, description = desc,
      in_pxd_sig = in_pxd,
      pxd_lfc    = if (in_pxd) round(as.numeric(pxd_row$log2FC[1]),2) else
                   if (nrow(pxd_any)>0) round(as.numeric(pxd_any$log2FC[1]),2) else NA,
      pxd_padj   = if (in_pxd) signif(as.numeric(pxd_row$padj[1]),3) else
                   if (nrow(pxd_any)>0) signif(as.numeric(pxd_any$padj[1]),3) else NA,
      in_gse137_sig = in_r137,
      rna137_lfc = if (in_r137) round(as.numeric(r137_row$log2FoldChange[1]),2) else
                   if (nrow(r137_any)>0) round(as.numeric(r137_any$log2FoldChange[1]),2) else NA,
      rna137_padj= if (in_r137) signif(as.numeric(r137_row$padj[1]),3) else
                   if (nrow(r137_any)>0) signif(as.numeric(r137_any$padj[1]),3) else NA,
      in_gse186_sig = in_r186,
      rna186_lfc = if (in_r186) round(as.numeric(r186_row$logFC[1]),2) else NA,
      rna186_padj= if (in_r186) signif(as.numeric(r186_row$adj.P.Val[1]),3) else NA,
      in_pxd_any = nrow(pxd_any) > 0,
      in_rna137_any = nrow(r137_any) > 0,
      in_concordance_294 = in_conc,
      concordance_direction = conc_dir,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, row)
  }
}

# Print results grouped by category
cat("VALIDATED CLINICAL MARKERS:\n")
clinical <- c("LTF","S100A8","S100A9","S100A12","LCN2","OLFM4","CHI3L1","SERPINA1")
for (g in clinical) {
  r <- results[results$gene == g, ]
  if (nrow(r) > 0) {
    cat(sprintf("\n  %-10s | %s\n", g, r$description))
    if (r$in_pxd_sig)     cat(sprintf("    PROTEOMICS (sig)  : LFC=%+.2f  padj=%.2e\n", r$pxd_lfc, r$pxd_padj))
    else if (r$in_pxd_any) cat(sprintf("    PROTEOMICS (ns)   : LFC=%+.2f  padj=%.2e\n", r$pxd_lfc, r$pxd_padj))
    if (r$in_gse137_sig)  cat(sprintf("    RNA ileum (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna137_lfc, r$rna137_padj))
    else if (r$in_rna137_any) cat(sprintf("    RNA ileum (ns)    : LFC=%+.2f  padj=%.2e\n", r$rna137_lfc, r$rna137_padj))
    if (r$in_gse186_sig)  cat(sprintf("    RNA blood (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna186_lfc, r$rna186_padj))
    if (r$in_concordance_294) cat(sprintf("    *** IN CONCORDANCE 294 *** direction: %s\n", r$concordance_direction))
  }
}

cat("\n\nGENETIC RISK GENES:\n")
genetic <- c("NOD2","ATG16L1","ATG16L2","LRRK2","PTPN2","CARD9","IRF4")
for (g in genetic) {
  r <- results[results$gene == g, ]
  if (nrow(r) > 0) {
    cat(sprintf("\n  %-10s | %s\n", g, r$description))
    if (r$in_pxd_sig)      cat(sprintf("    PROTEOMICS (sig)  : LFC=%+.2f  padj=%.2e\n", r$pxd_lfc, r$pxd_padj))
    if (r$in_gse137_sig)   cat(sprintf("    RNA ileum (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna137_lfc, r$rna137_padj))
    if (r$in_gse186_sig)   cat(sprintf("    RNA blood (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna186_lfc, r$rna186_padj))
    if (r$in_concordance_294) cat(sprintf("    *** IN CONCORDANCE 294 *** direction: %s\n", r$concordance_direction))
  }
}

cat("\n\nEPITHELIAL BARRIER GENES:\n")
barrier <- c("MUC2","MUC1","TFF3","CDH1","CLDN1","VIL1","WFDC2","FCGBP")
for (g in barrier) {
  r <- results[results$gene == g, ]
  if (nrow(r) > 0) {
    cat(sprintf("\n  %-10s | %s\n", g, r$description))
    if (r$in_pxd_sig)      cat(sprintf("    PROTEOMICS (sig)  : LFC=%+.2f  padj=%.2e\n", r$pxd_lfc, r$pxd_padj))
    if (r$in_gse137_sig)   cat(sprintf("    RNA ileum (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna137_lfc, r$rna137_padj))
    if (r$in_gse186_sig)   cat(sprintf("    RNA blood (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna186_lfc, r$rna186_padj))
    if (r$in_concordance_294) cat(sprintf("    *** IN CONCORDANCE 294 *** direction: %s\n", r$concordance_direction))
  }
}

cat("\n\nINFLAMMATORY CYTOKINES & DRUG TARGETS:\n")
cytokines <- c("CXCL8","TNF","IL1B","IL6","IL17A","IFNG","IL10","JAK1","JAK2","STAT3","IL23A")
for (g in cytokines) {
  r <- results[results$gene == g, ]
  if (nrow(r) > 0) {
    cat(sprintf("\n  %-10s | %s\n", g, r$description))
    if (r$in_pxd_sig)      cat(sprintf("    PROTEOMICS (sig)  : LFC=%+.2f  padj=%.2e\n", r$pxd_lfc, r$pxd_padj))
    if (r$in_gse137_sig)   cat(sprintf("    RNA ileum (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna137_lfc, r$rna137_padj))
    if (r$in_gse186_sig)   cat(sprintf("    RNA blood (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna186_lfc, r$rna186_padj))
    if (r$in_concordance_294) cat(sprintf("    *** IN CONCORDANCE 294 *** direction: %s\n", r$concordance_direction))
  }
}

cat("\n\nMETABOLISM / ECM / OTHER:\n")
other <- c("HMGCS2","ELN","COL1A1","COL1A2","MMP3","MMP9","NEAT1","VEGFA","HP","APOA1","PPARG","WFDC2")
for (g in other) {
  r <- results[results$gene == g, ]
  if (nrow(r) > 0) {
    cat(sprintf("\n  %-10s | %s\n", g, r$description))
    if (r$in_pxd_sig)      cat(sprintf("    PROTEOMICS (sig)  : LFC=%+.2f  padj=%.2e\n", r$pxd_lfc, r$pxd_padj))
    if (r$in_gse137_sig)   cat(sprintf("    RNA ileum (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna137_lfc, r$rna137_padj))
    if (r$in_gse186_sig)   cat(sprintf("    RNA blood (sig)   : LFC=%+.2f  padj=%.2e\n", r$rna186_lfc, r$rna186_padj))
    if (r$in_concordance_294) cat(sprintf("    *** IN CONCORDANCE 294 *** direction: %s\n", r$concordance_direction))
  }
}

# Check metabolomics separately
cat("\n\n================================================================\n")
cat("KNOWN IBD METABOLITES IN MTBLS8713\n")
cat("================================================================\n")
known_metabolites <- c("Deoxycholic_acid","Lithocholic_acid","Isolithocholic_acid",
                       "Taurocholic_acid","Chenodeoxycholic_acid","Cholic_acid",
                       "Pyruvic_acid","Butyric_acid","L_Tryptophan","L_Glutamine",
                       "Kynurenic_acid","Indole_3_propionic_acid",
                       "Short_chain_fatty_acids","Acetic_acid","Propionic_acid",
                       "L_Carnitine","Trimethylamine_N_oxide","Phenylacetylglutamine")
mtbls_de <- read.csv(file.path(BASE,"processed/MTBLS8713_metabolomics/MTBLS8713_DE_results.csv"), stringsAsFactors=FALSE)
for (met in known_metabolites) {
  r <- mtbls_de[mtbls_de$metabolite == met, ]
  if (nrow(r) > 0) {
    sig_flag <- if (!is.na(r$padj) && r$padj < 0.05) "SIGNIFICANT" else "not sig"
    cat(sprintf("  %-40s LFC=%+.2f  padj=%.2e  [%s]\n",
                met, r$log2FC, r$padj, sig_flag))
  }
}

cat("\n=== SUMMARY ===\n")
cat("Total known biomarkers found in any sig list:", nrow(results), "\n")
cat("Found in proteomics (sig):", sum(results$in_pxd_sig), "\n")
cat("Found in RNA ileum (sig):", sum(results$in_gse137_sig), "\n")
cat("Found in RNA blood (sig):", sum(results$in_gse186_sig), "\n")
cat("Found in concordance 294:", sum(results$in_concordance_294), "\n")
