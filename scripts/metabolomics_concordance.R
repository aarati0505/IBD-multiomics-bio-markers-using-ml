# ================================================================
# METABOLOMICS CONCORDANCE ANALYSES
#
# Analysis 1: Cross-dataset concordance
#   MTBLS8713 (named metabolites) vs HMP2 BIOM (annotated features)
#   Question: Do the same metabolites appear significant in both
#   independent cohorts?
#
# Analysis 2: Cross-omics concordance (same HMP2 patients)
#   HMP2 metabolomics vs HMP2 metaproteomics ECS
#   Question: Do samples with depleted enzymes also have depleted
#   corresponding metabolites?
#
# Output: figures/concordance/
#   metabolomics_cross_dataset_concordance.pdf
#   metabolomics_metaproteomics_concordance.pdf
#   metabolomics_concordance_table.csv
# ================================================================

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

BASE    <- "C:/Users/Aarati/Desktop/IBD datasets"
OUT_DIR <- file.path(BASE, "figures", "concordance")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("Metabolomics Concordance Analyses\n")
cat("================================================================\n\n")

# ----------------------------------------------------------------
# LOAD ALL REQUIRED DATA
# ----------------------------------------------------------------
cat("Loading data...\n")

# MTBLS8713
mtbls_de  <- read.csv(file.path(BASE,
              "processed/MTBLS8713_metabolomics/MTBLS8713_DE_results.csv"),
              stringsAsFactors = FALSE)
mtbls_sig <- read.csv(file.path(BASE,
              "processed/MTBLS8713_metabolomics/MTBLS8713_sig_features.csv"),
              stringsAsFactors = FALSE)

# HMP2 metabolomics
hmp2_de   <- read.csv(file.path(BASE,
              "processed/HMP2_metabolomics/HMP2_metabolomics_DE_results.csv"),
              stringsAsFactors = FALSE)
hmp2_sig  <- read.csv(file.path(BASE,
              "processed/HMP2_metabolomics/HMP2_metabolomics_sig_features.csv"),
              stringsAsFactors = FALSE)

# HMP2 metaproteomics ECS
ecs_de    <- read.csv(file.path(BASE,
              "processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_DE_results.csv"),
              stringsAsFactors = FALSE)
ecs_sig   <- read.csv(file.path(BASE,
              "processed/HMP2_metaproteomics/HMP2_metaproteomics_ecs_sig_features.csv"),
              stringsAsFactors = FALSE)

# C18n annotation file (maps HMP2 feature codes to named metabolites)
c18n_anno <- read.csv(
  file.path(BASE, "metabolomics/c18n_annotation.csv"),
  stringsAsFactors = FALSE
)
cat("  C18n annotation entries:", nrow(c18n_anno), "\n")

cat("  MTBLS8713:", nrow(mtbls_de), "features,", nrow(mtbls_sig), "significant\n")
cat("  HMP2 metabolomics:", nrow(hmp2_de), "features,", nrow(hmp2_sig), "significant\n")
cat("  HMP2 ECS:", nrow(ecs_de), "features,", nrow(ecs_sig), "significant\n")

# ================================================================
# ANALYSIS 1: CROSS-DATASET METABOLITE CONCORDANCE
# MTBLS8713 vs annotated HMP2 BIOM
# ================================================================
cat("\n================================================================\n")
cat("Analysis 1: Cross-dataset concordance (MTBLS8713 vs HMP2 BIOM)\n")
cat("================================================================\n")

# Clean metabolite names for matching
clean_name <- function(x) {
  x <- tolower(x)
  x <- gsub("_", " ", x)
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- trimws(x)
  x
}

# Annotate HMP2 features that have a known name
# c18n_anno has columns: Method, Compound (e.g. QI1629), m/z, RT, Metabolite
# HMP2 feature IDs look like "C18n_QI1629" — strip prefix for matching
hmp2_de$compound_code <- sub("^C18n_", "", hmp2_de$feature)
hmp2_de$compound_code <- sub("^C18n", "", hmp2_de$compound_code)

c18n_anno$compound_clean <- paste0("QI", gsub("QI", "", c18n_anno$Compound))
c18n_anno$metabolite_clean <- clean_name(c18n_anno[,5])  # last column = metabolite name

# Merge annotation onto HMP2 DE results
hmp2_anno <- hmp2_de %>%
  left_join(
    c18n_anno %>% select(Compound, metabolite_name = 5),
    by = c("compound_code" = "Compound")
  )

annotated_count <- sum(!is.na(hmp2_anno$metabolite_name))
cat("  HMP2 features with annotation:", annotated_count, "out of", nrow(hmp2_de), "\n")

# Now match MTBLS8713 metabolites to annotated HMP2 features
mtbls_de$metabolite_clean <- clean_name(mtbls_de$metabolite)
hmp2_anno$metabolite_clean <- ifelse(
  is.na(hmp2_anno$metabolite_name), NA,
  clean_name(hmp2_anno$metabolite_name)
)

# Join on cleaned metabolite name
cross_match <- inner_join(
  mtbls_de %>% select(metabolite, metabolite_clean,
                      mtbls_lfc = log2FC, mtbls_padj = padj),
  hmp2_anno %>% filter(!is.na(metabolite_clean)) %>%
    select(feature, metabolite_clean,
           hmp2_lfc = log2FC, hmp2_padj = padj),
  by = "metabolite_clean"
)

cat("  Matched metabolites between datasets:", nrow(cross_match), "\n")

if (nrow(cross_match) > 0) {
  # Concordance
  cross_match <- cross_match %>%
    filter(!is.na(mtbls_lfc), !is.na(hmp2_lfc)) %>%
    mutate(
      concordant = (mtbls_lfc > 0) == (hmp2_lfc > 0),
      both_sig   = !is.na(mtbls_padj) & !is.na(hmp2_padj) &
                   mtbls_padj < 0.05 & hmp2_padj < 0.05,
      direction  = case_when(
        mtbls_lfc > 0 & hmp2_lfc > 0 ~ "Both UP in IBD",
        mtbls_lfc < 0 & hmp2_lfc < 0 ~ "Both DOWN in IBD",
        mtbls_lfc > 0 & hmp2_lfc < 0 ~ "MTBLS up / HMP2 down",
        mtbls_lfc < 0 & hmp2_lfc > 0 ~ "MTBLS down / HMP2 up"
      )
    )

  cat("  Overall concordance:",
      round(mean(cross_match$concordant)*100, 1), "%\n")
  both_sig_cross <- cross_match %>% filter(both_sig)
  cat("  Both significant:", nrow(both_sig_cross), "\n")
  if (nrow(both_sig_cross) > 0) {
    cat("  Both sig concordant:", sum(both_sig_cross$concordant),
        paste0("(", round(mean(both_sig_cross$concordant)*100,1),"%)"), "\n")
    cat("  Direction breakdown:\n")
    print(table(both_sig_cross$direction))
    cat("\n  Top concordant both-sig metabolites:\n")
    print(both_sig_cross %>% filter(concordant) %>%
            arrange(mtbls_padj) %>%
            select(metabolite, mtbls_lfc, hmp2_lfc, mtbls_padj, hmp2_padj, direction) %>%
            mutate(across(where(is.numeric), ~round(., 3))) %>%
            head(10))
  }

  # Scatter plot
  p1 <- ggplot(cross_match,
               aes(x = mtbls_lfc, y = hmp2_lfc,
                   colour = both_sig & concordant)) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype="dashed", colour="grey50") +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50") +
    scale_colour_manual(
      values = c("TRUE" = "#D73027", "FALSE" = "grey70"),
      labels = c("TRUE" = "Both sig & concordant", "FALSE" = "Other"),
      name   = ""
    ) +
    labs(
      title    = "Cross-Dataset Metabolite Concordance",
      subtitle = "MTBLS8713 vs HMP2 BIOM (annotated features)",
      x        = "MTBLS8713 log2FC (IBD vs nonIBD)",
      y        = "HMP2 BIOM log2FC (IBD vs nonIBD)",
      caption  = paste0("n = ", nrow(cross_match),
                        " matched metabolites | ",
                        nrow(both_sig_cross), " both significant")
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # Add labels for both-sig concordant
  if (nrow(both_sig_cross) > 0 && nrow(both_sig_cross %>% filter(concordant)) > 0) {
    label_df <- both_sig_cross %>% filter(concordant) %>%
      arrange(mtbls_padj) %>% slice_head(n = 15)
    tryCatch({
      library(ggrepel)
      p1 <- p1 + ggrepel::geom_text_repel(
        data = label_df,
        aes(x = mtbls_lfc, y = hmp2_lfc, label = metabolite),
        size = 3, colour = "black", max.overlaps = 15,
        inherit.aes = FALSE
      )
    }, error = function(e) NULL)
  }

  ggsave(file.path(OUT_DIR, "metabolomics_cross_dataset_concordance.pdf"),
         p1, width = 9, height = 7)
  cat("  Saved: metabolomics_cross_dataset_concordance.pdf\n")

  # Save table
  write.csv(cross_match %>% arrange(desc(both_sig), mtbls_padj),
            file.path(OUT_DIR, "metabolomics_concordance_cross_dataset.csv"),
            row.names = FALSE)
  cat("  Saved: metabolomics_concordance_cross_dataset.csv\n")
} else {
  cat("  No matched metabolites found — annotation may need adjustment\n")
}

# ================================================================
# ANALYSIS 2: HMP2 METABOLOMICS vs HMP2 METAPROTEOMICS (ECS)
# Same patients, different measurement types
# Match via KEGG/enzyme-metabolite relationships
# ================================================================
cat("\n================================================================\n")
cat("Analysis 2: HMP2 metabolomics vs HMP2 metaproteomics (ECS)\n")
cat("Same patient cohort — fold change correlation\n")
cat("================================================================\n")

# Known enzyme-metabolite pairs from biochemistry
# Format: enzyme EC number pattern -> related metabolites in MTBLS8713
enzyme_metabolite_pairs <- data.frame(
  ec_pattern   = c("1.1.1.27", "4.1.1.32", "2.7.9.1", "1.4.1.13",
                   "1.4.1.14", "6.3.1.2",  "1.2.1.12", "1.1.1.3",
                   "4.2.1.3",  "1.1.1.37", "4.1.3.7"),
  metabolite   = c("Pyruvic_acid", "Pyruvic_acid", "Pyruvic_acid",
                   "L_Glutamine",  "L_Glutamine",  "L_Glutamine",
                   "Pyruvic_acid", "L_Threonine",
                   "Citric_acid",  "Malic_acid",   "Citric_acid"),
  enzyme_name  = c("Lactate dehydrogenase", "PEP carboxykinase",
                   "Pyruvate dikinase",    "Glutamate synthase (NADPH)",
                   "Glutamate synthase (NADH)", "Glutamine synthase",
                   "GAPDH",               "Threonine dehydrogenase",
                   "Aconitase",           "Malate dehydrogenase",
                   "Citrate synthase"),
  stringsAsFactors = FALSE
)

cat("  Testing", nrow(enzyme_metabolite_pairs),
    "known enzyme-metabolite relationships\n")

# Get ECS fold changes for these enzymes
# ECS feature IDs look like "1.1.1.27: Lactate dehydrogenase"
get_ecs_lfc <- function(ec_pattern) {
  rows <- ecs_de[grepl(ec_pattern, ecs_de$feature, fixed = TRUE), ]
  if (nrow(rows) == 0) return(NA)
  rows$log2FC[1]
}

get_mtbls_lfc <- function(met_name) {
  rows <- mtbls_de[mtbls_de$metabolite == met_name, ]
  if (nrow(rows) == 0) return(NA)
  rows$log2FC[1]
}

get_mtbls_padj <- function(met_name) {
  rows <- mtbls_de[mtbls_de$metabolite == met_name, ]
  if (nrow(rows) == 0) return(NA)
  rows$padj[1]
}

get_ecs_padj <- function(ec_pattern) {
  rows <- ecs_de[grepl(ec_pattern, ecs_de$feature, fixed = TRUE), ]
  if (nrow(rows) == 0) return(NA)
  rows$padj[1]
}

enzyme_metabolite_pairs$ecs_lfc    <- sapply(enzyme_metabolite_pairs$ec_pattern, get_ecs_lfc)
enzyme_metabolite_pairs$ecs_padj   <- sapply(enzyme_metabolite_pairs$ec_pattern, get_ecs_padj)
enzyme_metabolite_pairs$mtbls_lfc  <- sapply(enzyme_metabolite_pairs$metabolite, get_mtbls_lfc)
enzyme_metabolite_pairs$mtbls_padj <- sapply(enzyme_metabolite_pairs$metabolite, get_mtbls_padj)

pairs_complete <- enzyme_metabolite_pairs %>%
  filter(!is.na(ecs_lfc) & !is.na(mtbls_lfc))

cat("  Enzyme-metabolite pairs with data:", nrow(pairs_complete), "\n")

if (nrow(pairs_complete) > 0) {
  pairs_complete <- pairs_complete %>%
    mutate(
      concordant = (ecs_lfc > 0) == (mtbls_lfc > 0),
      both_sig   = !is.na(ecs_padj) & !is.na(mtbls_padj) &
                   ecs_padj < 0.05 & mtbls_padj < 0.05
    )

  cat("  Concordant pairs:", sum(pairs_complete$concordant),
      "out of", nrow(pairs_complete), "\n")
  cat("  Both significant:", sum(pairs_complete$both_sig), "\n")

  cat("\n  Full results:\n")
  print(pairs_complete %>%
          select(enzyme_name, metabolite, ecs_lfc, mtbls_lfc,
                 ecs_padj, mtbls_padj, concordant, both_sig) %>%
          mutate(across(where(is.numeric), ~round(.,3))))

  # Scatter plot
  p2 <- ggplot(pairs_complete,
               aes(x = ecs_lfc, y = mtbls_lfc,
                   colour = both_sig & concordant,
                   label  = paste0(enzyme_name, "\n(", metabolite, ")"))) +
    geom_point(size = 3.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", colour = "grey70") +
    scale_colour_manual(
      values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
      labels = c("TRUE" = "Both sig & concordant",
                 "FALSE" = "Not both sig"),
      name = ""
    ) +
    tryCatch({
      ggrepel::geom_text_repel(size = 2.8, colour = "black",
                               max.overlaps = 20, box.padding = 0.4)
    }, error = function(e) geom_text(size = 2.5, vjust = -0.7, colour = "black")) +
    labs(
      title    = "Enzyme–Metabolite Concordance",
      subtitle = "HMP2 Metaproteomics ECS vs MTBLS8713 Metabolomics",
      x        = "Enzyme log2FC in IBD (stool metaproteomics)",
      y        = "Metabolite log2FC in IBD (MTBLS8713)",
      caption  = "Dotted diagonal = perfect concordance line"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(OUT_DIR, "metabolomics_metaproteomics_concordance.pdf"),
         p2, width = 10, height = 8)
  cat("  Saved: metabolomics_metaproteomics_concordance.pdf\n")

  write.csv(pairs_complete,
            file.path(OUT_DIR, "metabolomics_enzyme_concordance_table.csv"),
            row.names = FALSE)
  cat("  Saved: metabolomics_enzyme_concordance_table.csv\n")
}

# ================================================================
# SUMMARY
# ================================================================
cat("\n================================================================\n")
cat("=== METABOLOMICS CONCORDANCE SUMMARY ===\n")
cat("================================================================\n\n")
cat("Analysis 1 — Cross-dataset (MTBLS8713 vs HMP2 BIOM):\n")
if (exists("cross_match") && nrow(cross_match) > 0) {
  cat("  Matched metabolites:", nrow(cross_match), "\n")
  cat("  Both sig & concordant:", sum(cross_match$both_sig & cross_match$concordant), "\n")
}
cat("\nAnalysis 2 — Enzyme-metabolite (ECS vs MTBLS8713):\n")
if (exists("pairs_complete") && nrow(pairs_complete) > 0) {
  cat("  Pairs tested:", nrow(pairs_complete), "\n")
  cat("  Concordant:", sum(pairs_complete$concordant), "\n")
  cat("  Both significant:", sum(pairs_complete$both_sig), "\n")
}
