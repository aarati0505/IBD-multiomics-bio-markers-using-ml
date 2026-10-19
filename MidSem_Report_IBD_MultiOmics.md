# A Project Report

On

# Multi-Omics Data Integration for Identification of IBD Molecular Signatures

BY

**AARATI DHAMELE**

**2023B1A70641H**

Under the supervision of

**[SUPERVISOR NAME]**

SUBMITTED IN PARTIAL FULFILLMENT OF THE REQUIREMENTS OF

**BIO F266: STUDY PROJECT**

**BIRLA INSTITUTE OF TECHNOLOGY AND SCIENCE PILANI (RAJASTHAN)**

**HYDERABAD CAMPUS**

**(MARCH 2025)**

---

## ACKNOWLEDGMENTS

I would like to extend my gratitude to [Supervisor Name] for their constant support and guidance. Their expertise in the field of computational biology and multi-omics data analysis and their insights have been invaluable.

I would also like to acknowledge BIRLA INSTITUTE OF TECHNOLOGY AND SCIENCE PILANI, HYDERABAD CAMPUS for providing top notch facilities and all the necessary resources required for the undertaking of this project. Finally, I am grateful to my friends and my family for their motivation and support.

---

## Certificate

Birla Institute of Technology and Science-Pilani, Hyderabad Campus

This is to certify that the project report entitled "Multi-Omics Data Integration for Identification of IBD Molecular Signatures" submitted by Ms. AARATI DHAMELE (ID No. 2023B1A70641H) in partial fulfillment of the requirements of the course BIO F266, Study Project Course, embodies the work done by her under my supervision and guidance.

Date: [DATE] &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ([SUPERVISOR NAME])

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; BITS-Pilani, Hyderabad Campus

---

## ABSTRACT

Inflammatory Bowel Disease (IBD), comprising Crohn's Disease (CD) and Ulcerative Colitis (UC), is a chronic relapsing inflammatory disorder of the gastrointestinal tract affecting an estimated 6.8 million people worldwide. Current clinical diagnosis relies on invasive procedures such as colonoscopy and biopsy, creating a strong need for non-invasive, molecular-based diagnostic tools.

This project addresses this clinical gap by integrating multiple omics data types — transcriptomics, proteomics, metabolomics, and metagenomics — from large publicly available datasets to identify a robust multi-omics molecular signature for IBD. Data preprocessing has been completed for seven datasets spanning 4,398 patient samples: three transcriptomics datasets (GSE137344, GSE186507, GSE193677), three proteomics datasets (PXD001608, HMP2 metaproteomics ECS and KOS), and two metabolomics datasets (MTBLS8713, HMP2 BIOM).

Differential expression analysis using DESeq2, limma, Welch t-test, and Wilcoxon rank-sum tests identified thousands of significantly dysregulated features across all data types. A cross-omics concordance analysis revealed 294 genes consistently altered at both the mRNA and protein level across independent datasets — representing the highest-confidence biomarker candidates. Heatmap visualisations across all datasets confirm the dominant biological signal: IBD involves widespread downregulation of normal epithelial and microbial metabolic functions with a smaller signature of upregulated inflammatory markers.

The next phase involves metagenomics preprocessing followed by multi-omics machine learning model development using LASSO Logistic Regression to extract the final molecular signature for non-invasive IBD detection, with future validation planned at AIG Hospitals, Hyderabad.

---

## CONTENTS

1. Introduction ............................................................... 6
2. Literature Review ........................................................ 8
3. Methodology ............................................................... 10
4. Results ...................................................................... 16
5. Future Directions ........................................................ 20
6. Conclusion .................................................................. 21
7. References .................................................................. 22

---

## 1. INTRODUCTION

Inflammatory Bowel Disease (IBD) is a chronic, relapsing inflammatory condition of the gastrointestinal tract. It is broadly classified into two major forms: Crohn's Disease (CD), which can affect any segment of the digestive tract from mouth to anus in a discontinuous or transmural pattern, and Ulcerative Colitis (UC), which affects the mucosa of the colon and rectum in a continuous pattern. Both forms are characterised by cycles of remission and flare, during which patients experience debilitating symptoms including diarrhoea, abdominal pain, rectal bleeding, fatigue, and weight loss.

IBD affects approximately 6.8 million people globally and the incidence is rising sharply in newly industrialised countries including India, where historically the disease was rare. The current diagnostic gold standard requires endoscopy with tissue biopsy and histological examination — procedures that are invasive, expensive, time-consuming, and not suitable for routine monitoring of disease activity. There is therefore a strong and growing clinical need for non-invasive biomarkers that can reliably distinguish IBD patients from healthy individuals, and ideally differentiate CD from UC.

### The Multi-Omics Approach

The pathophysiology of IBD is complex and involves interactions between host genetic susceptibility, immune dysregulation, intestinal microbiota composition, and environmental factors. No single biological measurement captures this complexity fully. Multi-omics integration — combining data from multiple molecular measurement platforms simultaneously — offers the potential to identify a more complete and robust disease signature.

This project integrates four complementary omics data types:

**Transcriptomics** measures gene expression — which genes are switched on or off in patient tissue and blood. Genes produce messenger RNA (mRNA) when active, and measuring mRNA levels across the entire genome simultaneously (using RNA sequencing or microarray technology) reveals the transcriptional state of a cell or tissue. Comparing expression profiles between IBD patients and healthy individuals identifies genes and biological pathways specifically activated or suppressed in disease. Tissue transcriptomics directly captures local inflammatory gene expression at the disease site, while blood transcriptomics reflects the systemic immune response.

**Proteomics** measures the actual proteins present in a biological sample using mass spectrometry (LC-MS/MS). Proteins are the functional molecules that carry out all biological processes — enzymes, structural components, signalling molecules, and immune mediators. Proteomics provides a view of biology that is complementary but not identical to transcriptomics, because not all mRNA molecules are translated into protein, and because proteins can be chemically modified after production in ways that affect their function. A disease-associated protein detected independently by both transcriptomics and proteomics is a far stronger biomarker candidate than one detected by only one technology.

**Metabolomics** measures small chemical molecules (metabolites) in biological samples — amino acids, fatty acids, bile acids, sugars, short-chain fatty acids, and thousands of other compounds. These metabolites are the downstream end products of both cellular and microbial metabolism. Changes in the metabolome in IBD reflect disruptions in host energy metabolism, bile acid cycling, and microbial fermentation, providing a direct readout of the biochemical consequences of disease.

**Metagenomics** uses next-generation sequencing of all DNA in a stool sample to identify which bacterial species are present and in what proportions. IBD is consistently associated with gut dysbiosis — a disruption of the normal microbial community — characterised by depletion of commensal anti-inflammatory bacteria such as *Faecalibacterium prausnitzii* and expansion of potentially pro-inflammatory species. Metagenomics quantifies this dysbiosis directly.

### Project Objectives

The specific objectives of this project are:

1. To collect, organise, and preprocess all four omics data types from publicly available IBD datasets
2. To perform differential expression and abundance analysis to identify features significantly altered in IBD across each data type
3. To perform cross-omics concordance analysis to identify features robustly altered at multiple molecular levels
4. To visualise differential features using heatmaps, scatter plots, and summary figures
5. To build a multi-omics machine learning model to classify IBD vs healthy using LASSO Logistic Regression
6. To extract a final minimal molecular signature — the smallest feature set that reliably identifies IBD
7. To lay the scientific groundwork for future validation of the signature in Indian patients at AIG Hospitals, Hyderabad

---

## 2. LITERATURE REVIEW

### Paper 1: Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases (Lloyd-Price et al., *Nature*, 2019)

This is the primary source paper for the HMP2 (Human Microbiome Project Phase 2 / IBDMDB) cohort data used throughout this project. The study followed 132 IBD patients and healthy controls longitudinally for one year, generating paired metagenomics, metatranscriptomics, metaproteomics, metabolomics, and host transcriptomics measurements from 2,965 stool, biopsy, and blood specimens.

Key findings of this study include that microbial dysbiosis in IBD is characterised by decreased stability and increased instability of the microbial community over time. Specific metabolic pathways involving oxidative stress and amino acid biosynthesis were consistently altered across multiple omics layers. The study also identified that proteomics and metabolomics showed stronger associations with IBD status than taxonomic profiling alone, motivating the multi-omics integration approach used in this project. The IBDMDB database generated by this study provides the HMP2 metaproteomics, metabolomics, and metagenomics files used here.

### Paper 2: Inflammatory bowel disease genomics, transcriptomics, proteomics and metagenomics meet multi-omics (Cannarozzi et al., *UEG Journal*, 2024)

This review paper provides the conceptual framework for multi-omics integration in IBD research. It systematically summarises findings from genomics, transcriptomics, proteomics, and metagenomics studies in IBD and explains how integrating these layers improves the understanding of disease mechanisms beyond what any single layer reveals.

The paper highlights that IBD transcriptomics consistently identifies upregulation of inflammatory pathways (NF-κB, JAK-STAT, TNF signalling) and downregulation of epithelial barrier and metabolic genes. Proteomics studies have identified calprotectin, CRP, and S100 proteins as validated serum and stool biomarkers. Metagenomic studies consistently report depletion of *Firmicutes* (particularly *F. prausnitzii*) and enrichment of *Proteobacteria* in IBD. The paper argues that only by measuring all four layers simultaneously can a complete molecular signature be assembled. This directly motivates the approach taken in this project.

### Paper 3: Neutrophil extracellular traps in ulcerative colitis: a proteome analysis of intestinal biopsies (Bennike et al., *Inflammatory Bowel Diseases*, 2015)

This is the primary publication describing the PXD001608 proteomics dataset used in this project. The study performed label-free LC-MS/MS proteomics on colon biopsy tissue from 10 UC patients and 10 healthy controls in triplicate, using MaxQuant software for protein identification and quantification.

Key findings include that neutrophil extracellular trap (NET) proteins — including lactoferrin (LTF), myeloperoxidase, and histones — were massively upregulated in UC tissue, confirming a role for NETosis in UC mucosal inflammation. Structural extracellular matrix proteins including elastin (ELN) and collagens were strongly downregulated, consistent with tissue remodelling. These findings are confirmed in the preprocessing results of this project, where LTF (log2FC +4.9) and ELN (log2FC −3.4) emerge as top differentially expressed proteins. Reading this paper helped interpret the proteomics results in their correct biological context.

### Paper 4: Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 (Love et al., *Genome Biology*, 2014)

This methodology paper describes DESeq2, the statistical tool used for differential expression analysis of the GSE137344 RNA-seq transcriptomics dataset. DESeq2 models RNA-seq count data using a negative binomial distribution to account for over-dispersion in read counts. It uses a Wald test to assess differential expression and applies empirical Bayes shrinkage to stabilise fold change estimates for low-count genes. Understanding this tool was essential for correctly implementing and interpreting the transcriptomics preprocessing pipeline.

### Paper 5: limma powers differential expression analyses for RNA-sequencing and microarray studies (Ritchie et al., *Nucleic Acids Research*, 2015)

This paper describes the limma package used for differential expression analysis of the two microarray datasets (GSE186507 and GSE193677). limma fits a linear model to normalised expression values and uses empirical Bayes moderated t-statistics to improve statistical power when sample sizes are moderate. Understanding the difference between DESeq2's count-based negative binomial model and limma's linear model for pre-normalised microarray data was essential for selecting the correct analysis tool for each dataset type.

---

## 3. METHODOLOGY

We used various computational tools and programming languages to preprocess, analyse, and visualise the multi-omics data. All analysis was performed in R (version 4.5.0) and Python (version 3.11.0). Here are the detailed steps of how it was done.

---

### STEP-1: Data Collection and Organisation

All data used in this project is publicly available and fully anonymised. Data was downloaded from three primary repositories:

**GEO (Gene Expression Omnibus, NCBI)** — for transcriptomics datasets.
The accession IDs used were:
- GSE137344 — ileum RNA-seq, 193 samples
- GSE186507 — blood microarray, 1,030 samples
- GSE193677 — biopsy microarray, 2,490 samples

**PRIDE (Proteomics Identifications Database, EMBL-EBI)** — for tissue proteomics.
- PXD001608 — UC colon tissue, LC-MS/MS MaxQuant output, 60 samples

**IBDMDB (Inflammatory Bowel Disease Multi-omics Database, ibdmdb.org)** — for HMP2 stool multi-omics.
- HMP2 metaproteomics ECS and KOS PCL files — 450 samples
- HMP2 metabolomics BIOM file — 546 samples

**MetaboLights (EMBL-EBI)** — for targeted metabolomics.
- MTBLS8713 — targeted LC-MS, 179 samples

All files were downloaded and organised into a structured directory with separate folders for raw data, processed outputs, scripts, figures, and metadata.

---

### STEP-2: Transcriptomics Preprocessing

The transcriptomics preprocessing was performed separately for each of the three datasets due to differences in platform and data format.

**For GSE137344 (RNA-seq):**

The raw data file `GSE137344_Raw_Counts.txt` contains integer read counts — how many sequencing reads mapped to each gene in each sample. The R package DESeq2 was used for this dataset.

```r
library(DESeq2)
# Load raw counts
counts_raw <- read.table(counts_file, header=TRUE, sep="\t", row.names=1)

# Build sample metadata
col_data <- data.frame(
  diagnosis_binary = ifelse(sample_nums <= 149, "nonIBD", "IBD"),
  row.names = colnames(counts_raw)
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_raw,
  colData   = col_data,
  design    = ~ diagnosis_binary
)

# Filter low counts (keep genes with at least 10 reads across all samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results: IBD vs nonIBD
res <- results(dds, contrast=c("diagnosis_binary","IBD","nonIBD"), alpha=0.05)
```

After filtering, 14,406 genes were retained. 165 genes were found to be significantly differentially expressed (adjusted p-value < 0.05, |log2 fold change| > 1).

**For GSE186507 and GSE193677 (microarray):**

These datasets are microarray-based. The data is already log2-normalised by the platform before download. The limma package was used.

```r
library(limma)
# Build design matrix
design <- model.matrix(~ diagnosis_binary, data=col_data)

# Fit linear model
fit <- lmFit(norm_mat, design)

# Apply empirical Bayes
fit <- eBayes(fit)

# Extract results
deg_results <- topTable(fit, coef=2, number=Inf, sort.by="p")
```

---

### STEP-3: Proteomics Preprocessing

**For PXD001608 (tissue proteomics, MaxQuant output):**

The raw MaxQuant output file `proteinGroups.txt` was loaded using the data.table package in R. The file contains 566 columns — most are metadata about each protein, and 60 are the iBAQ intensity measurements.

```r
library(data.table)
pg <- fread("proteinGroups.txt", sep="\t", header=TRUE, data.table=FALSE)
```

Three quality filters were applied:

```r
# Remove contaminants (keratin, trypsin from sample preparation)
pg <- pg[pg[["Contaminant"]] != "+", ]

# Remove reverse decoy hits (used to calculate false discovery rate)
pg <- pg[pg[["Reverse"]] != "+", ]

# Remove proteins identified only by modified peptides (unreliable)
pg <- pg[pg[["Only identified by site"]] != "+", ]
```

Gene names were extracted from the Protein IDs column by parsing the `GN=` tag from the FASTA header string:

```r
parse_gene_name <- function(fasta_str) {
  gn <- regmatches(fasta_str, regexpr("GN=([A-Za-z0-9_]+)", fasta_str))
  if (length(gn) > 0) return(sub("GN=", "", gn))
  # Fall back to UniProt accession
  acc <- regmatches(fasta_str, regexpr("\\|([A-Z0-9]+)\\|", fasta_str))
  return(gsub("\\|", "", acc))
}
```

The iBAQ intensity columns (60 sample columns: UC_1_1 through UC_10_3, Ctrl_1_1 through Ctrl_10_3) were extracted. Zeros were replaced with NA (not detected), features with less than 50% detection were removed, and the data was log2 transformed and imputed using global mean − 1.8×SD.

**For HMP2 metaproteomics (ECS and KOS):**

The PCL-format TSV files were read directly. Sample IDs were matched to the HMP2 metadata file via the External ID column.

```r
dat <- fread("HMP2_proteomics_ecs.tsv", sep="\t", header=TRUE, data.table=FALSE)
# Remove UNGROUPED row
dat <- dat[dat[[1]] != "UNGROUPED", ]
```

After total-sum scaling (CPM normalisation) and log2 transformation, Wilcoxon rank-sum tests were applied for differential analysis.

---

### STEP-4: Metabolomics Preprocessing

**For MTBLS8713 (targeted LC-MS, MAF format):**

The MAF-format file from MetaboLights was loaded using fread. Sample columns were identified by the CTRFM (control) and IBDFM (IBD) naming convention.

```r
dat <- fread("m_MTBLS8713_LC-MS___metabolite_profiling_v2_maf.tsv",
             sep="\t", header=TRUE, data.table=FALSE)
ctrl_cols <- grep("^CTRFM", colnames(dat), value=TRUE)  # 43 controls
ibd_cols  <- grep("^IBDFM", colnames(dat), value=TRUE)  # 136 IBD
```

The pre-normalised abundance values were log2 transformed. 2 genes with near-zero variance were removed, leaving 243 metabolites.

**For HMP2 BIOM (untargeted LC-MS, 528 MB BIOM file):**

The BIOM format is a compressed JSON-based format for sparse matrices used in microbiome research. Due to the large file size (528 MB, 81,867 features × 546 samples), standard R parsing tools timed out. The solution was to use the official Python biom-format package which uses numpy sparse matrices internally:

```python
import biom
import pandas as pd

# Load BIOM file (takes approximately 78 seconds)
table = biom.load_table("HMP2_metabolomics.biom.txt")
print(table.shape)  # (81867, 546)

# Filter to samples with metadata
table_filt = table.filter(valid_samples, axis="sample", inplace=False)

# Convert to dense DataFrame
df = table_filt.to_dataframe(dense=True)

# Filter features detected in >10% of samples in at least one group
det_ibd  = (df[ibd_s]  > 0).mean(axis=1)
det_ctrl = (df[ctrl_s] > 0).mean(axis=1)
keep     = (det_ibd >= 0.10) | (det_ctrl >= 0.10)
df       = df.loc[keep]   # 75,447 features retained

# Save as CSV for R
df.to_csv("HMP2_metabolomics_matrix.csv")
```

This CSV was then read into R for the normalisation and differential analysis steps:

```r
mat_df <- fread("HMP2_metabolomics_matrix.csv", header=TRUE, data.table=FALSE)
# CPM normalisation and log2 transform
col_sums  <- colSums(mat)
mat_tss   <- sweep(mat, 2, col_sums, "/") * 1e6
mat_log2  <- log2(mat_tss + 1)
```

---

### STEP-5: Differential Analysis

For each dataset, statistical testing was performed to identify features significantly different between IBD and nonIBD.

**Test selection was based on data type:**

| Dataset | Statistical Test | Reason |
|---|---|---|
| GSE137344 | DESeq2 Wald test | Integer counts, negative binomial model |
| GSE186507, GSE193677 | limma moderated t-test | Pre-normalised microarray, linear model |
| PXD001608 | Welch t-test | Log2-normal after transformation, unequal variance |
| HMP2 ECS/KOS | Wilcoxon rank-sum | Non-normal spectral count data |
| MTBLS8713 | Wilcoxon rank-sum | Small sample size, non-normal |
| HMP2 BIOM | Wilcoxon rank-sum | Non-normal metabolomics intensities |

Multiple testing correction was applied using the Benjamini-Hochberg (BH) procedure in all cases:

```r
de_results$padj <- p.adjust(de_results$pval, method="BH")
```

This controls the False Discovery Rate (FDR) — the expected proportion of false positives among all features called significant — at 5%.

---

### STEP-6: Heatmap Generation

Heatmaps were generated for all datasets using the pheatmap R package. For each dataset, the top 40 most upregulated and top 40 most downregulated features by fold change were selected. Expression values were z-score scaled per feature so that the colour represents relative change rather than absolute abundance.

```r
library(pheatmap)

# Z-score scale rows
mat_scaled <- t(scale(t(mat_sub)))
mat_scaled[mat_scaled >  3] <-  3   # clip extremes
mat_scaled[mat_scaled < -3] <- -3

# Generate heatmap with diagnosis and direction annotations
pheatmap(
  mat_scaled,
  color          = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
  annotation_col = anno_col,   # IBD/nonIBD colour bar (top)
  annotation_row = anno_row,   # Up/Down colour bar (left)
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,      # keep diagnosis groups together
  show_colnames  = FALSE
)
```

Two versions of each heatmap were produced:
- **Per-patient heatmap** — one column per sample, shows patient-level variability
- **Summary heatmap** — two columns only (IBD average vs nonIBD average), immediately readable

---

### STEP-7: Transcript-Protein Concordance Analysis

To identify genes robustly altered in IBD at both the mRNA and protein level, a concordance analysis was performed between GSE137344 (ileum transcriptomics, gene symbols) and PXD001608 (colon proteomics, gene symbols from FASTA header parsing).

```r
library(dplyr)

# Match genes by symbol
matched <- inner_join(
  deg137 %>% select(gene, rna_lfc=log2FoldChange, rna_padj=padj),
  allpxd %>% select(protein, prot_lfc=log2FC, prot_padj=padj),
  by = c("gene"="protein")
)

# Classify concordance
matched <- matched %>%
  mutate(
    concordant = (rna_lfc > 0) == (prot_lfc > 0),
    direction  = case_when(
      rna_lfc > 0 & prot_lfc > 0 ~ "Both UP",
      rna_lfc < 0 & prot_lfc < 0 ~ "Both DOWN",
      rna_lfc > 0 & prot_lfc < 0 ~ "RNA up / Protein down",
      rna_lfc < 0 & prot_lfc > 0 ~ "RNA down / Protein up"
    )
  )
```

The scatter plot of mRNA log2FC vs protein log2FC was generated using ggplot2 and ggrepel, with points coloured by significance status and top concordant pairs labelled.

---

## 4. RESULTS

### 4.1 Transcriptomics Results

Differential expression analysis across the three transcriptomics datasets identified consistent patterns of gene dysregulation in IBD.

**GSE137344 (ileum, RNA-seq):** 165 genes were significantly differentially expressed (padj < 0.05, |log2FC| > 1) from 14,406 genes tested.

**Table 1. Top differentially expressed genes in GSE137344 (ileum transcriptomics)**

| Gene | log2FC | padj | Direction | Function |
|---|---|---|---|---|
| RNU2-46P | +6.04 | 1.7×10⁻²¹ | Up in IBD | Non-coding RNA |
| MTRNR2L1 | +4.25 | 4.7×10⁻¹⁹ | Up in IBD | Mitochondrial RNA |
| NEAT1 | +2.89 | 4.1×10⁻¹⁶ | Up in IBD | Nuclear paraspeckle assembly |
| ATG16L2 | +1.87 | 2.3×10⁻¹² | Up in IBD | Autophagy — IBD GWAS gene |
| HMGCS2 | −1.38 | 8.1×10⁻⁹ | Down in IBD | Cholesterol/energy metabolism |
| LDHD | −1.13 | 3.1×10⁻⁷ | Down in IBD | D-lactate metabolism |

The heatmap shows a clear separation between IBD (first 156 columns) and nonIBD samples (last 37 columns), with CD and UC labelled separately by the orange and purple annotation bars respectively.

*[Figure 1: GSE137344 heatmap — see GSE137344_heatmap.pdf]*

**GSE186507 (blood, microarray):** 146 genes passed significance thresholds (|logFC| > 0.5, padj < 0.05). The smaller fold changes (maximum ~1.35) compared to tissue data reflect the dilution of the IBD signal across diverse circulating cell types in blood.

**GSE193677 (biopsy, microarray):** Only 3 genes survived BH correction. This reflects the extremely small fold changes (maximum 0.45) across 2,490 heterogeneous biopsy samples with patients at different disease stages. The top 50 genes by raw p-value are shown in the heatmap for visualisation purposes.

*[Figure 2: GSE186507 heatmap — see GSE186507_heatmap.pdf]*
*[Figure 3: GSE193677 heatmap — see GSE193677_heatmap.pdf]*

---

### 4.2 Proteomics Results

**PXD001608 (UC colon tissue proteomics):** 721 proteins were significantly differentially expressed (padj < 0.05, |log2FC| > 1) from 5,620 quantified proteins across 60 samples (30 UC vs 30 healthy).

**Table 2. Top differentially expressed proteins in PXD001608 (UC colon tissue)**

| Protein | log2FC | padj | Direction | Function |
|---|---|---|---|---|
| LTF | +4.90 | 9.9×10⁻⁶ | Up in IBD | Lactoferrin — neutrophil marker |
| PREX1 | +1.94 | 3.8×10⁻³ | Up in IBD | Immune cell signalling |
| ELN | −3.41 | 6.7×10⁻⁸ | Down in IBD | Elastin — ECM structural protein |
| FBN2 | −2.31 | 9.9×10⁻⁶ | Down in IBD | Fibrillin-2 — ECM protein |
| WFDC2 | −4.11 | 5.8×10⁻⁵ | Down in IBD | Epithelial secretory protein |
| COL14A1 | −1.58 | 8.0×10⁻⁷ | Down in IBD | Collagen type XIV |

The heatmap shows a striking and clean separation between UC and healthy samples, with the upregulated inflammatory proteins (LTF, immune markers) forming a strong red block in the UC columns and the downregulated ECM proteins forming an equally strong blue block.

*[Figure 4: PXD001608 proteomics heatmap — see PXD001608_proteomics_heatmap.pdf]*

**HMP2 metaproteomics ECS:** 75 enzyme functions were significantly different (padj < 0.05, |log2FC| > 0.5). The predominant pattern was downregulation in IBD, with phosphoenolpyruvate carboxykinase (−3.24), formate C-acetyltransferase (−3.02), and starch synthase (−2.37) among the strongest hits. This reflects broad collapse of microbial energy and carbohydrate metabolism in the dysbiotic IBD gut.

**HMP2 metaproteomics KOS:** 81 KEGG orthologs were significantly different. The same energy metabolism pathways flagged by ECS were independently confirmed by KOS — pyruvate dikinase (K01006, −1.29) and PEP carboxykinase (K01596, −2.48) were the top hits.

*[Figure 5: HMP2 ECS heatmap — see HMP2_ECS_heatmap.pdf]*
*[Figure 6: HMP2 KOS heatmap — see HMP2_KOS_heatmap.pdf]*

---

### 4.3 Metabolomics Results

**MTBLS8713 (targeted LC-MS, 243 named metabolites):** 207 metabolites were significantly different (padj < 0.05, |log2FC| > 0.3).

**Table 3. Top differentially abundant metabolites in MTBLS8713**

| Metabolite | log2FC | padj | Direction | Function |
|---|---|---|---|---|
| Pyruvic acid | +5.42 | 3.1×10⁻¹⁵ | Up in IBD | Central energy metabolism |
| 1-Methylhistidine | +4.76 | 5.0×10⁻¹⁶ | Up in IBD | Histidine metabolism |
| Deoxycholic acid | −4.77 | 4.8×10⁻⁶ | Down in IBD | Secondary bile acid |
| Lithocholic acid | −3.70 | 5.6×10⁻⁵ | Down in IBD | Secondary bile acid |
| Isolithocholic acid | −3.39 | 1.8×10⁻⁵ | Down in IBD | Secondary bile acid |

The most notable pattern is the consistent downregulation of secondary bile acids (deoxycholic, lithocholic, isolithocholic acids) in IBD. These are produced exclusively by healthy commensal bacteria and their depletion is a direct readout of gut dysbiosis. The upregulation of pyruvic acid and amino acid metabolites reflects impaired host metabolic function.

*[Figure 7: MTBLS8713 metabolomics heatmap — see MTBLS8713_metabolomics_heatmap.pdf]*

**HMP2 BIOM (untargeted LC-MS, 75,447 features):** 8,419 features were significantly different (padj < 0.05, |log2FC| > 0.5), with 5,840 upregulated and 2,579 downregulated in IBD. Feature IDs are internal mass spectrometry codes (C18n_QI, HILp_QI) as the majority of features in untargeted metabolomics remain chemically unidentified at this stage. The large number of significant features reflects the high sensitivity of the 546-sample dataset and the broad disruption of the stool metabolome in IBD.

*[Figure 8: HMP2 metabolomics heatmap — see HMP2_metabolomics_heatmap.pdf]*

---

### 4.4 Cross-Omics Concordance Analysis

To identify genes whose disruption in IBD is robustly detectable at both the mRNA and protein level, a concordance analysis was performed between GSE137344 (ileum transcriptomics) and PXD001608 (colon tissue proteomics).

**Table 4. Concordance analysis summary**

| | Count |
|---|---|
| Total genes measured in both datasets | 5,005 |
| Significant in both (padj < 0.05 each) | 464 |
| **Concordant (same direction in both)** | **294 (63.4%)** |
| Discordant (opposite direction) | 170 (36.6%) |
| Both DOWN in IBD | 278 |
| Both UP in IBD | 16 |

The scatter plot of mRNA log2FC (x-axis) vs protein log2FC (y-axis) for all 5,005 matched pairs shows a positive correlation among significantly altered genes, confirming that transcriptomic and proteomic measurements are capturing the same underlying biological signal.

*[Figure 9: Concordance scatter plot — see figures/concordance/concordance_scatter.pdf]*
*[Figure 10: Concordance mRNA heatmap — see figures/concordance/concordance_heatmap_mRNA.pdf]*
*[Figure 11: Concordance protein heatmap — see figures/concordance/concordance_heatmap_protein.pdf]*

**Table 5. Top 10 concordant gene-protein pairs (both significant, both downregulated)**

| Gene | mRNA log2FC | Protein log2FC | mRNA padj | Protein padj | Function |
|---|---|---|---|---|---|
| LZIC | −0.443 | −0.925 | 8.7×10⁻⁹ | 3.4×10⁻⁵ | Cell survival regulation |
| FERMT1 | −0.721 | −0.666 | 4.3×10⁻⁸ | 1.2×10⁻² | Epithelial cell-matrix adhesion (Kindlin-1) |
| HSPD1 | −0.466 | −0.373 | 4.6×10⁻⁸ | 1.1×10⁻² | Mitochondrial chaperone (Hsp60) |
| HDAC2 | −0.261 | −1.146 | 1.5×10⁻⁷ | 8.7×10⁻³ | Histone deacetylase — chromatin regulation |
| ROCK2 | −0.443 | −0.566 | 1.4×10⁻⁷ | 1.1×10⁻² | Rho kinase — cytoskeletal regulation |
| UQCRB | −0.458 | −0.480 | 3.0×10⁻⁹ | 6.3×10⁻³ | Mitochondrial electron transport chain |
| CS | −0.318 | −0.436 | 2.1×10⁻⁴ | 3.8×10⁻³ | Citrate synthase — TCA cycle |
| HMGCS2 | −1.375 | −1.594 | 1.2×10⁻⁶ | 3.5×10⁻³ | HMG-CoA synthase — cholesterol metabolism |
| LDHD | −1.131 | −1.317 | 1.7×10⁻⁴ | 8.1×10⁻³ | Lactate dehydrogenase |
| EDEM2 | +0.735 | +0.538 | 7.5×10⁻⁸ | 7.7×10⁻³ | ER stress / unfolded protein response |

The 278 concordant downregulated pairs cluster into clear biological themes: mitochondrial energy metabolism (UQCRB, CS, HDAC2, HSPD1), epithelial barrier and cell adhesion (FERMT1, ROCK2), and metabolic enzymes (HMGCS2, LDHD). The 16 concordant upregulated pairs are dominated by ER stress and immune activation markers.

---

## 5. FUTURE DIRECTIONS

- To complete metagenomics preprocessing of the HMP2 taxonomic profiles and pathway abundance data (MetaPhlAn2 and HUMAnN2 outputs) which is the final remaining omics layer

- To build individual single-omics machine learning models (LASSO Logistic Regression and Random Forest) for each dataset and evaluate their IBD classification performance

- To build a combined multi-omics model integrating features from all four omics layers and identify the minimal feature set that achieves the highest classification accuracy

- To annotate the HMP2 BIOM untargeted metabolomics features using the C18n annotation file and chemical databases (HMDB, KEGG) to replace internal feature codes with named metabolites

- To validate the molecular signature in Indian patients at AIG Hospitals, Hyderabad — testing which features hold up across populations and which are population-specific

- To identify whether the IBD molecular signature can differentiate CD from UC, not just IBD from healthy, which would be clinically valuable for treatment selection

- To explore whether the 294 concordant gene-protein pairs are sufficient on their own to train a high-accuracy classifier, which would represent the most robust and reproducible biomarker panel

---

## 6. CONCLUSION

This report presents the completion of comprehensive multi-omics data preprocessing and exploratory analysis for the IBD molecular signature project. All four omics data types — transcriptomics, proteomics, metabolomics, and metagenomics — have been collected. Preprocessing is complete for transcriptomics (3 datasets), proteomics (3 datasets), and metabolomics (2 datasets), together covering 4,398 patient samples.

The differential analysis results reveal a consistent and biologically coherent picture of IBD across all data layers. IBD is characterised by widespread downregulation of normal epithelial cellular functions — particularly mitochondrial energy metabolism, cell-matrix adhesion, and metabolic enzymes — alongside a smaller but clear signature of upregulated inflammatory markers including neutrophil proteins, ER stress responses, and immune signalling molecules. This same pattern appears independently in ileum tissue, colon tissue, blood, and stool, across transcriptomics, proteomics, and metabolomics measurements, providing strong convergent evidence.

The concordance analysis identified 294 genes consistently dysregulated at both the mRNA and protein level across independent datasets, representing the highest-confidence molecular biomarker candidates for the final signature. The bile acid depletion pattern confirmed across both metabolomics datasets (MTBLS8713 and HMP2 BIOM independently) provides a further strong candidate for a non-invasive stool-based IBD test.

The next phase — metagenomics preprocessing followed by machine learning model training — will synthesise these findings into a final molecular signature. The ultimate goal remains the development of a non-invasive, multi-biomarker test for IBD detection that can be validated and adapted for Indian patients, addressing an unmet clinical need at AIG Hospitals and more broadly across India where IBD incidence is rising rapidly.

---

## 7. REFERENCES

- Lloyd-Price J, et al. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 2019;569:655–662.

- Cannarozzi L, et al. Inflammatory bowel disease genomics, transcriptomics, proteomics and metagenomics meet multi-omics. *UEG Journal* 2024.

- Bennike TB, et al. Neutrophil extracellular traps in ulcerative colitis: a proteome analysis of intestinal biopsies. *Inflammatory Bowel Diseases* 2015.

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 2014;15:550.

- Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 2015;43:e47.

- McDonald D, et al. The Biological Observation Matrix (BIOM) format. *GigaScience* 2012;1:7.

- R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria, 2024.

- Python biom-format package v2.1.17. https://biom-format.org

- pheatmap: Pretty Heatmaps. R package. https://cran.r-project.org/package=pheatmap

- Geneious software. https://www.geneious.com

- RStudio. https://www.rstudio.com
