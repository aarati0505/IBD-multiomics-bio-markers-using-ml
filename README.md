# IBD Multi-Omics Molecular Signature Project

**BIO F266: Study Project — BITS Pilani, Hyderabad Campus**  
**Student:** Aarati Dhamele (2023B1A70641H)  
**Supervisor:** [Supervisor Name]

---

## Overview

This project integrates multiple types of biological data — transcriptomics, proteomics, metabolomics, and metagenomics — from publicly available IBD (Inflammatory Bowel Disease) datasets to identify a robust molecular signature that can distinguish IBD patients from healthy controls.

The end goal is a minimal panel of molecular markers (genes, proteins, metabolites, microbial features) that together reliably identify IBD — potentially enabling a non-invasive diagnostic test from blood or stool, with future validation planned in Indian patients at AIG Hospitals, Hyderabad.

---

## What is IBD?

Inflammatory Bowel Disease includes two forms:
- **Crohn's Disease (CD)** — can affect any part of the digestive tract
- **Ulcerative Colitis (UC)** — confined to the colon and rectum

Current diagnosis requires invasive colonoscopy and biopsy. This project aims to find molecular markers detectable from simpler sample types.

---

## Datasets Used

All data is publicly available and fully anonymised.

| Dataset | Omics Type | Source | Samples | IBD | Healthy |
|---|---|---|---|---|---|
| GSE137344 | Transcriptomics — RNA-seq (ileum) | GEO | 193 | 156 | 37 |
| GSE186507 | Transcriptomics — RNA-seq (blood) | GEO | 1,030 | 821 | 209 |
| GSE193677 | Transcriptomics — RNA-seq (colon biopsy) | GEO | 2,490 | 2,029 | 461 |
| PXD001608 | Proteomics — LC-MS/MS (colon tissue) | PRIDE | 60 | 30 | 30 |
| HMP2 MPX ECS | Metaproteomics — stool enzyme functions | IBDMDB | 450 | 328 | 122 |
| HMP2 MPX KOS | Metaproteomics — stool KEGG orthologs | IBDMDB | 450 | 328 | 122 |
| MTBLS8713 | Metabolomics — targeted LC-MS (stool) | MetaboLights | 179 | 136 | 43 |
| HMP2 BIOM | Metabolomics — untargeted LC-MS (stool) | IBDMDB | 546 | 411 | 135 |
| HMP2 MGX | Metagenomics — taxonomic + pathway profiles | IBDMDB | 96 | 72 | 24 |

**Total: 4,398+ samples across 9 datasets**

---

## Project Structure

```
IBD-MultiOmics/
│
├── scripts/                          # All analysis scripts
│   ├── preprocess_GSE137344.R        # Ileum RNA-seq preprocessing (DESeq2)
│   ├── preprocess_GSE186507.R        # Blood RNA-seq preprocessing (limma)
│   ├── preprocess_GSE193677.R        # Biopsy RNA-seq preprocessing (limma)
│   ├── preprocess_proteomics_MaxQuant.R   # Tissue proteomics (PXD001608)
│   ├── preprocess_HMP2_metaproteomics.R   # Stool metaproteomics (ECS + KOS)
│   ├── preprocess_MTBLS8713.R        # Targeted metabolomics
│   ├── preprocess_HMP2_metabolomics.R     # Untargeted metabolomics (BIOM)
│   ├── parse_biom.py                 # Python parser for 528MB BIOM file
│   ├── transcriptome_proteome_concordance.R   # Cross-omics concordance
│   ├── concordance_colon_matched.R   # Tissue-matched concordance
│   ├── metabolomics_concordance.R    # Metabolomics concordance analysis
│   ├── generate_heatmaps.R           # Per-patient heatmaps (all datasets)
│   ├── generate_summary_heatmaps.R   # Group-average summary heatmaps
│   ├── generate_metabolomics_heatmaps.R   # Metabolomics heatmaps
│   ├── generate_HMP2_metabolomics_heatmaps.R
│   ├── check_known_biomarkers.R      # Validate against literature biomarkers
│   └── build_ensembl_map.py          # Ensembl ID -> gene symbol conversion
│
├── figures/
│   ├── heatmaps/                     # Per-patient heatmaps (PDF)
│   │   ├── summary/                  # Group-average summary heatmaps (PDF)
│   │   └── *.pdf
│   └── concordance/                  # Concordance scatter plots and heatmaps
│       └── *.pdf
│
├── data/
│   ├── hmp2_metadata_clean.csv       # Cleaned HMP2 patient metadata
│   └── hmp2_omics_per_patient.csv
│
├── processed/                        # Processed outputs (not included in repo)
│   ├── GSE137344/
│   ├── GSE186507/
│   ├── GSE193677/
│   ├── proteomics_MaxQuant/
│   ├── HMP2_metaproteomics/
│   ├── HMP2_metabolomics/
│   └── MTBLS8713_metabolomics/
│
├── literature/                       # Key reference papers
├── IBD_MultiOmics_Project_Summary.md # Plain-language project summary
└── README.md
```

---

## Analysis Pipeline

Each dataset follows a six-step preprocessing pipeline:

```
Raw data
  ↓ Quality filtering      — remove contaminants, low-detection features
  ↓ Normalisation          — DESeq2 / TMM / Total Sum Scaling (CPM)
  ↓ Log₂ transform         — compress intensity range
  ↓ Imputation             — fill missing values (mean − 1.8 SD)
  ↓ Differential analysis  — DESeq2 / limma / Welch t-test / Wilcoxon
  ↓ BH correction          — control false discovery rate at 5%
  ↓ ML matrix              — top features → samples × features table
```

Statistical test selection by data type:

| Data type | Test | Reason |
|---|---|---|
| RNA-seq raw counts | DESeq2 Wald test | Negative binomial model for count data |
| Microarray (normalised) | limma moderated t-test | Linear model for continuous data |
| Tissue proteomics | Welch t-test | Log-normal after transformation |
| Stool metaproteomics | Wilcoxon rank-sum | Non-parametric, non-normal counts |
| Metabolomics | Wilcoxon rank-sum | Non-normal intensities |

---

## Key Results

### Significant features identified

| Dataset | Features tested | Significant | Direction |
|---|---|---|---|
| GSE137344 (ileum) | 14,406 genes | 165 | 151 up, 14 down in IBD |
| GSE186507 (blood) | 24,835 genes | 146 | mixed |
| PXD001608 (colon) | 5,620 proteins | 721 | 509 up, 212 down in IBD |
| HMP2 ECS | 177 enzymes | 75 | mostly down in IBD |
| HMP2 KOS | 246 KOs | 81 | mostly down in IBD |
| MTBLS8713 | 243 metabolites | 207 | 186 up, 21 down in IBD |
| HMP2 BIOM | 75,447 features | 8,419 | 5,840 up, 2,579 down in IBD |

### Cross-omics concordance

Of 5,005 genes measured in both transcriptomics (GSE137344) and proteomics (PXD001608):
- **464** were statistically significant in both datasets independently (padj < 0.05)
- **294 were concordant** — same direction of change in both mRNA and protein (63.4%)
- 278 concordant pairs are **Both DOWN** in IBD — widespread loss of normal cellular function
- 16 concordant pairs are **Both UP** — ER stress and immune activation

Top concordant genes include **FERMT1** (epithelial adhesion), **HDAC2** (chromatin regulation), **ROCK2** (cytoskeletal organisation), **HMGCS2** (energy metabolism), and **S100A9** (calprotectin subunit — clinically validated IBD marker).

### Validated clinical biomarkers found

The following clinically validated IBD biomarkers were independently identified in this analysis:

- **LTF (Lactoferrin)** — +4.90 log2FC in proteomics — validated stool IBD test
- **S100A8/S100A9 (Calprotectin)** — significant in proteomics — FDA-approved stool IBD test
- **LCN2 (Lipocalin-2)** — significant in proteomics — validated stool IBD marker
- **MMP9** — +3.28 in proteomics — neutrophil gelatinase, known IBD marker
- **ATG16L2** — +1.81 in RNA-seq — autophagy gene, IBD GWAS candidate
- **NEAT1** — +1.40 in RNA-seq — lncRNA associated with IBD inflammation
- **Deoxycholic acid / Lithocholic acid** — strongly depleted in both metabolomics datasets — secondary bile acid loss, one of the most replicated IBD metabolomics findings

### Biological interpretation

The consistent finding across all data types is that **IBD is predominantly a disease of loss**:
- Host tissue loses structural integrity (ECM proteins ELN, COL1A1, COL1A2 depleted)
- Host cells lose normal metabolic function (HMGCS2, UQCRB, mitochondrial genes depleted)
- Gut microbiome loses metabolic activity (pyruvate, glutamine processing enzymes depleted)
- Secondary bile acids disappear (microbial conversion pathway disrupted)

The smaller set of upregulated features — LTF, S100A8/A9, MMP9, pyruvic acid — represents the inflammatory and stress response.

---

## Technologies Used

**R packages:** DESeq2, limma, pheatmap, ggplot2, ggrepel, dplyr, data.table, Matrix, biomformat, RColorBrewer

**Python packages:** biom-format, pandas, numpy, requests

**Data sources:** NCBI GEO, EMBL-EBI PRIDE, MetaboLights, IBDMDB (ibdmdb.org)

---

## How to Reproduce

### Requirements
- R ≥ 4.0 with Bioconductor packages: DESeq2, limma
- R CRAN packages: pheatmap, ggplot2, ggrepel, dplyr, data.table, Matrix, RColorBrewer
- Python ≥ 3.8 with: biom-format, pandas, numpy

### Data download
Raw data files are not included in this repository due to size (total ~6 GB).  
Download instructions and accession numbers:

| Dataset | Repository | Accession |
|---|---|---|
| GSE137344 | NCBI GEO | GSE137344 |
| GSE186507 | NCBI GEO | GSE186507 |
| GSE193677 | NCBI GEO | GSE193677 |
| PXD001608 | EMBL-EBI PRIDE | PXD001608 |
| HMP2 multi-omics | IBDMDB | ibdmdb.org |
| MTBLS8713 | MetaboLights | MTBLS8713 |

### Run order
```
# Transcriptomics
Rscript scripts/preprocess_GSE137344.R
Rscript scripts/preprocess_GSE186507.R
Rscript scripts/preprocess_GSE193677.R

# Proteomics
Rscript scripts/preprocess_proteomics_MaxQuant.R
Rscript scripts/preprocess_HMP2_metaproteomics.R

# Metabolomics
python scripts/parse_biom.py                    # parse HMP2 BIOM first
Rscript scripts/preprocess_HMP2_metabolomics.R
Rscript scripts/preprocess_MTBLS8713.R

# Concordance analysis
Rscript scripts/transcriptome_proteome_concordance.R
Rscript scripts/concordance_colon_matched.R

# Visualisation
Rscript scripts/generate_heatmaps.R
Rscript scripts/generate_summary_heatmaps.R
Rscript scripts/generate_metabolomics_heatmaps.R
Rscript scripts/check_known_biomarkers.R
```

---

## Future Work

- Metagenomics preprocessing (HMP2 taxonomic profiles + HUMAnN2 pathway abundances)
- Multi-omics machine learning model (LASSO Logistic Regression + Random Forest)
- Molecular signature extraction — minimal feature set for IBD classification
- Validation in Indian IBD patients at AIG Hospitals, Hyderabad

---

## Citation

If you use this code or analysis, please cite the primary data sources:

- Lloyd-Price J, et al. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 2019;569:655–662.
- Bennike TB, et al. Neutrophil extracellular traps in ulcerative colitis. *Inflammatory Bowel Diseases* 2015.

---

## License

This project uses only publicly available, anonymised research data.  
Code is available under the MIT License.

---

*BITS Pilani, Hyderabad Campus — BIO F266 Study Project — 2025*
