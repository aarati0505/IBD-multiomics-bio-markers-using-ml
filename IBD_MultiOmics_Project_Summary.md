# IBD Multi-Omics ML Pipeline — Simple Summary

---

## What Is This Project?

We are building a machine learning model that can look at a patient's biology
and predict whether they have IBD (Inflammatory Bowel Disease) or are healthy.

IBD has two types:
- Crohn's Disease (CD) — can affect any part of the digestive tract
- Ulcerative Colitis (UC) — affects the large intestine

The model looks at 4 types of biological data from the same patient:
1. Genes that are switched on or off (Transcriptomics)
2. Proteins present in the body (Proteomics)
3. Small chemical molecules (Metabolomics)
4. Gut bacteria (Metagenomics)

The end goal is a small list of biological markers — genes, proteins, molecules,
bacteria — that together reliably identify IBD. This is called a Molecular Signature.

---

## Why This Matters

- IBD affects millions of people worldwide
- Diagnosis currently requires invasive procedures like colonoscopy
- A blood or stool test based on this signature could detect IBD non-invasively
- This project will eventually be adapted for Indian patients via AIG Hospitals

---

## The Data We Are Using

### Where the data comes from
All data is from publicly available American/Western research studies.
Indian data does not yet exist at the scale needed — collecting it is the next phase.

| Data Type | Dataset Name | What It Is | Number of Samples |
|---|---|---|---|
| Gene expression (biopsy) | GSE193677 | Intestinal tissue from IBD + healthy adults | 2,490 samples |
| Gene expression (blood) | GSE186507 | Blood samples from IBD + healthy adults | 1,030 samples |
| Gene expression (ileum) | GSE137344 | Small intestine tissue | 195 samples |
| Tissue proteomics | PXD001608 | Proteins found in UC colon biopsy tissue (MaxQuant) | 60 samples (30 UC / 30 Ctrl) |
| Stool metaproteomics (ECS) | HMP2 MPX ECS | Enzyme functions detected in stool — microbial activity | 450 samples |
| Stool metaproteomics (KOS) | HMP2 MPX KOS | KEGG orthologs detected in stool — pathway activity | 450 samples |
| Gut bacteria | HMP2 metagenomics | Species of bacteria in stool | 1,638 samples |
| Metabolites | HMP2 metabolomics | Chemical molecules in stool | 546 samples |
| Patient info | HMP2 metadata | Links all data to diagnosis | 5,533 entries |

---

## Step 1 — Data Preprocessing (COMPLETED for Transcriptomics)

### What is preprocessing?
Raw biological data is messy and noisy. Preprocessing cleans it and
converts it into a format a computer can learn from.

Think of it like this: before you cook a meal, you wash vegetables,
chop them into equal sizes, and remove bad parts. Preprocessing does
the same for biological data.

### What we did for gene expression data (Transcriptomics)

**GSE137344 (195 samples)**
- Loaded 35,857 genes measured across 195 patients
- Removed 21,451 genes that were barely detected (noise)
- Kept 14,406 meaningful genes
- Normalized the data so all patients are on the same scale
- Found 165 genes that are significantly different between IBD and healthy
- Selected the top 500 most informative genes for the model
- Saved as a clean table: 193 patients × 500 genes ✅

**GSE186507 (1,030 samples — blood)**
- Loaded 27,595 genes across 1,030 patients
- Downloaded patient diagnosis labels (CD / UC / Control) from public database
- 821 IBD patients, 209 healthy controls
- Removed low-information genes, kept 24,835
- Found differentially expressed genes between IBD and healthy
- Selected top 500 genes for model
- Saved as clean table: 1,030 patients × 500 genes ✅

**GSE193677 (2,490 samples — biopsy)**
- Loaded 25,969 genes across 2,490 patients
- Downloaded diagnosis labels: 2,029 IBD, 461 healthy
- Removed low-information genes, kept 23,372
- Selected top 500 genes for model
- Saved as clean table: 2,490 patients × 500 genes ✅

### Output of preprocessing
Three clean tables where:
- Every row = one patient
- Every column = one gene
- Every cell = how active that gene is in that patient
- Every row has a label: IBD or Healthy

---

## Step 2 — Preprocessing Remaining Data (IN PROGRESS)

### Proteomics ✅ COMPLETED

**PXD001608 — UC Colon Biopsy Tissue (MaxQuant iBAQ)**
- 6,465 protein groups loaded → 6,203 after contaminant/decoy removal
- 5,620 proteins retained after detection-rate filter (>50% in at least one group)
- 60 samples: 30 UC / 30 healthy controls (triplicates for 20 patients)
- 721 significant differentially expressed proteins (padj < 0.05, |log2FC| > 1)
- Top hits: LTF (+4.9 in UC), ELN (−3.4), COL14A1 (−1.6), FBN2 (−2.3)
- Saved: processed/proteomics_MaxQuant/ ✅

**HMP2 Stool Metaproteomics — ECS (Enzyme Functions)**
- 910 enzyme functions → 177 retained after detection-rate filter (>30%)
- 450 samples: 328 IBD (204 CD + 124 UC) / 122 nonIBD
- 75 significant features (padj < 0.05, |log2FC| > 0.5)
- Top hits: Phosphoenolpyruvate carboxykinase (−3.2), Formate C-acetyltransferase (−3.0)
- Key pattern: broad downregulation of microbial metabolic enzymes in IBD
- Saved: processed/HMP2_metaproteomics/ ✅

**HMP2 Stool Metaproteomics — KOS (KEGG Orthologs)**
- 1,823 KEGG orthologs → 246 retained after detection-rate filter (>30%)
- 450 samples: 328 IBD / 122 nonIBD
- 81 significant features (padj < 0.05, |log2FC| > 0.5)
- Top hits: K01006 pyruvate dikinase (−1.3), K01596 PEP carboxykinase (−2.5)
- Cross-validated with ECS: same metabolic pathways flagged independently
- Saved: processed/HMP2_metaproteomics/ ✅

- Metagenomics (bacteria) — clean species abundance data 🔲 Next
- Metabolomics (molecules) — clean metabolite data 🔲 Upcoming

---

## Step 3 — Building the ML Model (UPCOMING)

### How the model learns
1. We show it thousands of examples: "this patient has IBD" / "this patient is healthy"
2. The model finds patterns — combinations of genes/proteins/bacteria that
   consistently appear in IBD patients but not in healthy ones
3. We test it on patients it has never seen before
4. If it classifies correctly >85% of the time, the model is considered good

### What model we will use
LASSO Logistic Regression — a type of model that:
- Automatically selects only the most important features
- Ignores irrelevant ones (sets them to zero)
- Is transparent — you can see exactly which features it picked and why
- Works well with biological data

We will also use Random Forest as a backup to double-check the results.

---

## Step 4 — Extracting the Molecular Signature

After training, the model will have selected a small set of features.
These are the ones that matter most for identifying IBD.

The expected signature looks something like:

| Layer | Example Markers |
|---|---|
| Genes | S100A8, S100A9, NEAT1, ATG16L2 |
| Proteins | Calprotectin, CRP, CXCL8 |
| Bacteria | ↓ Faecalibacterium prausnitzii, ↑ E. coli |
| Metabolites | ↓ Butyrate, ↓ Short-chain fatty acids |

Features that appear across ALL datasets = strongest candidates for the signature.

---

## Step 5 — Indian Validation (AIG Hospitals)

- The signature discovered on US data will be tested on Indian patients
- Features that hold up = universal IBD markers
- Features that change = population-specific differences
- New Indian-specific markers may emerge
- Final output = an IBD molecular signature adapted for Indian patients

---

## Current Status

| Step | Status |
|---|---|
| Collect transcriptomics data | ✅ Done |
| Collect proteomics data | ✅ Done |
| Collect metagenomics data | ✅ Done |
| Collect metabolomics data | ✅ Done |
| Preprocess transcriptomics | ✅ Done — 3 datasets (GSE137344, GSE186507, GSE193677) |
| Preprocess proteomics — tissue (PXD001608) | ✅ Done — 5,620 proteins, 60 samples, 721 DEPs |
| Preprocess proteomics — stool ECS (HMP2) | ✅ Done — 177 features, 450 samples, 75 sig. |
| Preprocess proteomics — stool KOS (HMP2) | ✅ Done — 246 features, 450 samples, 81 sig. |
| Preprocess metagenomics | 🔲 Next |
| Preprocess metabolomics | 🔲 Upcoming |
| Build ML model | 🔲 Upcoming |
| Extract molecular signature | 🔲 Upcoming |
| Indian validation (AIG) | 🔲 Future phase |

---

## Key Numbers So Far

- Total gene expression samples processed: **3,713**
- Genes selected for ML per transcriptomics dataset: **500**
- IBD patients across transcriptomics: **3,015** | Healthy: **698**
- Tissue proteomics (PXD001608): **60 samples**, 5,620 proteins, **721 DEPs**
- Stool metaproteomics ECS: **450 samples**, 177 enzyme features, **75 significant**
- Stool metaproteomics KOS: **450 samples**, 246 KEGG features, **81 significant**
- Total proteomics samples processed: **510**

---

*This project uses only publicly available, anonymized research data.
No patient names or personal information are used at any stage.*
