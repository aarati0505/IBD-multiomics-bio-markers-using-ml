# ============================================================
# PREPROCESSING: HMP2 Metadata
# Master linking file — connects patients to all omics
# Input : hmp2_metadata_2018-08-20.csv
# Output: data/hmp2_metadata_clean.rds
# ============================================================

library(dplyr)
library(tidyr)

dest <- "C:/Users/Aarati/Desktop/IBD datasets"
dir.create(file.path(dest, "data"), showWarnings = FALSE)

cat("=== HMP2 Metadata Preprocessing ===\n\n")

# ----------------------------------------------------------
# 1. LOAD
# ----------------------------------------------------------
cat("1. Loading metadata...\n")
meta <- read.csv(file.path(dest, "hmp2_metadata_2018-08-20.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
cat("   Dimensions:", nrow(meta), "rows x", ncol(meta), "columns\n")

# ----------------------------------------------------------
# 2. DIAGNOSIS DISTRIBUTION
# ----------------------------------------------------------
cat("\n2. Diagnosis breakdown:\n")
print(table(meta$diagnosis))

# ----------------------------------------------------------
# 3. DATA TYPE BREAKDOWN
# ----------------------------------------------------------
cat("\n3. Data types available:\n")
print(sort(table(meta$data_type), decreasing = TRUE))

# ----------------------------------------------------------
# 4. SELECT KEY CLINICAL COLUMNS
# ----------------------------------------------------------
cat("\n4. Selecting key columns...\n")

key_cols <- c(
  "External ID",          # sample ID linking to omics files
  "Participant ID",       # patient ID
  "data_type",            # metagenomics / metabolomics / proteomics etc.
  "diagnosis",            # CD / UC / nonIBD
  "week_num",             # time point
  "consent_age",          # patient age
  "sex",                  # sex
  "race",                 # race/ethnicity
  "BMI",                  # BMI
  "site_name",            # collection site
  "fecalcal_ng_ml",       # fecal calprotectin (inflammation marker)
  "hbi",                  # Harvey-Bradshaw Index (CD activity)
  "sccai",                # SCCAI (UC activity)
  "is_inflamed"           # inflamed tissue flag
)

# Keep only columns that actually exist
key_cols_exist <- key_cols[key_cols %in% colnames(meta)]
meta_clean <- meta[, key_cols_exist]

# ----------------------------------------------------------
# 5. CLEAN DIAGNOSIS COLUMN
# ----------------------------------------------------------
cat("\n5. Cleaning diagnosis labels...\n")
meta_clean$diagnosis_binary <- case_when(
  meta_clean$diagnosis == "nonIBD" ~ "nonIBD",
  meta_clean$diagnosis %in% c("CD", "UC") ~ "IBD",
  TRUE ~ NA_character_
)
meta_clean <- meta_clean %>% filter(!is.na(diagnosis_binary))

cat("   Binary label counts:\n")
print(table(meta_clean$diagnosis_binary))

# ----------------------------------------------------------
# 6. SEPARATE BY DATA TYPE
# ----------------------------------------------------------
# Create one clean table per omics type, each linking sample_id to diagnosis
cat("\n6. Splitting by data type...\n")

meta_metagenomics     <- meta_clean %>% filter(data_type == "metagenomics")
meta_metabolomics     <- meta_clean %>% filter(data_type == "metabolomics")
meta_proteomics       <- meta_clean %>% filter(data_type == "proteomics")
meta_transcriptomics  <- meta_clean %>% filter(data_type == "host_transcriptomics")
meta_16S              <- meta_clean %>% filter(data_type %in% c("stool_16S","biopsy_16S"))

cat("   Metagenomics samples:    ", nrow(meta_metagenomics), "\n")
cat("   Metabolomics samples:    ", nrow(meta_metabolomics), "\n")
cat("   Proteomics samples:      ", nrow(meta_proteomics), "\n")
cat("   Transcriptomics samples: ", nrow(meta_transcriptomics), "\n")
cat("   16S samples:             ", nrow(meta_16S), "\n")

# ----------------------------------------------------------
# 7. FIND PATIENTS WITH MULTIPLE OMICS (integration candidates)
# ----------------------------------------------------------
cat("\n7. Finding patients with multiple omics layers...\n")

# Count how many omics types each participant has
omics_per_patient <- meta_clean %>%
  group_by(`Participant ID`) %>%
  summarise(
    n_omics_types = n_distinct(data_type),
    has_metagenomics    = any(data_type == "metagenomics"),
    has_metabolomics    = any(data_type == "metabolomics"),
    has_proteomics      = any(data_type == "proteomics"),
    has_transcriptomics = any(data_type == "host_transcriptomics"),
    diagnosis           = first(diagnosis),
    .groups = "drop"
  ) %>%
  arrange(desc(n_omics_types))

cat("   Patients by number of omics layers:\n")
print(table(omics_per_patient$n_omics_types))

multi_omics_patients <- omics_per_patient %>% filter(n_omics_types >= 3)
cat("   Patients with 3+ omics types:", nrow(multi_omics_patients), "\n")
cat("   Their diagnosis breakdown:\n")
print(table(multi_omics_patients$diagnosis))

# ----------------------------------------------------------
# 8. SAVE
# ----------------------------------------------------------
cat("\n8. Saving outputs...\n")

saveRDS(list(
  meta_all             = meta_clean,
  meta_metagenomics    = meta_metagenomics,
  meta_metabolomics    = meta_metabolomics,
  meta_proteomics      = meta_proteomics,
  meta_transcriptomics = meta_transcriptomics,
  omics_per_patient    = omics_per_patient,
  multi_omics_patients = multi_omics_patients
), file.path(dest, "data/hmp2_metadata_clean.rds"))

write.csv(meta_clean,           file.path(dest, "data/hmp2_metadata_clean.csv"),         row.names = FALSE)
write.csv(omics_per_patient,    file.path(dest, "data/hmp2_omics_per_patient.csv"),      row.names = FALSE)
write.csv(multi_omics_patients, file.path(dest, "data/hmp2_multi_omics_patients.csv"),   row.names = FALSE)

cat("   Saved: data/hmp2_metadata_clean.rds\n")
cat("   Saved: data/hmp2_metadata_clean.csv\n")
cat("   Saved: data/hmp2_omics_per_patient.csv\n")
cat("\n=== DONE ===\n")
