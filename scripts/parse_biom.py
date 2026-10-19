"""
HMP2 Metabolomics BIOM Parser
Loads the BIOM file using the official biom-format package,
saves the sparse matrix as a filtered CSV for R.
"""

import biom
import pandas as pd
import numpy as np
import os
import time

BASE     = r"C:\Users\Aarati\Desktop\IBD datasets"
BIOM_PATH = os.path.join(BASE, "metabolomics", "HMP2_metabolomics.biom.txt")
META_PATH = os.path.join(BASE, "data", "hmp2_metadata_clean.csv")
OUT_DIR   = os.path.join(BASE, "processed", "HMP2_metabolomics")
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("HMP2 Metabolomics BIOM Parser")
print(f"File size: {os.path.getsize(BIOM_PATH)/1e6:.1f} MB")
print("=" * 60)

# ----------------------------------------------------------------
# 1. LOAD BIOM FILE
# ----------------------------------------------------------------
print("\nStep 1: Loading BIOM file (may take 1-2 minutes)...")
t0 = time.time()
table = biom.load_table(BIOM_PATH)
print(f"  Done in {time.time()-t0:.1f} sec")
print(f"  Shape: {table.shape[0]} features x {table.shape[1]} samples")

# ----------------------------------------------------------------
# 2. ATTACH METADATA — filter to metabolomics rows
# ----------------------------------------------------------------
print("\nStep 2: Loading and attaching sample metadata...")
meta = pd.read_csv(META_PATH)
mbx_meta = (meta[meta["data_type"] == "metabolomics"]
            .drop_duplicates(subset="External ID")
            .set_index("External ID"))

sample_ids = list(table.ids(axis="sample"))
print(f"  BIOM sample IDs: {len(sample_ids)}")
matched    = mbx_meta.reindex(sample_ids)
has_meta   = matched["diagnosis_binary"].notna()
print(f"  Samples with metadata: {has_meta.sum()} / {len(sample_ids)}")

print("  Diagnosis breakdown:")
print(matched["diagnosis_binary"].value_counts().to_string())

# ----------------------------------------------------------------
# 3. CONVERT TO DENSE DATAFRAME
#    features x samples — only samples with metadata
# ----------------------------------------------------------------
print("\nStep 3: Converting to dense DataFrame...")
t1 = time.time()

valid_samples = [s for s in sample_ids if has_meta[s]]

# Filter table to valid samples only before converting (saves RAM)
table_filt = table.filter(valid_samples, axis="sample", inplace=False)
print(f"  Filtered to {table_filt.shape[1]} samples")

# Convert to dense numpy array (features x samples)
# Use to_dataframe for convenience — returns features as rows
df = table_filt.to_dataframe(dense=True)
print(f"  Dense DataFrame: {df.shape[0]} x {df.shape[1]}")
print(f"  Done in {time.time()-t1:.1f} sec")

# ----------------------------------------------------------------
# 4. FILTER FEATURES BY DETECTION RATE
#    Keep features detected (>0) in >10% of samples in either group
# ----------------------------------------------------------------
print("\nStep 4: Filtering features by detection rate (>10%)...")
print(f"  Before: {df.shape[0]} features")

diag = matched.loc[valid_samples, "diagnosis_binary"]
ibd_s  = diag[diag == "IBD"].index.tolist()
ctrl_s = diag[diag == "nonIBD"].index.tolist()

det_ibd  = (df[ibd_s]  > 0).mean(axis=1)
det_ctrl = (df[ctrl_s] > 0).mean(axis=1)
keep     = (det_ibd >= 0.10) | (det_ctrl >= 0.10)
df       = df.loc[keep]
print(f"  After: {df.shape[0]} features")

# ----------------------------------------------------------------
# 5. SAVE MATRIX AS CSV
# ----------------------------------------------------------------
print("\nStep 5: Saving filtered matrix as CSV...")
t2 = time.time()
matrix_path = os.path.join(OUT_DIR, "HMP2_metabolomics_matrix.csv")
df.to_csv(matrix_path)
size_mb = os.path.getsize(matrix_path) / 1e6
print(f"  Saved: HMP2_metabolomics_matrix.csv ({size_mb:.1f} MB)")
print(f"  Done in {time.time()-t2:.1f} sec")

# ----------------------------------------------------------------
# 6. SAVE SAMPLE METADATA CSV
# ----------------------------------------------------------------
print("\nStep 6: Saving sample metadata...")
col_data = matched.loc[valid_samples, [
    "Participant ID", "diagnosis", "diagnosis_binary",
    "week_num", "sex", "site_name"
]].copy()
col_data.index.name = "sample_id"
col_data.to_csv(os.path.join(OUT_DIR, "HMP2_metabolomics_col_data.csv"))
print("  Saved: HMP2_metabolomics_col_data.csv")

# ----------------------------------------------------------------
# 7. SUMMARY
# ----------------------------------------------------------------
print("\n" + "=" * 60)
print("=== DONE ===")
print(f"Features retained: {df.shape[0]}")
print(f"Samples:           {df.shape[1]}")
print(f"IBD:               {len(ibd_s)}")
print(f"nonIBD:            {len(ctrl_s)}")
print(f"\nOutputs saved to: {OUT_DIR}")
print("=" * 60)
