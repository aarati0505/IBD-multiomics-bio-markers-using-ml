"""
Build an Ensembl ID -> Gene Symbol mapping for GSE193677
by downloading from Ensembl REST API for all ENSG IDs in the dataset
"""
import pandas as pd
import requests
import json
import time
import os

BASE = r"C:\Users\Aarati\Desktop\IBD datasets"

# Load GSE193677 DEG results to get all ENSG IDs
deg193 = pd.read_csv(os.path.join(BASE, "processed/GSE193677/GSE193677_DEG_results.csv"))
ensg_ids = deg193["gene"].dropna().tolist()
print(f"Total ENSG IDs to map: {len(ensg_ids)}")

# Ensembl REST API — batch lookup, 1000 at a time
def get_gene_symbols(ensg_list):
    """Query Ensembl REST API for gene symbols"""
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    results = {}
    
    # Process in batches of 1000
    batch_size = 1000
    for i in range(0, len(ensg_list), batch_size):
        batch = ensg_list[i:i+batch_size]
        payload = json.dumps({"ids": batch})
        
        try:
            r = requests.post(url, headers=headers, data=payload, timeout=60)
            if r.status_code == 200:
                data = r.json()
                for ensg, info in data.items():
                    if info and isinstance(info, dict):
                        symbol = info.get("display_name", None)
                        results[ensg] = symbol
                print(f"  Batch {i//batch_size + 1}: mapped {len(results)} so far")
            else:
                print(f"  Batch {i//batch_size + 1}: HTTP {r.status_code}")
        except Exception as e:
            print(f"  Batch {i//batch_size + 1}: Error - {e}")
        
        time.sleep(0.5)  # be polite to the API
    
    return results

print("Querying Ensembl REST API...")
symbol_map = get_gene_symbols(ensg_ids)
print(f"Successfully mapped: {sum(v is not None for v in symbol_map.values())} / {len(ensg_ids)}")

# Save as CSV
map_df = pd.DataFrame([
    {"ensembl_id": k, "gene_symbol": v}
    for k, v in symbol_map.items()
    if v is not None
])
out_path = os.path.join(BASE, "processed/GSE193677/GSE193677_ensembl_to_symbol.csv")
map_df.to_csv(out_path, index=False)
print(f"Saved: {out_path}")
print(f"Rows: {len(map_df)}")
print("First 10:")
print(map_df.head(10).to_string())
