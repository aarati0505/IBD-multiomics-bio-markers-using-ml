import pandas as pd

df = pd.read_excel(
    r"C:\Users\Aarati\Desktop\IBD datasets\metabolomics\C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"
)
print("Shape:", df.shape)
print("Columns:", list(df.columns))
print(df.head(5).to_string())

# Save as CSV for R
df.to_csv(
    r"C:\Users\Aarati\Desktop\IBD datasets\metabolomics\c18n_annotation.csv",
    index=False
)
print("\nSaved: c18n_annotation.csv")
