import pandas as pd

path = r"C:\Users\Aarati\Desktop\IBD datasets\metabolomics\C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"
df = pd.read_excel(path)
print("Shape:", df.shape)
print("Columns:", list(df.columns))
print("\nFirst 15 rows:")
print(df.head(15).to_string())
print("\nUnique compound names:", df.iloc[:,1].nunique() if df.shape[1] > 1 else "N/A")
