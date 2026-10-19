rds <- readRDS("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE137344/GSE137344_normalized.rds")
cd  <- rds$col_data
cat("col_data columns:", paste(names(cd), collapse=", "), "\n")
cat("Total samples:", nrow(cd), "\n\n")
print(table(cd$diagnosis, cd$diagnosis_binary))

# Check if there is a more detailed diagnosis column
cat("\nDiagnosis values:\n")
print(table(cd$diagnosis))

# Show first 20 samples with their diagnosis
cat("\nFirst 20 samples:\n")
print(head(cd[order(cd$diagnosis_binary), ], 20))
