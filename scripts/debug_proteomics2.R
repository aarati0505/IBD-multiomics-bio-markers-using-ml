library(data.table)
pg <- fread("C:/Users/Aarati/Desktop/IBD datasets/proteomics/PXD001608_UC_colon_biopsy/proteinGroups.txt",
            sep = "\t", nrows = 3, data.table = FALSE)

# Gene-related columns
cat("Gene columns:\n")
print(grep("Gene|gene|name|Name", colnames(pg), value = TRUE))

# iBAQ columns and their values
ibaq_cols <- grep("^iBAQ", colnames(pg), value = TRUE)
cat("\nFirst 3 iBAQ col names:", paste(ibaq_cols[1:3], collapse = " | "), "\n")
cat("Row 1 iBAQ values:", paste(pg[1, ibaq_cols[1:3]], collapse = " | "), "\n")
cat("Row 2 iBAQ values:", paste(pg[2, ibaq_cols[1:3]], collapse = " | "), "\n")
cat("Class of iBAQ col:", class(pg[[ibaq_cols[1]]]), "\n")

# Check if there's a header row issue
cat("\nFirst column name:", colnames(pg)[1], "\n")
cat("First value of col 1:", pg[1,1], "\n")
