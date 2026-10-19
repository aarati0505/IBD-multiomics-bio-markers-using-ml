library(data.table)
pg <- fread("C:/Users/Aarati/Desktop/IBD datasets/proteomics/PXD001608_UC_colon_biopsy/proteinGroups.txt",
            sep = "\t", header = TRUE, data.table = FALSE)

cat("Rows before filter:", nrow(pg), "\n")

# Filter
if ("Potential contaminant" %in% colnames(pg)) pg <- pg[pg[["Potential contaminant"]] != "+", ]
if ("Reverse" %in% colnames(pg))               pg <- pg[pg[["Reverse"]] != "+", ]
if ("Only identified by site" %in% colnames(pg)) pg <- pg[pg[["Only identified by site"]] != "+", ]

cat("Rows after filter:", nrow(pg), "\n")

ibaq_cols <- grep("^iBAQ ", colnames(pg), value = TRUE)
cat("iBAQ cols:", length(ibaq_cols), "\n")
cat("First col name:", ibaq_cols[1], "\n")
cat("Class:", class(pg[[ibaq_cols[1]]]), "\n")
cat("First 5 values:", paste(head(pg[[ibaq_cols[1]]], 5), collapse = ", "), "\n")
cat("Gene names rows:", length(pg[["Gene names"]]), "\n")
cat("Rows match:", nrow(pg) == length(pg[["Gene names"]]), "\n")
