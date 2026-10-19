d193 <- read.csv("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE193677/GSE193677_DEG_results.csv", stringsAsFactors=FALSE)
pxd  <- read.csv("C:/Users/Aarati/Desktop/IBD datasets/processed/proteomics_MaxQuant/proteomics_MaxQuant_DE_results.csv", stringsAsFactors=FALSE)

cat("GSE193677 gene examples:", paste(head(d193$gene, 8), collapse=", "), "\n")
cat("PXD001608 protein examples:", paste(head(pxd$protein, 8), collapse=", "), "\n")
cat("GSE193677 starts with ENSG:", grepl("^ENSG", d193$gene[1]), "\n")

overlap <- intersect(d193$gene, pxd$protein)
cat("Direct overlap:", length(overlap), "\n")

# Check the normalized matrix row names for GSE193677
rds193 <- readRDS("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE193677/GSE193677_normalized.rds")
cat("GSE193677 norm matrix rownames (first 5):", paste(rownames(rds193$normalized_matrix)[1:5], collapse=", "), "\n")
