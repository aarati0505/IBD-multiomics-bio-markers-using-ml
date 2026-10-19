rds <- readRDS("C:/Users/Aarati/Desktop/IBD datasets/processed/MTBLS8713_metabolomics/MTBLS8713_processed.rds")
cat("Names:", paste(names(rds), collapse=", "), "\n")
cat("norm_matrix dim:", dim(rds$normalized_matrix), "\n")
cat("col_data cols:", paste(names(rds$col_data), collapse=", "), "\n")
print(table(rds$col_data$diagnosis_binary))
cat("sig features:", nrow(rds$sig_features), "\n")
cat("DE cols:", paste(names(rds$de_results), collapse=", "), "\n")

# Check sig features up vs down
up   <- sum(rds$sig_features$log2FC > 0)
down <- sum(rds$sig_features$log2FC < 0)
cat("Up in IBD:", up, "| Down in IBD:", down, "\n")
cat("Top 5 up:\n")
print(head(rds$sig_features[rds$sig_features$log2FC > 0,
  c("metabolite","log2FC","padj")], 5))
cat("Top 5 down:\n")
print(head(rds$sig_features[rds$sig_features$log2FC < 0,
  c("metabolite","log2FC","padj")], 5))
