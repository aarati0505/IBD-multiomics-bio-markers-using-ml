rds <- readRDS("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE137344/GSE137344_normalized.rds")
cat("GSE137344 class:", class(rds), "\n")
cat("names:", paste(names(rds), collapse=", "), "\n")
cat("norm_mat dim:", dim(rds$normalized_matrix), "\n")
cat("col_data cols:", paste(names(rds$col_data), collapse=", "), "\n")
print(table(rds$col_data$diagnosis_binary))

rds2 <- readRDS("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE186507/GSE186507_normalized.rds")
cat("\nGSE186507 class:", class(rds2), "\n")
if (is.list(rds2)) {
  cat("names:", paste(names(rds2), collapse=", "), "\n")
  cat("norm_mat dim:", dim(rds2$normalized_matrix), "\n")
} else {
  cat("dim:", dim(rds2), "\n")
  cat("rownames sample:", paste(rownames(rds2)[1:3], collapse=", "), "\n")
  cat("colnames sample:", paste(colnames(rds2)[1:3], collapse=", "), "\n")
}

rds3 <- readRDS("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE193677/GSE193677_normalized.rds")
cat("\nGSE193677 class:", class(rds3), "\n")
if (is.list(rds3)) {
  cat("names:", paste(names(rds3), collapse=", "), "\n")
  cat("norm_mat dim:", dim(rds3$normalized_matrix), "\n")
} else {
  cat("dim:", dim(rds3), "\n")
  cat("rownames sample:", paste(rownames(rds3)[1:3], collapse=", "), "\n")
}
