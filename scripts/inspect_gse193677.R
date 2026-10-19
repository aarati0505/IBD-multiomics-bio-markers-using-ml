d <- read.csv("C:/Users/Aarati/Desktop/IBD datasets/processed/GSE193677/GSE193677_DEG_results.csv",
              stringsAsFactors = FALSE)
cat("Columns:", paste(names(d), collapse=", "), "\n")
cat("Rows:", nrow(d), "\n")
cat("First 3 rows:\n")
print(head(d, 3))

# Check all p-value related columns
for (col in names(d)) {
  vals <- suppressWarnings(as.numeric(d[[col]]))
  n_lt05 <- sum(!is.na(vals) & vals < 0.05)
  if (n_lt05 > 0) cat(col, "< 0.05:", n_lt05, "\n")
}
