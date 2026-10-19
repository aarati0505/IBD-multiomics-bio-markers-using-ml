pkgs <- c("biomaRt", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "EnsDb.Hsapiens.v86")
for (p in pkgs) {
  path <- system.file(package = p)
  cat(p, ":", ifelse(nchar(path) > 0, "INSTALLED", "NOT INSTALLED"), "\n")
}
