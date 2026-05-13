#!/usr/bin/env Rscript

gandal_rdata <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"

env <- new.env()
load(gandal_rdata, envir = env)

cat("===== objects =====\n")
print(ls(env))

if ("datMeta_model" %in% ls(env)) {
  meta <- get("datMeta_model", env)
} else {
  stop("datMeta_model not found")
}

cat("===== metadata dim =====\n")
print(dim(meta))

cat("===== metadata columns =====\n")
print(colnames(meta))

cat("===== first few rows =====\n")
print(head(meta[, seq_len(min(20, ncol(meta))), drop = FALSE]))

cat("===== value summary for character/factor-like columns =====\n")
for (cc in colnames(meta)) {
  v <- meta[[cc]]
  if (is.character(v) || is.factor(v) || length(unique(v)) <= 20) {
    cat("\n---", cc, "---\n")
    print(head(sort(table(v, useNA = "ifany"), decreasing = TRUE), 20))
  }
}
