#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
PROJECT <- file.path(BASE, "neurotrace_algorithm_project")
STEP <- file.path(PROJECT, "step08_algorithm_strengthening")
OUT <- file.path(STEP, "results")
FIG <- file.path(STEP, "figures")
PARSED <- file.path(PROJECT, "cross_disease_GEO/parsed")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

MODULE_FILE <- file.path(BASE, "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")

DATASETS <- data.table(
  dataset = c("GSE53987", "GSE12649", "GSE21138"),
  expression_path = file.path(PARSED, c(
    "GSE53987/GSE53987.expression_symbol.tsv.gz",
    "GSE12649/GSE12649.expression_symbol.tsv.gz",
    "GSE21138/GSE21138.expression_symbol.tsv.gz"
  )),
  sample_path = file.path(PARSED, c(
    "GSE53987/GSE53987.samples.standardized.tsv",
    "GSE12649/GSE12649.samples.standardized.tsv",
    "GSE21138/GSE21138.samples.standardized.tsv"
  )),
  preferred_region = c("PFC_BA46", "PFC_BA46", "PFC_BA46")
)

MAIN_TOP_N <- c(200, 500)
TARGET_DX <- c("SCZ", "BD", "MDD")

safe_fwrite <- function(x, file) {
  x <- as.data.table(x)
  if (grepl(".gz", file, fixed = TRUE) && endsWith(file, ".gz")) {
    tmp <- sub(".gz$", "", file)
    fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
    system(sprintf("gzip -f %s", shQuote(tmp)))
  } else {
    fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
  }
}

clean_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

read_expr_symbol <- function(path) {
  dt <- fread(path)
  gene_col <- names(dt)[1]
  setnames(dt, gene_col, "feature_id")
  dt[, gene_key := clean_gene(feature_id)]
  dt <- dt[gene_key != "" & !is.na(gene_key)]
  dt <- dt[!duplicated(gene_key)]

  sample_cols <- setdiff(names(dt), c("feature_id", "gene_key"))
  mat <- as.matrix(dt[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- dt$gene_key
  colnames(mat) <- sample_cols
  mat
}

row_zscore <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

make_model <- function(md, target_dx) {
  md <- copy(md)
  md[, target_bin := ifelse(dx_multi == target_dx, 1, 0)]

  covars <- c()

  for (cc in c("age_cov", "rin_cov", "pmi_cov")) {
    if (cc %in% names(md)) {
      vals <- suppressWarnings(as.numeric(md[[cc]]))
      if (sum(is.finite(vals)) >= 5 && sd(vals, na.rm = TRUE) > 0) {
        med <- median(vals, na.rm = TRUE)
        vals[!is.finite(vals)] <- med
        md[, (cc) := as.numeric(scale(vals))]
        covars <- c(covars, cc)
      }
    }
  }

  if ("sex_cov" %in% names(md)) {
    md[is.na(sex_cov) | sex_cov == "", sex_cov := "Unknown"]
    if (length(unique(md$sex_cov)) >= 2) covars <- c(covars, "sex_cov")
  }

  rhs <- c("target_bin", covars)
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  X <- model.matrix(form, data = md)

  qrX <- qr(X)
  keep <- sort(qrX$pivot[seq_len(qrX$rank)])
  X <- X[, keep, drop = FALSE]

  if (!"target_bin" %in% colnames(X)) stop("target_bin missing from model matrix")

  list(md = md, X = X, covars = covars)
}

score_weighted <- function(mat_z, genes, weights) {
  genes <- clean_gene(genes)
  common <- intersect(rownames(mat_z), genes)
  if (length(common) < 5) {
    return(list(score = rep(NA_real_, ncol(mat_z)), n_common = length(common), common = common))
  }
  w <- weights[match(common, genes)]
  score <- as.numeric(crossprod(w, mat_z[common, , drop = FALSE]) / (sum(abs(w)) + 1e-12))
  names(score) <- colnames(mat_z)
  list(score = score, n_common = length(common), common = common)
}

fit_beta <- function(score, md, X) {
  y <- as.numeric(score)
  ok <- is.finite(y) & complete.cases(X)
  y <- y[ok]
  X2 <- X[ok, , drop = FALSE]
  md2 <- md[ok]

  if (length(unique(md2$target_bin)) < 2) {
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, n = length(y), n_case = sum(md2$target_bin == 1), n_control = sum(md2$target_bin == 0)))
  }

  fit <- lm.fit(x = X2, y = y)
  df <- nrow(X2) - ncol(X2)
  rss <- sum(fit$residuals^2)
  sigma2 <- rss / df
  XtX_inv <- solve(crossprod(X2))
  idx <- which(colnames(X2) == "target_bin")

  beta <- fit$coefficients[idx]
  se <- sqrt(sigma2 * XtX_inv[idx, idx])
  tval <- beta / se
  pval <- 2 * pt(-abs(tval), df = df)

  c(beta = as.numeric(beta), se = as.numeric(se), t = as.numeric(tval), p = as.numeric(pval), n = length(y), n_case = sum(md2$target_bin == 1), n_control = sum(md2$target_bin == 0))
}

auc_rank <- function(y, score) {
  y <- as.integer(y)
  score <- as.numeric(score)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)

  pos <- y == 1
  n_pos <- sum(pos)
  n_neg <- sum(!pos)
  rk <- rank(score, ties.method = "average")
  auc <- (sum(rk[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  as.numeric(auc)
}

message("[", timestamp(), "] Step08E cross-disease NTM scoring started")

if (!file.exists(MODULE_FILE)) stop("Missing NTM module file: ", MODULE_FILE)

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

input_rows <- list()
match_rows <- list()
score_rows <- list()
model_rows <- list()

for (i in seq_len(nrow(DATASETS))) {
  ds <- DATASETS$dataset[i]
  expr_path <- DATASETS$expression_path[i]
  sample_path <- DATASETS$sample_path[i]
  preferred_region <- DATASETS$preferred_region[i]

  if (!file.exists(expr_path)) stop("Missing expression file: ", expr_path)
  if (!file.exists(sample_path)) stop("Missing sample file: ", sample_path)

  message("[", timestamp(), "] Loading ", ds)

  mat <- read_expr_symbol(expr_path)
  md <- fread(sample_path)

  # Select preferred region when available
  if ("region_standard" %in% names(md) && any(md$region_standard == preferred_region, na.rm = TRUE)) {
    md <- md[region_standard == preferred_region]
  }

  md <- md[dx_multi %in% c("Control", TARGET_DX)]
  common_samples <- intersect(colnames(mat), md$sample_id)

  mat <- mat[, common_samples, drop = FALSE]
  md <- md[match(common_samples, sample_id)]
  stopifnot(all(md$sample_id == colnames(mat)))

  # log transform if needed
  if (quantile(mat, 0.99, na.rm = TRUE) > 50) {
    mat <- log2(mat + 1)
  }

  finite_rate <- rowMeans(is.finite(mat))
  mat <- mat[finite_rate >= 0.95, , drop = FALSE]
  mat_z <- row_zscore(mat)

  input_rows[[length(input_rows) + 1]] <- data.table(
    dataset = ds,
    expression_path = expr_path,
    sample_path = sample_path,
    preferred_region = preferred_region,
    n_genes = nrow(mat_z),
    n_samples = ncol(mat_z),
    dx_summary = paste(names(table(md$dx_multi)), as.integer(table(md$dx_multi)), sep = ":", collapse = ";")
  )

  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      msub <- modules[get("program") == prog & top_n == N]
      setorder(msub, rank_within_program)
      msub <- msub[!duplicated(gene_key)]
      msub <- msub[gene_key %in% rownames(mat_z)]

      sc <- score_weighted(mat_z, msub$gene_key, as.numeric(msub$weight))

      match_rows[[length(match_rows) + 1]] <- data.table(
        dataset = ds,
        program = prog,
        top_n = N,
        n_module_genes = uniqueN(modules[get("program") == prog & top_n == N, gene_key]),
        n_common_genes = sc$n_common,
        overlap_rate = sc$n_common / uniqueN(modules[get("program") == prog & top_n == N, gene_key])
      )

      score_dt <- data.table(
        dataset = ds,
        sample_id = names(sc$score),
        program = prog,
        top_n = N,
        ntm_score = as.numeric(sc$score)
      )
      score_dt <- merge(score_dt, md[, .(sample_id, dx_multi, region_standard)], by = "sample_id", all.x = TRUE)
      score_rows[[length(score_rows) + 1]] <- score_dt

      for (target_dx in intersect(TARGET_DX, unique(md$dx_multi))) {
        md_sub <- md[dx_multi %in% c("Control", target_dx)]
        if (md_sub[dx_multi == target_dx, .N] < 5 || md_sub[dx_multi == "Control", .N] < 5) next

        model <- make_model(md_sub, target_dx)
        md2 <- model$md
        X <- model$X

        sc_sub <- score_dt[match(md2$sample_id, sample_id)]
        ft <- fit_beta(sc_sub$ntm_score, md2, X)
        auc <- auc_rank(ifelse(md2$dx_multi == target_dx, 1, 0), sc_sub$ntm_score)

        model_rows[[length(model_rows) + 1]] <- data.table(
          dataset = ds,
          target_dx = target_dx,
          program = prog,
          top_n = N,
          n_common_genes = sc$n_common,
          overlap_rate = sc$n_common / uniqueN(modules[get("program") == prog & top_n == N, gene_key]),
          beta_target_vs_Control = unname(ft["beta"]),
          se = unname(ft["se"]),
          t = unname(ft["t"]),
          p_value = unname(ft["p"]),
          auc_target_vs_Control = auc,
          n = unname(ft["n"]),
          n_target = unname(ft["n_case"]),
          n_Control = unname(ft["n_control"]),
          covariates = paste(model$covars, collapse = ";")
        )
      }
    }
  }
}

inputs <- rbindlist(input_rows, fill = TRUE)
matches <- rbindlist(match_rows, fill = TRUE)
scores <- rbindlist(score_rows, fill = TRUE)
models <- rbindlist(model_rows, fill = TRUE)

if (nrow(models) > 0) {
  models[, fdr := p.adjust(p_value, method = "BH")]
  models[, direction_positive := beta_target_vs_Control > 0]
}

summary_by_dx <- models[, .(
  n_tests = .N,
  n_direction_positive = sum(direction_positive, na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  min_fdr = min(fdr, na.rm = TRUE),
  mean_beta = mean(beta_target_vs_Control, na.rm = TRUE),
  median_auc = median(auc_target_vs_Control, na.rm = TRUE),
  median_overlap_rate = median(overlap_rate, na.rm = TRUE)
), by = .(target_dx, program, top_n)]

summary_by_dx[, direction_concordance_rate := n_direction_positive / n_tests]

safe_fwrite(inputs, file.path(OUT, "79_step08E_cross_disease_input_audit.tsv"))
safe_fwrite(matches, file.path(OUT, "80_step08E_cross_disease_NTM_gene_match.tsv"))
safe_fwrite(scores, file.path(OUT, "81_step08E_cross_disease_NTM_scores.tsv.gz"))
safe_fwrite(models, file.path(OUT, "82_step08E_cross_disease_NTM_models.tsv"))
safe_fwrite(summary_by_dx, file.path(OUT, "83_step08E_cross_disease_NTM_summary.tsv"))

# Disease-level summary
disease_summary <- models[, .(
  n_module_tests = .N,
  n_positive = sum(direction_positive, na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  min_fdr = min(fdr, na.rm = TRUE),
  mean_beta = mean(beta_target_vs_Control, na.rm = TRUE),
  median_auc = median(auc_target_vs_Control, na.rm = TRUE)
), by = target_dx]
disease_summary[, positive_rate := n_positive / n_module_tests]
safe_fwrite(disease_summary, file.path(OUT, "84_step08E_cross_disease_disease_level_summary.tsv"))

# Plot
pdf(file.path(FIG, "18_step08E_cross_disease_NTM_scores.pdf"), width = 12, height = 8)
par(mfrow = c(3, 2), mar = c(5, 4, 3, 1))
for (prog in unique(scores$program)) {
  for (N in MAIN_TOP_N) {
    ss <- scores[program == prog & top_n == N]
    if (nrow(ss) == 0) next
    boxplot(ntm_score ~ dx_multi + dataset, data = ss,
            las = 2, main = paste0(prog, " top", N),
            xlab = "", ylab = "NTM score")
  }
}
dev.off()

summary_lines <- c(
  "# NeuroTRACE Step08E cross-disease NTM scoring summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step08E tests whether ASD-derived NeuroTRACE-native modules generalize to SCZ/BD/MDD adult cortex datasets.",
  "",
  "## Inputs",
  paste(capture.output(print(inputs)), collapse = "\n"),
  "",
  "## Gene matching",
  paste(capture.output(print(matches)), collapse = "\n"),
  "",
  "## Cross-disease model results",
  paste(capture.output(print(models)), collapse = "\n"),
  "",
  "## Summary by disease/module",
  paste(capture.output(print(summary_by_dx)), collapse = "\n"),
  "",
  "## Disease-level summary",
  paste(capture.output(print(disease_summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "Positive beta indicates that the ASD-derived NTM score is higher in the target disease than controls. If SCZ/BD/MDD are also positive, the NTM captures broader psychiatric or cortical dysregulation rather than ASD-only specificity."
)

writeLines(summary_lines, file.path(OUT, "85_step08E_cross_disease_NTM_scoring_summary.md"))

message("[", timestamp(), "] Step08E done")
message("[", timestamp(), "] Results: ", OUT)
