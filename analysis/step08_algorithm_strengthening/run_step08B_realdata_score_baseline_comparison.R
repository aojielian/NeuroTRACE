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
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

MODULE_FILE <- file.path(BASE, "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")
CHOSEN_INPUTS <- file.path(BASE, "neurotrace_algorithm_project/step06_cross_cohort_validation/results/32_step06B_chosen_external_inputs.tsv")

MAIN_TOP_N <- c(200, 500)

SCORING_METHODS <- c(
  "neurotrace_weighted_signed",
  "sign_only",
  "unsigned_mean",
  "abs_weight_unsigned"
)

safe_fwrite <- function(x, file) {
  x <- as.data.table(x)
  if (grepl("\\.gz$", file)) {
    tmp <- sub("\\.gz$", "", file)
    fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
    system(sprintf("gzip -f %s", shQuote(tmp)))
  } else {
    fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
  }
}

clean_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

standardize_dx <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("control|ctl|normal|unaffected", y)] <- "Control"
  out[grepl("^asd$|autism|autistic|case", y)] <- "ASD"
  out
}

detect_col <- function(dt, candidates, regex = NULL) {
  cn <- names(dt)
  low <- tolower(cn)
  for (cc in candidates) {
    hit <- which(low == tolower(cc))
    if (length(hit) > 0) return(cn[hit[1]])
  }
  if (!is.null(regex)) {
    hit <- grep(regex, low)
    if (length(hit) > 0) return(cn[hit[1]])
  }
  NA_character_
}

read_expr <- function(path) {
  dt <- fread(path)
  gene_col <- names(dt)[1]
  setnames(dt, gene_col, "feature_id")
  dt[, gene_key := clean_gene(feature_id)]
  dt <- dt[gene_key != ""]
  dt <- dt[!duplicated(gene_key)]

  sample_cols <- setdiff(names(dt), c("feature_id", "gene_key"))
  mat <- as.matrix(dt[, ..sample_cols])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- dt$gene_key
  colnames(mat) <- sample_cols
  mat
}

read_meta <- function(path) {
  md <- fread(path)
  sample_col <- detect_col(md, c("sample_id", "local_sample_id", "geo_accession", "sample", "id"))
  dx_col <- detect_col(md, c("diagnosis", "disease_status", "disease status", "disease status:ch1", "diagnosis:ch1", "dx", "group"), regex = "diagnosis|disease|dx|group")

  if (is.na(sample_col)) stop("Cannot detect sample_id column in metadata: ", path)
  if (is.na(dx_col)) stop("Cannot detect diagnosis column in metadata: ", path)

  setnames(md, sample_col, "sample_id_raw", skip_absent = TRUE)
  setnames(md, dx_col, "diagnosis_raw", skip_absent = TRUE)

  md[, sample_id := as.character(sample_id_raw)]
  md[, dx := standardize_dx(diagnosis_raw)]

  age_col <- detect_col(md, c("age", "age:ch1"))
  sex_col <- detect_col(md, c("sex", "Sex", "sex:ch1"))
  rin_col <- detect_col(md, c("rin", "RIN", "rin:ch1"))
  pmi_col <- detect_col(md, c("pmi", "PMI", "pmi:ch1", "postmortem.interval:ch1"))

  if (!is.na(age_col)) suppressWarnings(md[, age_cov := as.numeric(get(age_col))])
  if (!is.na(rin_col)) suppressWarnings(md[, rin_cov := as.numeric(get(rin_col))])
  if (!is.na(pmi_col)) suppressWarnings(md[, pmi_cov := as.numeric(get(pmi_col))])
  if (!is.na(sex_col)) md[, sex_cov := as.factor(as.character(get(sex_col)))]

  md
}

make_model <- function(md) {
  md <- copy(md)
  md <- md[dx %in% c("ASD", "Control")]
  md[, dx_bin := ifelse(dx == "ASD", 1, 0)]

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

  rhs <- c("dx_bin", covars)
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  X <- model.matrix(form, data = md)

  qrX <- qr(X)
  keep <- sort(qrX$pivot[seq_len(qrX$rank)])
  X <- X[, keep, drop = FALSE]
  if (!"dx_bin" %in% colnames(X)) stop("dx_bin missing from model matrix")

  list(md = md, X = X, covars = covars)
}

row_zscore <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

score_weighted <- function(mat_z, genes, weights) {
  genes <- clean_gene(genes)
  common <- intersect(rownames(mat_z), genes)
  if (length(common) < 5) return(list(score = rep(NA_real_, ncol(mat_z)), n_common = length(common)))

  w <- weights[match(common, genes)]
  score <- as.numeric(crossprod(w, mat_z[common, , drop = FALSE]) / (sum(abs(w)) + 1e-12))
  names(score) <- colnames(mat_z)

  list(score = score, n_common = length(common))
}

make_method_weights <- function(msub, method) {
  w <- as.numeric(msub$weight)

  if (method == "neurotrace_weighted_signed") {
    return(w)
  }

  if (method == "sign_only") {
    return(sign(w))
  }

  if (method == "unsigned_mean") {
    return(rep(1, length(w)))
  }

  if (method == "abs_weight_unsigned") {
    return(abs(w))
  }

  stop("Unknown scoring method: ", method)
}

fit_beta <- function(score, md, X) {
  y <- as.numeric(score)
  ok <- is.finite(y) & complete.cases(X)
  y <- y[ok]
  X2 <- X[ok, , drop = FALSE]
  md2 <- md[ok]

  if (length(unique(md2$dx)) < 2) {
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
  }

  fit <- lm.fit(x = X2, y = y)
  df <- nrow(X2) - ncol(X2)
  rss <- sum(fit$residuals^2)
  sigma2 <- rss / df
  XtX_inv <- solve(crossprod(X2))
  idx <- which(colnames(X2) == "dx_bin")
  beta <- fit$coefficients[idx]
  se <- sqrt(sigma2 * XtX_inv[idx, idx])
  tval <- beta / se
  pval <- 2 * pt(-abs(tval), df = df)

  c(beta = as.numeric(beta), se = as.numeric(se), t = as.numeric(tval), p = as.numeric(pval))
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

message("[", timestamp(), "] Step08B real-data score baseline comparison started")

for (f in c(MODULE_FILE, CHOSEN_INPUTS)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

chosen <- fread(CHOSEN_INPUTS)

model_rows <- list()
score_rows <- list()
match_rows <- list()

for (i in seq_len(nrow(chosen))) {
  coh <- chosen$cohort[i]
  message("[", timestamp(), "] Loading cohort: ", coh)

  mat <- read_expr(chosen$expression_path[i])
  md <- read_meta(chosen$metadata_path[i])

  common_samples <- intersect(colnames(mat), md$sample_id)
  mat <- mat[, common_samples, drop = FALSE]
  md <- md[match(common_samples, sample_id)]
  stopifnot(all(md$sample_id == colnames(mat)))

  md <- md[dx %in% c("ASD", "Control")]
  mat <- mat[, md$sample_id, drop = FALSE]

  if (quantile(mat, 0.99, na.rm = TRUE) > 50) {
    mat <- log2(mat + 1)
  }

  finite_rate <- rowMeans(is.finite(mat))
  mat <- mat[finite_rate >= 0.95, , drop = FALSE]
  mat_z <- row_zscore(mat)

  model <- make_model(md)
  md2 <- model$md
  X <- model$X
  mat_z <- mat_z[, md2$sample_id, drop = FALSE]

  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      msub <- modules[get("program") == prog & top_n == N]
      setorder(msub, rank_within_program)
      msub <- msub[!duplicated(gene_key)]
      msub <- msub[gene_key %in% rownames(mat_z)]

      for (method in SCORING_METHODS) {
        weights <- make_method_weights(msub, method)
        sc <- score_weighted(mat_z, msub$gene_key, weights)
        ft <- fit_beta(sc$score[md2$sample_id], md2, X)

        auc <- auc_rank(ifelse(md2$dx == "ASD", 1, 0), sc$score[md2$sample_id])

        model_rows[[length(model_rows) + 1]] <- data.table(
          cohort = coh,
          program = prog,
          top_n = N,
          scoring_method = method,
          n_common_genes = sc$n_common,
          beta_ASD_vs_Control = unname(ft["beta"]),
          se = unname(ft["se"]),
          t = unname(ft["t"]),
          p_value = unname(ft["p"]),
          auc_ASD_vs_Control = auc,
          n = nrow(md2),
          n_ASD = sum(md2$dx == "ASD"),
          n_Control = sum(md2$dx == "Control"),
          covariates = paste(model$covars, collapse = ";")
        )

        score_rows[[length(score_rows) + 1]] <- data.table(
          cohort = coh,
          sample_id = md2$sample_id,
          dx = md2$dx,
          program = prog,
          top_n = N,
          scoring_method = method,
          score = as.numeric(sc$score[md2$sample_id])
        )
      }

      match_rows[[length(match_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        n_module_genes = uniqueN(modules[get("program") == prog & top_n == N, gene_key]),
        n_common_genes = nrow(msub),
        overlap_rate = nrow(msub) / uniqueN(modules[get("program") == prog & top_n == N, gene_key])
      )
    }
  }
}

models <- rbindlist(model_rows, fill = TRUE)
scores <- rbindlist(score_rows, fill = TRUE)
matches <- rbindlist(match_rows, fill = TRUE)

models[, fdr := p.adjust(p_value, method = "BH")]
models[, direction_positive := beta_ASD_vs_Control > 0]

# rank methods within cohort/program/top_n by p and beta
models[, rank_by_p := frank(p_value, ties.method = "min"), by = .(cohort, program, top_n)]
models[, rank_by_beta := frank(-beta_ASD_vs_Control, ties.method = "min"), by = .(cohort, program, top_n)]
models[, rank_by_auc := frank(-auc_ASD_vs_Control, ties.method = "min"), by = .(cohort, program, top_n)]

safe_fwrite(models, file.path(OUT, "48_step08B_realdata_baseline_score_models.tsv"))
safe_fwrite(scores, file.path(OUT, "49_step08B_realdata_baseline_sample_scores.tsv.gz"))
safe_fwrite(matches, file.path(OUT, "50_step08B_realdata_baseline_gene_match.tsv"))

# method summary
method_summary <- models[, .(
  n_tests = .N,
  n_direction_positive = sum(direction_positive, na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  median_p = median(p_value, na.rm = TRUE),
  mean_beta = mean(beta_ASD_vs_Control, na.rm = TRUE),
  median_auc = median(auc_ASD_vs_Control, na.rm = TRUE),
  n_rank1_by_p = sum(rank_by_p == 1, na.rm = TRUE),
  n_rank1_by_beta = sum(rank_by_beta == 1, na.rm = TRUE),
  n_rank1_by_auc = sum(rank_by_auc == 1, na.rm = TRUE)
), by = scoring_method]

method_summary[, direction_concordance_rate := n_direction_positive / n_tests]
setorder(method_summary, -n_rank1_by_p, -mean_beta)

safe_fwrite(method_summary, file.path(OUT, "51_step08B_realdata_baseline_method_summary.tsv"))

# pairwise comparison: neurotrace_weighted_signed vs each baseline
nt <- models[scoring_method == "neurotrace_weighted_signed",
             .(cohort, program, top_n,
               nt_beta = beta_ASD_vs_Control,
               nt_p = p_value,
               nt_auc = auc_ASD_vs_Control)]

pair_rows <- list()
for (method in setdiff(SCORING_METHODS, "neurotrace_weighted_signed")) {
  bb <- models[scoring_method == method,
               .(cohort, program, top_n,
                 baseline_beta = beta_ASD_vs_Control,
                 baseline_p = p_value,
                 baseline_auc = auc_ASD_vs_Control)]
  cmp <- merge(nt, bb, by = c("cohort", "program", "top_n"))
  cmp[, baseline_method := method]
  cmp[, delta_beta := nt_beta - baseline_beta]
  cmp[, delta_auc := nt_auc - baseline_auc]
  cmp[, nt_better_p := nt_p < baseline_p]
  cmp[, nt_better_beta := nt_beta > baseline_beta]
  cmp[, nt_better_auc := nt_auc > baseline_auc]
  pair_rows[[length(pair_rows) + 1]] <- cmp
}
pairwise <- rbindlist(pair_rows, fill = TRUE)
safe_fwrite(pairwise, file.path(OUT, "52_step08B_neurotrace_vs_baseline_pairwise.tsv"))

pair_summary <- pairwise[, .(
  n = .N,
  n_nt_better_p = sum(nt_better_p, na.rm = TRUE),
  n_nt_better_beta = sum(nt_better_beta, na.rm = TRUE),
  n_nt_better_auc = sum(nt_better_auc, na.rm = TRUE),
  mean_delta_beta = mean(delta_beta, na.rm = TRUE),
  mean_delta_auc = mean(delta_auc, na.rm = TRUE)
), by = baseline_method]
pair_summary[, frac_nt_better_p := n_nt_better_p / n]
pair_summary[, frac_nt_better_beta := n_nt_better_beta / n]
pair_summary[, frac_nt_better_auc := n_nt_better_auc / n]
safe_fwrite(pair_summary, file.path(OUT, "53_step08B_neurotrace_vs_baseline_pairwise_summary.tsv"))

# plots
pdf(file.path(FIG, "16_step08B_realdata_baseline_effects.pdf"), width = 11, height = 7)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
boxplot(beta_ASD_vs_Control ~ scoring_method, data = models, las = 2,
        main = "Beta by scoring method", ylab = "Beta ASD vs Control")
boxplot(-log10(p_value) ~ scoring_method, data = models, las = 2,
        main = "-log10(P) by scoring method", ylab = "-log10(P)")
boxplot(auc_ASD_vs_Control ~ scoring_method, data = models, las = 2,
        main = "AUC by scoring method", ylab = "AUC")
barplot(method_summary$n_rank1_by_p, names.arg = method_summary$scoring_method,
        las = 2, main = "Rank-1 by P-value count", ylab = "Count")
dev.off()

summary_lines <- c(
  "# NeuroTRACE Step08B real-data score baseline comparison summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step08B compares the NeuroTRACE signed weighted module score against simpler scoring baselines in external adult ASD cohorts.",
  "",
  "## Scoring methods",
  "- neurotrace_weighted_signed: original signed continuous NTM weights",
  "- sign_only: direction-only score using sign(weight)",
  "- unsigned_mean: unweighted module mean without direction",
  "- abs_weight_unsigned: absolute-weight score ignoring direction",
  "",
  "## Method summary",
  paste(capture.output(print(method_summary)), collapse = "\n"),
  "",
  "## NeuroTRACE vs baseline pairwise summary",
  paste(capture.output(print(pair_summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "If neurotrace_weighted_signed ranks higher than unsigned or sign-only baselines, this supports the value of signed continuous NTM weights. If sign_only performs similarly, the main information is directional structure rather than weight magnitude."
)
writeLines(summary_lines, file.path(OUT, "54_step08B_realdata_score_baseline_comparison_summary.md"))

message("[", timestamp(), "] Step08B done")
message("[", timestamp(), "] Results: ", OUT)
