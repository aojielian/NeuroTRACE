#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
PROJECT <- file.path(BASE, "neurotrace_algorithm_project")
STEP <- file.path(PROJECT, "step06_cross_cohort_validation")
OUT <- file.path(STEP, "results")
FIG <- file.path(STEP, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

MODULE_FILE <- file.path(BASE, "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")
CHOSEN_INPUTS <- file.path(OUT, "32_step06B_chosen_external_inputs.tsv")
OBS_MODELS <- file.path(OUT, "35_step06B_external_NTM_score_models.tsv")

N_PERM <- 1000
MAIN_TOP_N <- c(200, 500)
set.seed(20260508)

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
  common <- intersect(rownames(mat_z), genes)
  if (length(common) < 5) return(rep(NA_real_, ncol(mat_z)))
  w <- weights[match(common, genes)]
  as.numeric(crossprod(w, mat_z[common, , drop = FALSE]) / (sum(abs(w)) + 1e-12))
}

fit_beta <- function(score, md, X) {
  y <- as.numeric(score)
  ok <- is.finite(y) & complete.cases(X)
  y <- y[ok]
  X2 <- X[ok, , drop = FALSE]
  md2 <- md[ok]

  if (length(unique(md2$dx)) < 2) return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))

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

make_random_module <- function(module_genes, module_weights, universe_genes) {
  n <- length(module_genes)
  if (length(universe_genes) < n) stop("Universe smaller than module size.")
  sampled <- sample(universe_genes, n, replace = FALSE)

  # 保留权重分布和符号结构，只随机替换 gene identity
  weights <- sample(module_weights, n, replace = FALSE)
  list(genes = sampled, weights = weights)
}

empirical_p_greater <- function(obs, null) {
  null <- null[is.finite(null)]
  if (!is.finite(obs) || length(null) == 0) return(NA_real_)
  (sum(null >= obs) + 1) / (length(null) + 1)
}

message("[", timestamp(), "] Step06C matched-random specificity started")

for (f in c(MODULE_FILE, CHOSEN_INPUTS, OBS_MODELS)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

chosen <- fread(CHOSEN_INPUTS)
obs <- fread(OBS_MODELS)

all_perm_rows <- list()
summary_rows <- list()

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

  universe_genes <- rownames(mat_z)

  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      msub <- modules[get("program") == prog & top_n == N]
      module_genes <- unique(msub$gene_key)
      module_genes <- module_genes[module_genes %in% universe_genes]

      # 对重复基因保留 rank 最前
      msub2 <- msub[gene_key %in% module_genes]
      setorder(msub2, rank_within_program)
      msub2 <- msub2[!duplicated(gene_key)]

      module_genes <- msub2$gene_key
      module_weights <- as.numeric(msub2$weight)

      obs_score <- score_weighted(mat_z, module_genes, module_weights)
      obs_fit <- fit_beta(obs_score, md2, X)
      obs_beta <- unname(obs_fit["beta"])

      null_beta <- numeric(N_PERM)
      for (b in seq_len(N_PERM)) {
        rm <- make_random_module(module_genes, module_weights, universe_genes)
        sc <- score_weighted(mat_z, rm$genes, rm$weights)
        ft <- fit_beta(sc, md2, X)
        null_beta[b] <- unname(ft["beta"])
      }

      emp_p <- empirical_p_greater(obs_beta, null_beta)
      z <- (obs_beta - mean(null_beta, na.rm = TRUE)) / (sd(null_beta, na.rm = TRUE) + 1e-12)
      pct <- mean(null_beta <= obs_beta, na.rm = TRUE)

      summary_rows[[length(summary_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        n_common_genes = length(module_genes),
        observed_beta = obs_beta,
        observed_p = unname(obs_fit["p"]),
        null_beta_mean = mean(null_beta, na.rm = TRUE),
        null_beta_sd = sd(null_beta, na.rm = TRUE),
        empirical_p_greater = emp_p,
        observed_beta_z_vs_null = z,
        observed_percentile_vs_null = pct,
        n_perm = N_PERM
      )

      all_perm_rows[[length(all_perm_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        perm_id = seq_len(N_PERM),
        null_beta = null_beta,
        observed_beta = obs_beta
      )
    }
  }
}

perm_dt <- rbindlist(all_perm_rows, fill = TRUE)
sum_dt <- rbindlist(summary_rows, fill = TRUE)
sum_dt[, empirical_fdr := p.adjust(empirical_p_greater, method = "BH")]

safe_fwrite(perm_dt, file.path(OUT, "38_step06C_matched_random_null_betas.tsv.gz"))
safe_fwrite(sum_dt, file.path(OUT, "39_step06C_matched_random_specificity_summary.tsv"))

meta <- sum_dt[, .(
  n_cohorts = .N,
  n_empirical_p_lt_0.05 = sum(empirical_p_greater < 0.05, na.rm = TRUE),
  min_empirical_p = min(empirical_p_greater, na.rm = TRUE),
  mean_observed_beta_z_vs_null = mean(observed_beta_z_vs_null, na.rm = TRUE),
  median_observed_percentile = median(observed_percentile_vs_null, na.rm = TRUE),
  all_direction_positive = all(observed_beta > 0, na.rm = TRUE)
), by = .(program, top_n)]

safe_fwrite(meta, file.path(OUT, "40_step06C_matched_random_meta_summary.tsv"))

# Plot
pdf(file.path(FIG, "14_step06C_matched_random_beta_null.pdf"), width = 10, height = 7)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (prog in unique(sum_dt$program)) {
  for (N in MAIN_TOP_N) {
    pp <- perm_dt[program == prog & top_n == N]
    ss <- sum_dt[program == prog & top_n == N]
    if (nrow(pp) == 0) next
    hist(pp$null_beta, breaks = 40, main = paste0(prog, " top", N),
         xlab = "Random-module beta", col = "grey")
    abline(v = ss$observed_beta, col = "red", lwd = 2)
  }
}
dev.off()

summary_lines <- c(
  "# NeuroTRACE Step06C matched-random specificity summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step06C tests whether external-cohort NTM score associations are stronger than matched random gene modules of the same size and weight/sign distribution.",
  "",
  "## Configuration",
  paste0("- Permutations per cohort/module/top_n: ", N_PERM),
  "- One-sided empirical P tests whether observed beta is greater than random-module beta.",
  "",
  "## Specificity summary",
  paste(capture.output(print(sum_dt)), collapse = "\n"),
  "",
  "## Meta summary",
  paste(capture.output(print(meta)), collapse = "\n"),
  "",
  "## Interpretation",
  "A module passes specificity if observed beta is positive and has low empirical P against matched random modules. If external direction is positive but empirical P is not significant, the result should be described as directional external validation rather than specificity-confirmed validation."
)
writeLines(summary_lines, file.path(OUT, "41_step06C_matched_random_specificity_summary.md"))

message("[", timestamp(), "] Step06C done")
message("[", timestamp(), "] Results: ", OUT)
