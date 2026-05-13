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

INPUTS <- data.table(
  cohort = c(
    "GSE102741", "GSE102741",
    "GSE64018", "GSE64018"
  ),
  input_role = c(
    "primary_full", "fallback_symbol10682",
    "primary_adjfpkm", "fallback_symbol10682"
  ),
  expression_path = c(
    file.path(BASE, "step03_GSE102741/standardized/GSE102741.expression.tsv.gz"),
    file.path(BASE, "step05_project_to_universe/GSE102741/GSE102741.symbol10682.expression.tsv.gz"),
    file.path(BASE, "step03_GSE64018/standardized/GSE64018.adjfpkm.expression.tsv.gz"),
    file.path(BASE, "step05_project_to_universe/GSE64018/GSE64018.adjfpkm.symbol10682.expression.tsv.gz")
  ),
  metadata_path = c(
    file.path(BASE, "step03_GSE102741/standardized/GSE102741.samples.tsv.gz"),
    file.path(BASE, "step05_project_to_universe/GSE102741/GSE102741.symbol10682.samples.tsv.gz"),
    file.path(BASE, "step03_GSE64018/standardized/GSE64018.adjfpkm.samples.tsv.gz"),
    file.path(BASE, "step05_project_to_universe/GSE64018/GSE64018.adjfpkm.symbol10682.samples.tsv.gz")
  )
)

MAIN_TOP_N <- c(200, 500)

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
  if (ncol(dt) < 3) stop("Expression file has too few columns: ", path)

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

  list(expr_dt = dt, mat = mat, sample_cols = sample_cols)
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

  # covariates
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

score_module <- function(mat_z, module_sub) {
  genes <- clean_gene(module_sub$gene_symbol_fixed)
  weights <- as.numeric(module_sub$weight)
  names(weights) <- genes
  common <- intersect(rownames(mat_z), names(weights))

  if (length(common) < 5) {
    return(list(score = rep(NA_real_, ncol(mat_z)), n_common = length(common), genes_common = common))
  }

  w <- weights[common]
  score <- as.numeric(crossprod(w, mat_z[common, , drop = FALSE]) / (sum(abs(w)) + 1e-12))
  names(score) <- colnames(mat_z)

  list(score = score, n_common = length(common), genes_common = common)
}

fit_score_model <- function(score, md, X) {
  y <- as.numeric(score)
  ok <- is.finite(y) & complete.cases(X)
  y <- y[ok]
  X2 <- X[ok, , drop = FALSE]
  md2 <- md[ok]

  if (length(unique(md2$dx)) < 2) {
    return(list(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, n = length(y), n_asd = sum(md2$dx == "ASD"), n_control = sum(md2$dx == "Control")))
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

  list(
    beta = as.numeric(beta),
    se = as.numeric(se),
    t = as.numeric(tval),
    p = as.numeric(pval),
    n = length(y),
    n_asd = sum(md2$dx == "ASD"),
    n_control = sum(md2$dx == "Control")
  )
}

message("[", timestamp(), "] Step06B external NTM module validation started")

if (!file.exists(MODULE_FILE)) stop("Missing module file: ", MODULE_FILE)
modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

safe_fwrite(INPUTS, file.path(OUT, "30_step06B_external_input_freeze_candidates.tsv"))

# Choose best input per cohort by module gene overlap
choice_rows <- list()
prepared <- list()

for (coh in unique(INPUTS$cohort)) {
  sub <- INPUTS[cohort == coh]
  best <- NULL

  for (i in seq_len(nrow(sub))) {
    row <- sub[i]
    if (!file.exists(row$expression_path) || !file.exists(row$metadata_path)) next

    expr <- read_expr(row$expression_path)
    md <- read_meta(row$metadata_path)

    sample_common <- intersect(colnames(expr$mat), md$sample_id)
    module_genes <- unique(modules$gene_key)
    gene_common <- intersect(rownames(expr$mat), module_genes)

    candidate <- list(
      cohort = coh,
      input_role = row$input_role,
      expression_path = row$expression_path,
      metadata_path = row$metadata_path,
      n_expr_genes = nrow(expr$mat),
      n_expr_samples = ncol(expr$mat),
      n_meta_rows = nrow(md),
      n_sample_common = length(sample_common),
      n_module_gene_common = length(gene_common),
      expr = expr,
      md = md
    )

    choice_rows[[length(choice_rows) + 1]] <- data.table(
      cohort = coh,
      input_role = row$input_role,
      expression_path = row$expression_path,
      metadata_path = row$metadata_path,
      n_expr_genes = nrow(expr$mat),
      n_expr_samples = ncol(expr$mat),
      n_meta_rows = nrow(md),
      n_sample_common = length(sample_common),
      n_module_gene_common = length(gene_common)
    )

    if (is.null(best) ||
        candidate$n_module_gene_common > best$n_module_gene_common ||
        (candidate$n_module_gene_common == best$n_module_gene_common && candidate$n_sample_common > best$n_sample_common)) {
      best <- candidate
    }
  }

  if (is.null(best)) stop("No usable input found for ", coh)
  prepared[[coh]] <- best
}

choice_dt <- rbindlist(choice_rows, fill = TRUE)
safe_fwrite(choice_dt, file.path(OUT, "31_step06B_input_choice_audit.tsv"))

chosen <- rbindlist(lapply(names(prepared), function(coh) {
  x <- prepared[[coh]]
  data.table(
    cohort = coh,
    chosen_input_role = x$input_role,
    expression_path = x$expression_path,
    metadata_path = x$metadata_path,
    n_expr_genes = x$n_expr_genes,
    n_expr_samples = x$n_expr_samples,
    n_meta_rows = x$n_meta_rows,
    n_sample_common = x$n_sample_common,
    n_module_gene_common = x$n_module_gene_common
  )
}))
safe_fwrite(chosen, file.path(OUT, "32_step06B_chosen_external_inputs.tsv"))

score_rows <- list()
match_rows <- list()
model_rows <- list()

for (coh in names(prepared)) {
  message("[", timestamp(), "] Scoring cohort: ", coh)

  obj <- prepared[[coh]]
  mat <- obj$expr$mat
  md <- copy(obj$md)

  common_samples <- intersect(colnames(mat), md$sample_id)
  mat <- mat[, common_samples, drop = FALSE]
  md <- md[match(common_samples, sample_id)]
  stopifnot(all(md$sample_id == colnames(mat)))

  md <- md[dx %in% c("ASD", "Control")]
  mat <- mat[, md$sample_id, drop = FALSE]

  # log transform if count-like
  if (quantile(mat, 0.99, na.rm = TRUE) > 50) {
    mat <- log2(mat + 1)
  }

  finite_rate <- rowMeans(is.finite(mat))
  mat <- mat[finite_rate >= 0.95, , drop = FALSE]
  mat_z <- row_zscore(mat)

  model <- make_model(md)
  md2 <- model$md
  X <- model$X

  # ensure same order after model filtering
  mat_z <- mat_z[, md2$sample_id, drop = FALSE]

  for (prog in unique(modules$program)) {
    for (N in sort(unique(modules[get("program") == prog]$top_n))) {
      msub <- modules[get("program") == prog & top_n == N]
      sc <- score_module(mat_z, msub)

      match_rows[[length(match_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        n_module_genes = uniqueN(msub$gene_key),
        n_common_genes = sc$n_common,
        overlap_rate = sc$n_common / uniqueN(msub$gene_key),
        common_genes = paste(sc$genes_common, collapse = ";")
      )

      score_dt <- data.table(
        cohort = coh,
        sample_id = names(sc$score),
        program = prog,
        top_n = N,
        ntm_score = as.numeric(sc$score)
      )
      score_dt <- merge(score_dt, md2[, .(sample_id, dx)], by = "sample_id", all.x = TRUE)
      score_rows[[length(score_rows) + 1]] <- score_dt

      fit <- fit_score_model(sc$score[md2$sample_id], md2, X)

      model_rows[[length(model_rows) + 1]] <- data.table(
        cohort = coh,
        chosen_input_role = obj$input_role,
        program = prog,
        top_n = N,
        n_common_genes = sc$n_common,
        overlap_rate = sc$n_common / uniqueN(msub$gene_key),
        beta_ASD_vs_Control = fit$beta,
        se = fit$se,
        t = fit$t,
        p_value = fit$p,
        n = fit$n,
        n_ASD = fit$n_asd,
        n_Control = fit$n_control,
        covariates = paste(model$covars, collapse = ";")
      )
    }
  }
}

matches <- rbindlist(match_rows, fill = TRUE)
scores <- rbindlist(score_rows, fill = TRUE)
models <- rbindlist(model_rows, fill = TRUE)
models[, fdr := p.adjust(p_value, method = "BH")]
models[, direction_match_expected := beta_ASD_vs_Control > 0]

safe_fwrite(matches, file.path(OUT, "33_step06B_NTM_gene_match_by_external_cohort.tsv"))
safe_fwrite(scores, file.path(OUT, "34_step06B_external_sample_NTM_scores.tsv.gz"))
safe_fwrite(models, file.path(OUT, "35_step06B_external_NTM_score_models.tsv"))

# meta-style summaries
summary <- models[, .(
  n_cohorts = .N,
  n_direction_positive = sum(beta_ASD_vs_Control > 0, na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  min_fdr = min(fdr, na.rm = TRUE),
  mean_beta = mean(beta_ASD_vs_Control, na.rm = TRUE),
  median_overlap_rate = median(overlap_rate, na.rm = TRUE)
), by = .(program, top_n)]
summary[, direction_concordance_rate := n_direction_positive / n_cohorts]
safe_fwrite(summary, file.path(OUT, "36_step06B_external_validation_summary.tsv"))

# plots
pdf(file.path(FIG, "13_step06B_external_NTM_scores_boxplot.pdf"), width = 10, height = 6)
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))
for (prog in unique(scores$program)) {
  for (N in MAIN_TOP_N) {
    ss <- scores[program == prog & top_n == N]
    if (nrow(ss) == 0) next
    boxplot(ntm_score ~ dx + cohort, data = ss,
            las = 2, main = paste0(prog, " top", N),
            ylab = "NTM score", xlab = "")
  }
}
dev.off()

summary_lines <- c(
  "# NeuroTRACE Step06B external adult NTM module validation summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step06B validates Gandal-derived NeuroTRACE-native modules in external adult ASD bulk cohorts GSE102741 and GSE64018. It does not use DevMap/DevBridge score tables.",
  "",
  "## Chosen inputs",
  paste(capture.output(print(chosen)), collapse = "\n"),
  "",
  "## Model outputs",
  paste(capture.output(print(models)), collapse = "\n"),
  "",
  "## Cross-cohort summary",
  paste(capture.output(print(summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "Positive beta indicates that the Gandal-derived signed NTM score is higher in ASD in the external cohort. Directional concordance across GSE102741 and GSE64018 is more important than single-cohort significance at this stage, because the external cohorts are small."
)
writeLines(summary_lines, file.path(OUT, "37_step06B_external_NTM_validation_summary.md"))

message("[", timestamp(), "] Step06B done")
message("[", timestamp(), "] Results: ", OUT)
