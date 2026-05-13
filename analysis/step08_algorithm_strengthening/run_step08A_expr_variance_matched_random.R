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
OBS_MODELS <- file.path(BASE, "neurotrace_algorithm_project/step06_cross_cohort_validation/results/35_step06B_external_NTM_score_models.tsv")

N_PERM <- 1000
MAIN_TOP_N <- c(200, 500)
N_EXPR_BINS <- 5
N_VAR_BINS <- 5
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

make_bins <- function(x, n_bins = 5) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_integer_, length(x))
  if (sum(ok) == 0) return(out)
  rk <- rank(x[ok], ties.method = "average")
  out[ok] <- pmin(n_bins, pmax(1, ceiling(rk / max(rk) * n_bins)))
  out
}

score_weighted <- function(mat_z, genes, weights) {
  genes <- clean_gene(genes)
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

sample_matched_genes <- function(module_genes, module_weights, gene_stats, universe_genes) {
  module_genes <- clean_gene(module_genes)
  universe_genes <- clean_gene(universe_genes)

  mstats <- gene_stats[gene %chin% module_genes]
  # data.table::setorder 不能直接使用 match(gene, module_genes) 表达式；
  # 先显式生成排序列，确保模块基因顺序与权重顺序一致。
  mstats[, module_order := match(gene, module_genes)]
  setorder(mstats, module_order)

  chosen <- character()
  tiers <- character()

  for (i in seq_len(nrow(mstats))) {
    g <- mstats$gene[i]
    mb <- mstats$expr_bin[i]
    vb <- mstats$var_bin[i]

    pool1 <- gene_stats[
      expr_bin == mb &
        var_bin == vb &
        !gene %chin% module_genes &
        !gene %chin% chosen,
      gene
    ]

    pool2 <- gene_stats[
      expr_bin == mb &
        !gene %chin% module_genes &
        !gene %chin% chosen,
      gene
    ]

    pool3 <- gene_stats[
      var_bin == vb &
        !gene %chin% module_genes &
        !gene %chin% chosen,
      gene
    ]

    pool4 <- gene_stats[
      !gene %chin% module_genes &
        !gene %chin% chosen,
      gene
    ]

    if (length(pool1) > 0) {
      chosen_gene <- sample(pool1, 1)
      tier <- "expr_var_bin"
    } else if (length(pool2) > 0) {
      chosen_gene <- sample(pool2, 1)
      tier <- "expr_bin_only"
    } else if (length(pool3) > 0) {
      chosen_gene <- sample(pool3, 1)
      tier <- "var_bin_only"
    } else if (length(pool4) > 0) {
      chosen_gene <- sample(pool4, 1)
      tier <- "any_nonmodule"
    } else {
      chosen_gene <- sample(gene_stats[!gene %chin% module_genes, gene], 1)
      tier <- "replacement_fallback"
    }

    chosen <- c(chosen, chosen_gene)
    tiers <- c(tiers, tier)
  }

  # Preserve module weight/sign distribution, but randomly assign weights to matched replacement genes
  weights <- sample(module_weights, length(chosen), replace = FALSE)

  list(
    genes = chosen,
    weights = weights,
    match_tiers = tiers
  )
}

empirical_p_greater <- function(obs, null) {
  null <- null[is.finite(null)]
  if (!is.finite(obs) || length(null) == 0) return(NA_real_)
  (sum(null >= obs) + 1) / (length(null) + 1)
}

message("[", timestamp(), "] Step08A expression/variance-matched random specificity started")

for (f in c(MODULE_FILE, CHOSEN_INPUTS, OBS_MODELS)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

chosen_inputs <- fread(CHOSEN_INPUTS)
obs_models <- fread(OBS_MODELS)

null_rows <- list()
summary_rows <- list()
tier_rows <- list()
gene_stat_rows <- list()

for (i in seq_len(nrow(chosen_inputs))) {
  coh <- chosen_inputs$cohort[i]
  message("[", timestamp(), "] Loading cohort: ", coh)

  mat <- read_expr(chosen_inputs$expression_path[i])
  md <- read_meta(chosen_inputs$metadata_path[i])

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

  gene_stats <- data.table(
    gene = rownames(mat),
    expr_mean = rowMeans(mat, na.rm = TRUE),
    expr_var = apply(mat, 1, var, na.rm = TRUE)
  )
  gene_stats <- gene_stats[is.finite(expr_mean) & is.finite(expr_var)]
  gene_stats[, expr_bin := make_bins(expr_mean, N_EXPR_BINS)]
  gene_stats[, var_bin := make_bins(expr_var, N_VAR_BINS)]

  gene_stat_rows[[length(gene_stat_rows) + 1]] <- data.table(
    cohort = coh,
    n_genes = nrow(gene_stats),
    n_expr_bins = uniqueN(gene_stats$expr_bin),
    n_var_bins = uniqueN(gene_stats$var_bin),
    expression_path = chosen_inputs$expression_path[i],
    metadata_path = chosen_inputs$metadata_path[i]
  )

  mat <- mat[gene_stats$gene, , drop = FALSE]
  mat_z <- row_zscore(mat)

  model <- make_model(md)
  md2 <- model$md
  X <- model$X
  mat_z <- mat_z[, md2$sample_id, drop = FALSE]

  universe_genes <- rownames(mat_z)

  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      msub <- modules[get("program") == prog & top_n == N]
      setorder(msub, rank_within_program)
      msub <- msub[!duplicated(gene_key)]
      msub <- msub[gene_key %in% universe_genes]

      module_genes <- msub$gene_key
      module_weights <- as.numeric(msub$weight)

      if (length(module_genes) < 20) {
        warning("Low overlap for ", coh, " ", prog, " top", N, ": ", length(module_genes))
      }

      message("[", timestamp(), "] Running matched-random: cohort=", coh,
              " program=", prog, " top_n=", N,
              " n_common_genes=", length(module_genes),
              " N_PERM=", N_PERM)

      obs_score <- score_weighted(mat_z, module_genes, module_weights)
      obs_fit <- fit_beta(obs_score, md2, X)
      obs_beta <- unname(obs_fit["beta"])

      null_beta <- numeric(N_PERM)
      tier_count <- data.table(
        expr_var_bin = integer(N_PERM),
        expr_bin_only = integer(N_PERM),
        var_bin_only = integer(N_PERM),
        any_nonmodule = integer(N_PERM),
        replacement_fallback = integer(N_PERM)
      )

      for (b in seq_len(N_PERM)) {
        if (b == 1 || b %% 100 == 0 || b == N_PERM) {
          message("[", timestamp(), "]   permutation ", b, "/", N_PERM,
                  " cohort=", coh, " program=", prog, " top_n=", N)
          flush.console()
        }

        rm <- sample_matched_genes(module_genes, module_weights, gene_stats, universe_genes)
        sc <- score_weighted(mat_z, rm$genes, rm$weights)
        ft <- fit_beta(sc, md2, X)
        null_beta[b] <- unname(ft["beta"])

        tt <- table(factor(rm$match_tiers, levels = names(tier_count)))

        # data.table 多列赋值在这里容易把多个值误塞入一个列；
        # 改为逐列赋值，保证每个 matching tier 都正确记录。
        for (tier_name in names(tt)) {
          tier_count[b, (tier_name) := as.integer(tt[[tier_name]])]
        }
      }

      emp_p <- empirical_p_greater(obs_beta, null_beta)
      null_mean <- mean(null_beta, na.rm = TRUE)
      null_sd <- sd(null_beta, na.rm = TRUE)
      z <- (obs_beta - null_mean) / (null_sd + 1e-12)
      pct <- mean(null_beta <= obs_beta, na.rm = TRUE)

      summary_rows[[length(summary_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        n_common_genes = length(module_genes),
        observed_beta = obs_beta,
        observed_p = unname(obs_fit["p"]),
        null_beta_mean = null_mean,
        null_beta_sd = null_sd,
        empirical_p_greater = emp_p,
        observed_beta_z_vs_exprvar_null = z,
        observed_percentile_vs_exprvar_null = pct,
        n_perm = N_PERM,
        exact_expr_var_match_rate_mean = mean(tier_count$expr_var_bin / length(module_genes), na.rm = TRUE),
        expr_bin_or_better_match_rate_mean = mean((tier_count$expr_var_bin + tier_count$expr_bin_only) / length(module_genes), na.rm = TRUE)
      )

      null_rows[[length(null_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        perm_id = seq_len(N_PERM),
        null_beta = null_beta,
        observed_beta = obs_beta
      )

      tier_dt <- copy(tier_count)
      tier_dt[, `:=`(
        cohort = coh,
        program = prog,
        top_n = N,
        perm_id = seq_len(N_PERM)
      )]
      tier_rows[[length(tier_rows) + 1]] <- tier_dt
    }
  }
}

gene_stats_summary <- rbindlist(gene_stat_rows, fill = TRUE)
null_dt <- rbindlist(null_rows, fill = TRUE)
summary_dt <- rbindlist(summary_rows, fill = TRUE)
tier_dt <- rbindlist(tier_rows, fill = TRUE)

summary_dt[, empirical_fdr := p.adjust(empirical_p_greater, method = "BH")]

meta_dt <- summary_dt[, .(
  n_cohorts = .N,
  n_empirical_p_lt_0.05 = sum(empirical_p_greater < 0.05, na.rm = TRUE),
  min_empirical_p = min(empirical_p_greater, na.rm = TRUE),
  mean_observed_beta_z_vs_exprvar_null = mean(observed_beta_z_vs_exprvar_null, na.rm = TRUE),
  median_observed_percentile = median(observed_percentile_vs_exprvar_null, na.rm = TRUE),
  mean_exact_expr_var_match_rate = mean(exact_expr_var_match_rate_mean, na.rm = TRUE),
  all_direction_positive = all(observed_beta > 0, na.rm = TRUE)
), by = .(program, top_n)]

safe_fwrite(gene_stats_summary, file.path(OUT, "42_step08A_external_gene_stat_bin_summary.tsv"))
safe_fwrite(null_dt, file.path(OUT, "43_step08A_exprvar_matched_random_null_betas.tsv.gz"))
safe_fwrite(summary_dt, file.path(OUT, "44_step08A_exprvar_matched_random_specificity_summary.tsv"))
safe_fwrite(meta_dt, file.path(OUT, "45_step08A_exprvar_matched_random_meta_summary.tsv"))
safe_fwrite(tier_dt, file.path(OUT, "46_step08A_exprvar_matching_tier_counts.tsv.gz"))

pdf(file.path(FIG, "15_step08A_exprvar_matched_random_nulls.pdf"), width = 10, height = 7)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (prog in unique(summary_dt$program)) {
  for (N in MAIN_TOP_N) {
    pp <- null_dt[program == prog & top_n == N]
    ss <- summary_dt[program == prog & top_n == N]
    if (nrow(pp) == 0) next
    hist(pp$null_beta, breaks = 40,
         main = paste0(prog, " top", N),
         xlab = "Expression/variance-matched random beta",
         col = "grey")
    abline(v = ss$observed_beta, col = "red", lwd = 2)
  }
}
dev.off()

summary_lines <- c(
  "# NeuroTRACE Step08A expression/variance-matched random specificity summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step08A strengthens Step06C by matching random modules not only by module size and weight/sign distribution, but also by external-cohort gene expression abundance and expression variance bins.",
  "",
  "## Configuration",
  paste0("- Permutations per cohort/module/top_n: ", N_PERM),
  paste0("- Expression bins: ", N_EXPR_BINS),
  paste0("- Variance bins: ", N_VAR_BINS),
  "- One-sided empirical P tests whether observed beta is greater than expression/variance-matched random-module beta.",
  "",
  "## Gene statistic bin summary",
  paste(capture.output(print(gene_stats_summary)), collapse = "\n"),
  "",
  "## Specificity summary",
  paste(capture.output(print(summary_dt)), collapse = "\n"),
  "",
  "## Meta summary",
  paste(capture.output(print(meta_dt)), collapse = "\n"),
  "",
  "## Interpretation",
  "A module passes stronger specificity if observed beta remains positive and has low empirical P against expression/variance-matched random modules. This is a stronger control than Step06C and is more appropriate for a methods-oriented manuscript."
)
writeLines(summary_lines, file.path(OUT, "47_step08A_exprvar_matched_random_specificity_summary.md"))

message("[", timestamp(), "] Step08A done")
message("[", timestamp(), "] Results: ", OUT)
