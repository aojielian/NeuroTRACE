#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
}

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

safe_dir_create <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

safe_write_tsv <- function(df, path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "wt")
    on.exit(close(con), add = TRUE)
    write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  } else {
    write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
}

find_col <- function(df, candidates) {
  cn <- colnames(df)
  low <- tolower(cn)
  for (cc in candidates) {
    hit <- which(low == tolower(cc))
    if (length(hit) > 0) return(cn[hit[1]])
  }
  return(NA_character_)
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- sub("^gene:", "", x)
  x <- sub("^GENE:", "", x)
  trimws(x)
}

harmonize_gene_symbols <- function(x, label = "gene_ids") {
  x0 <- as.character(x)
  x1 <- clean_gene(x0)

  ens_core <- rep(NA_character_, length(x1))
  has_ens <- grepl("ENSG[0-9]+", x1)

  if (any(has_ens)) {
    m <- regexpr("ENSG[0-9]+", x1[has_ens])
    ens_core[has_ens] <- regmatches(x1[has_ens], m)
  }

  if (any(!is.na(ens_core))) {
    if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)) {

      valid_keys <- tryCatch(
        AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL"),
        error = function(e) character(0)
      )

      query_keys <- unique(ens_core[!is.na(ens_core)])
      valid_query <- intersect(query_keys, valid_keys)

      if (length(valid_query) > 0) {
        mapped <- tryCatch(
          AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            keys = valid_query,
            keytype = "ENSEMBL",
            column = "SYMBOL",
            multiVals = "first"
          ),
          error = function(e) NULL
        )

        if (!is.null(mapped)) {
          pos <- which(!is.na(ens_core) & ens_core %in% names(mapped))
          vals <- as.character(mapped[ens_core[pos]])
          ok <- !is.na(vals) & vals != ""

          if (any(ok)) {
            x1[pos[ok]] <- vals[ok]
          }

          msg(label, ": mapped", sum(ok), "/", sum(!is.na(ens_core)),
              "Ensembl-like IDs to gene symbols")
        }
      } else {
        msg("WARNING:", label, "has Ensembl-like IDs but no valid ENSEMBL keys in current org.Hs.eg.db")
      }
    } else {
      msg("WARNING:", label, "has Ensembl-like IDs but AnnotationDbi/org.Hs.eg.db unavailable")
    }
  }

  clean_gene(x1)
}

collapse_duplicate_gene_columns <- function(expr) {
  genes <- colnames(expr)
  keep <- !is.na(genes) & genes != "" & genes != "NA" & genes != "nan"
  expr <- expr[, keep, drop = FALSE]
  genes <- genes[keep]

  if (!any(duplicated(genes))) return(expr)

  msg("Collapsing duplicated gene symbols:", sum(duplicated(genes)), "duplicated columns")
  split_idx <- split(seq_along(genes), genes)

  out <- matrix(NA_real_, nrow = nrow(expr), ncol = length(split_idx))
  rownames(out) <- rownames(expr)
  colnames(out) <- names(split_idx)

  for (j in seq_along(split_idx)) {
    ii <- split_idx[[j]]
    if (length(ii) == 1) {
      out[, j] <- expr[, ii]
    } else {
      out[, j] <- rowMeans(expr[, ii, drop = FALSE], na.rm = TRUE)
    }
  }

  out
}

parse_dx_from_DxReg <- function(x) {
  raw <- as.character(x)
  low <- tolower(raw)

  dx <- rep(NA_character_, length(raw))
  dx[grepl("^asd", low) | grepl("asd_", low)] <- "ASD"
  dx[grepl("^ctl", low) | grepl("ctl_", low) | grepl("^control", low)] <- "Control"

  dx
}

aggregate_subject_level <- function(expr, meta, subject_col, dx_col) {
  if (!subject_col %in% colnames(meta)) {
    stop("Subject column not found: ", subject_col)
  }
  if (!dx_col %in% colnames(meta)) {
    stop("Diagnosis column not found: ", dx_col)
  }

  meta$dx_bootstrap_tmp <- parse_dx_from_DxReg(meta[[dx_col]])

  keep <- meta$dx_bootstrap_tmp %in% c("ASD", "Control")
  meta <- meta[keep, , drop = FALSE]
  expr <- expr[keep, , drop = FALSE]

  sid <- as.character(meta[[subject_col]])
  sid[is.na(sid) | sid == ""] <- rownames(meta)[is.na(sid) | sid == ""]

  groups <- unique(sid)

  msg("Rows before subject aggregation:", nrow(expr))
  msg("Unique subjects:", length(groups))

  expr_agg <- matrix(NA_real_, nrow = length(groups), ncol = ncol(expr))
  rownames(expr_agg) <- groups
  colnames(expr_agg) <- colnames(expr)

  meta_agg <- data.frame(row.names = groups)

  for (g in groups) {
    ii <- which(sid == g)

    expr_agg[g, ] <- colMeans(expr[ii, , drop = FALSE], na.rm = TRUE)

    subm <- meta[ii, , drop = FALSE]

    for (cc in colnames(meta)) {
      v <- subm[[cc]]

      if (is.numeric(v) || is.integer(v)) {
        meta_agg[g, cc] <- mean(v, na.rm = TRUE)
      } else {
        vv <- as.character(v)
        vv <- vv[!is.na(vv) & vv != ""]
        meta_agg[g, cc] <- ifelse(length(vv) == 0, NA, vv[1])
      }
    }

    dx_vals <- unique(meta$dx_bootstrap_tmp[ii])
    dx_vals <- dx_vals[!is.na(dx_vals)]
    meta_agg[g, "dx_bootstrap_tmp"] <- ifelse(length(dx_vals) == 0, NA, dx_vals[1])
  }

  meta_agg$dx_binary <- ifelse(meta_agg$dx_bootstrap_tmp == "ASD", 1, 0)

  list(expr = expr_agg, meta = meta_agg)
}

load_gandal_subject_data <- function(gandal_rdata, dx_col = "DxReg", subject_col = "Subject") {
  env <- new.env()
  load(gandal_rdata, envir = env)

  obj_names <- ls(env)
  msg("Objects in RData:", paste(obj_names, collapse = ", "))

  if (!"datExpr" %in% obj_names) stop("datExpr not found in Gandal RData")
  if (!"datMeta_model" %in% obj_names) stop("datMeta_model not found in Gandal RData")

  expr <- as.matrix(get("datExpr", env))
  mode(expr) <- "numeric"

  meta <- as.data.frame(get("datMeta_model", env))

  msg("Raw datExpr dim:", paste(dim(expr), collapse = " x "))
  msg("Raw datMeta_model dim:", paste(dim(meta), collapse = " x "))

  # In this Gandal object, datExpr rows correspond to datMeta_model rows.
  if (nrow(expr) != nrow(meta)) {
    if (ncol(expr) == nrow(meta)) {
      expr <- t(expr)
    } else {
      stop("Cannot align datExpr and datMeta_model by row count.")
    }
  }

  rownames(expr) <- paste0("sample_", seq_len(nrow(expr)))
  rownames(meta) <- rownames(expr)

  colnames(expr) <- harmonize_gene_symbols(colnames(expr), label = "Gandal datExpr columns")
  expr <- collapse_duplicate_gene_columns(expr)

  dat <- aggregate_subject_level(expr, meta, subject_col = subject_col, dx_col = dx_col)

  expr_subj <- dat$expr
  meta_subj <- dat$meta

  msg("Subject-level expression dim:", paste(dim(expr_subj), collapse = " x "))
  msg("Subject-level metadata dim:", paste(dim(meta_subj), collapse = " x "))
  msg("Subject-level diagnosis table:")
  print(table(meta_subj$dx_bootstrap_tmp, useNA = "ifany"))

  list(expr = expr_subj, meta = meta_subj)
}

select_covariates <- function(meta, core_only = TRUE) {
  covs <- c()

  core_list <- list(
    c("Age", "age"),
    c("Sex", "sex", "gender"),
    c("PMI", "pmi"),
    c("RIN", "rin")
  )

  extra_list <- list(
    c("SeqBatch", "batch", "Batch"),
    c("Ancestry", "ancestry")
  )

  cov_list <- core_list
  if (!core_only) cov_list <- c(core_list, extra_list)

  for (cand in cov_list) {
    cc <- find_col(meta, cand)
    if (!is.na(cc)) covs <- c(covs, cc)
  }

  covs <- unique(covs)

  keep_covs <- c()
  for (cc in covs) {
    v <- meta[[cc]]
    miss <- mean(is.na(v))
    if (miss > 0.3) next
    if (length(unique(v[!is.na(v)])) <= 1) next
    keep_covs <- c(keep_covs, cc)
  }

  keep_covs
}

build_design <- function(meta, covariates) {
  dat <- meta
  dat$dx_binary <- as.numeric(dat$dx_binary)

  rhs <- c("dx_binary")

  for (cc in covariates) {
    if (is.numeric(dat[[cc]]) || is.integer(dat[[cc]])) {
      rhs <- c(rhs, cc)
    } else {
      dat[[cc]] <- as.factor(dat[[cc]])
      rhs <- c(rhs, cc)
    }
  }

  f <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  X <- model.matrix(f, data = dat)

  keep <- apply(X, 2, function(z) sd(z) > 0)
  keep[1] <- TRUE
  X <- X[, keep, drop = FALSE]

  if (!"dx_binary" %in% colnames(X)) stop("dx_binary missing from design matrix")

  list(X = X, formula = paste(deparse(f), collapse = " "))
}

fit_fast <- function(Y, meta, covariates) {
  # Y: samples x genes
  design <- build_design(meta, covariates)
  X <- design$X

  ok <- complete.cases(X) & is.finite(meta$dx_binary)
  X <- X[ok, , drop = FALSE]
  Y <- Y[ok, , drop = FALSE]

  n <- nrow(X)
  p <- ncol(X)

  if (n <= p + 2) stop("Too few samples for model design")

  qrX <- qr(X)
  coef <- qr.coef(qrX, Y)
  fitted <- X %*% coef
  resid <- Y - fitted

  df <- n - qrX$rank
  rss <- colSums(resid^2, na.rm = TRUE)
  sigma2 <- rss / df

  XtX_inv <- tryCatch(
    chol2inv(qr.R(qrX)),
    error = function(e) solve(crossprod(X))
  )

  dx_idx <- which(colnames(X) == "dx_binary")[1]

  beta <- coef[dx_idx, ]
  se <- sqrt(pmax(sigma2, 0) * XtX_inv[dx_idx, dx_idx])
  tval <- beta / se
  pval <- 2 * pt(-abs(tval), df = df)

  signed_score <- sign(beta) * (-log10(pmax(pval, 1e-300)))

  data.frame(
    gene = colnames(Y),
    beta = as.numeric(beta),
    se = as.numeric(se),
    t = as.numeric(tval),
    p_value = as.numeric(pval),
    signed_score = as.numeric(signed_score),
    stringsAsFactors = FALSE
  )
}

derive_modules <- function(stats, thresholds) {
  stats <- stats[is.finite(stats$beta) &
                   is.finite(stats$p_value) &
                   is.finite(stats$signed_score), , drop = FALSE]

  up_order <- stats[order(stats$signed_score, decreasing = TRUE), ]
  down_order <- stats[order(stats$signed_score, decreasing = FALSE), ]
  signed_order <- stats[order(abs(stats$signed_score), decreasing = TRUE), ]

  out <- list()

  for (n in thresholds) {
    out[[paste0("NTM1_ASD_up|top", n)]] <- up_order$gene[seq_len(min(n, nrow(up_order)))]
    out[[paste0("NTM2_ASD_down|top", n)]] <- down_order$gene[seq_len(min(n, nrow(down_order)))]
    out[[paste0("NTM3_ASD_signed|top", n)]] <- signed_order$gene[seq_len(min(n, nrow(signed_order)))]
  }

  out
}

load_frozen_modules <- function(module_weight_file, thresholds) {
  mw <- read.table(module_weight_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")

  gene_col <- find_col(mw, c("gene_symbol_fixed", "gene_symbol", "gene", "symbol"))
  program_col <- find_col(mw, c("program", "module_family", "module"))
  topn_col <- find_col(mw, c("top_n", "topn"))

  if (any(is.na(c(gene_col, program_col, topn_col)))) {
    stop("Cannot detect required columns in module weight file.")
  }

  mw$gene_clean <- clean_gene(mw[[gene_col]])
  mw$program_clean <- as.character(mw[[program_col]])
  mw$top_n_clean <- as.numeric(mw[[topn_col]])

  modules <- list()

  for (prog in c("NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed")) {
    for (n in thresholds) {
      sub <- mw[mw$program_clean == prog & mw$top_n_clean == n, , drop = FALSE]
      modules[[paste0(prog, "|top", n)]] <- unique(sub$gene_clean)
    }
  }

  modules
}

compare_modules <- function(boot_modules, frozen_modules, full_stats, boot_stats, iter) {
  records <- list()
  idx <- 1

  all_common <- intersect(full_stats$gene, boot_stats$gene)

  rho <- suppressWarnings(cor(
    full_stats$signed_score[match(all_common, full_stats$gene)],
    boot_stats$signed_score[match(all_common, boot_stats$gene)],
    method = "spearman",
    use = "complete.obs"
  ))

  for (key in names(frozen_modules)) {
    frozen <- unique(frozen_modules[[key]])
    boot <- unique(boot_modules[[key]])

    inter <- intersect(frozen, boot)
    union <- union(frozen, boot)

    program <- sub("\\|top.*$", "", key)
    top_n <- as.integer(sub("^.*\\|top", "", key))

    if (length(inter) > 0) {
      fb <- full_stats$beta[match(inter, full_stats$gene)]
      bb <- boot_stats$beta[match(inter, boot_stats$gene)]
      sign_conc <- mean(sign(fb) == sign(bb), na.rm = TRUE)
    } else {
      sign_conc <- NA_real_
    }

    records[[idx]] <- data.frame(
      bootstrap_id = iter,
      module = program,
      top_n = top_n,
      frozen_n = length(frozen),
      bootstrap_n = length(boot),
      overlap_n = length(inter),
      union_n = length(union),
      jaccard = ifelse(length(union) > 0, length(inter) / length(union), NA_real_),
      frozen_recall = ifelse(length(frozen) > 0, length(inter) / length(frozen), NA_real_),
      bootstrap_precision = ifelse(length(boot) > 0, length(inter) / length(boot), NA_real_),
      sign_concordance_overlap = sign_conc,
      full_vs_boot_signed_score_spearman = rho,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }

  do.call(rbind, records)
}

main <- function() {
  gandal_rdata <- get_arg("--gandal_rdata", "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData")
  module_weight_file <- get_arg("--module_weights", "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")
  outdir <- get_arg("--outdir", "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12A_v4_full")
  dx_col <- get_arg("--dx_col", "DxReg")
  subject_col <- get_arg("--subject_col", "Subject")
  core_covariates_only <- toupper(get_arg("--core_covariates_only", "TRUE")) %in% c("TRUE", "T", "1", "YES")
  n_boot <- as.integer(get_arg("--n_boot", "300"))
  seed <- as.integer(get_arg("--seed", "20260509"))

  safe_dir_create(outdir)

  thresholds <- c(50, 100, 200, 500)

  msg("Step12A v4 started")
  msg("Gandal RData:", gandal_rdata)
  msg("Module weights:", module_weight_file)
  msg("Diagnosis column:", dx_col)
  msg("Subject column:", subject_col)
  msg("n_boot:", n_boot)

  dat <- load_gandal_subject_data(gandal_rdata, dx_col = dx_col, subject_col = subject_col)

  expr <- dat$expr
  meta <- dat$meta

  covariates <- select_covariates(meta, core_only = core_covariates_only)
  msg("Covariates:", ifelse(length(covariates) == 0, "none", paste(covariates, collapse = ", ")))

  gene_missing <- colMeans(is.na(expr))
  gene_var <- apply(expr, 2, var, na.rm = TRUE)
  keep_gene <- gene_missing < 0.2 & is.finite(gene_var) & gene_var > 0
  expr <- expr[, keep_gene, drop = FALSE]

  msg("Expression after gene filtering:", paste(dim(expr), collapse = " x "))

  full_stats <- fit_fast(expr, meta, covariates)
  full_stats_out <- file.path(outdir, "00_step12A_full_discovery_model_stats.tsv.gz")
  safe_write_tsv(full_stats, full_stats_out)

  frozen_modules <- load_frozen_modules(module_weight_file, thresholds)
  full_modules <- derive_modules(full_stats, thresholds)

  full_compare <- compare_modules(full_modules, frozen_modules, full_stats, full_stats, iter = 0)
  full_compare_out <- file.path(outdir, "00b_step12A_full_model_vs_frozen_modules.tsv")
  safe_write_tsv(full_compare, full_compare_out)

  set.seed(seed)

  asd_ids <- rownames(meta)[meta$dx_bootstrap_tmp == "ASD"]
  ctl_ids <- rownames(meta)[meta$dx_bootstrap_tmp == "Control"]

  if (length(asd_ids) < 5 || length(ctl_ids) < 5) {
    stop("Too few ASD or Control subjects after aggregation.")
  }

  all_records <- list()

  freq_counter <- list()
  for (key in names(frozen_modules)) {
    freq_counter[[key]] <- setNames(rep(0L, length(frozen_modules[[key]])), frozen_modules[[key]])
  }

  for (b in seq_len(n_boot)) {
    if (b %% 25 == 0) msg("Bootstrap", b, "/", n_boot)

    boot_ids <- c(
      sample(asd_ids, length(asd_ids), replace = TRUE),
      sample(ctl_ids, length(ctl_ids), replace = TRUE)
    )

    expr_b <- expr[boot_ids, , drop = FALSE]
    meta_b <- meta[boot_ids, , drop = FALSE]

    rownames(expr_b) <- paste0("boot", b, "_", seq_len(nrow(expr_b)))
    rownames(meta_b) <- rownames(expr_b)

    boot_stats <- tryCatch(
      fit_fast(expr_b, meta_b, covariates),
      error = function(e) {
        msg("WARNING bootstrap", b, "failed:", conditionMessage(e))
        return(NULL)
      }
    )

    if (is.null(boot_stats)) next

    boot_modules <- derive_modules(boot_stats, thresholds)
    cmp <- compare_modules(boot_modules, frozen_modules, full_stats, boot_stats, iter = b)
    all_records[[length(all_records) + 1]] <- cmp

    for (key in names(frozen_modules)) {
      selected <- intersect(frozen_modules[[key]], boot_modules[[key]])
      if (length(selected) > 0) {
        freq_counter[[key]][selected] <- freq_counter[[key]][selected] + 1L
      }
    }
  }

  stability <- do.call(rbind, all_records)

  stability_out <- file.path(outdir, "01_step12A_bootstrap_module_stability_long.tsv")
  safe_write_tsv(stability, stability_out)

  summary <- aggregate(
    cbind(jaccard, frozen_recall, bootstrap_precision, sign_concordance_overlap, full_vs_boot_signed_score_spearman) ~ module + top_n,
    data = stability,
    FUN = function(x) {
      c(
        median = median(x, na.rm = TRUE),
        mean = mean(x, na.rm = TRUE),
        q05 = quantile(x, 0.05, na.rm = TRUE),
        q95 = quantile(x, 0.95, na.rm = TRUE)
      )
    }
  )

  flat <- data.frame(module = summary$module, top_n = summary$top_n)

  for (col in setdiff(colnames(summary), c("module", "top_n"))) {
    mat <- summary[[col]]
    if (is.matrix(mat)) {
      for (nm in colnames(mat)) {
        flat[[paste0(col, "_", nm)]] <- mat[, nm]
      }
    }
  }

  summary_out <- file.path(outdir, "02_step12A_bootstrap_module_stability_summary.tsv")
  safe_write_tsv(flat, summary_out)

  freq_list <- list()
  ii <- 1

  for (key in names(freq_counter)) {
    program <- sub("\\|top.*$", "", key)
    top_n <- as.integer(sub("^.*\\|top", "", key))
    cnt <- freq_counter[[key]]

    freq_list[[ii]] <- data.frame(
      module = program,
      top_n = top_n,
      gene = names(cnt),
      selected_n = as.integer(cnt),
      selected_fraction = as.integer(cnt) / n_boot,
      stringsAsFactors = FALSE
    )
    ii <- ii + 1
  }

  freq_df <- do.call(rbind, freq_list)
  freq_out <- file.path(outdir, "03_step12A_frozen_gene_selection_frequency.tsv.gz")
  safe_write_tsv(freq_df, freq_out)

  decision <- flat
  decision$decision <- "REVIEW"
  decision$criterion <- "top200/top500: median frozen recall >=0.30 and median sign concordance >=0.90; top50/top100 exploratory"

  decision$decision[
    decision$top_n %in% c(200, 500) &
      decision$frozen_recall_median >= 0.30 &
      decision$sign_concordance_overlap_median >= 0.90
  ] <- "PASS"

  decision$decision[decision$top_n %in% c(50, 100)] <- "INFO"

  decision_out <- file.path(outdir, "04_step12A_decision_table.tsv")
  safe_write_tsv(decision, decision_out)

  diag <- data.frame(
    gandal_rdata = gandal_rdata,
    dx_col = dx_col,
    subject_col = subject_col,
    n_subjects = nrow(meta),
    n_ASD = sum(meta$dx_bootstrap_tmp == "ASD"),
    n_Control = sum(meta$dx_bootstrap_tmp == "Control"),
    n_genes = ncol(expr),
    covariates = paste(covariates, collapse = ","),
    n_boot = n_boot,
    seed = seed,
    stringsAsFactors = FALSE
  )

  diag_out <- file.path(outdir, "00_step12A_input_diagnostics.tsv")
  safe_write_tsv(diag, diag_out)

  md_out <- file.path(outdir, "05_step12A_overall_summary.md")
  con <- file(md_out, open = "wt")

  writeLines("# Step12A discovery native module bootstrap stability", con)
  writeLines("", con)

  writeLines("## Purpose", con)
  writeLines("", con)
  writeLines("This analysis evaluates whether the native NeuroTRACE modules are stable under subject-level discovery-cohort bootstrap resampling and whether topN module thresholds introduce excessive researcher degrees of freedom.", con)
  writeLines("", con)

  writeLines("## Input diagnostics", con)
  writeLines("", con)
  writeLines(sprintf("- Gandal RData: `%s`", gandal_rdata), con)
  writeLines(sprintf("- Diagnosis column: `%s`", dx_col), con)
  writeLines(sprintf("- Subject column: `%s`", subject_col), con)
  writeLines(sprintf("- ASD subjects: %s", sum(meta$dx_bootstrap_tmp == "ASD")), con)
  writeLines(sprintf("- Control subjects: %s", sum(meta$dx_bootstrap_tmp == "Control")), con)
  writeLines(sprintf("- Genes after filtering: %s", ncol(expr)), con)
  writeLines(sprintf("- Covariates: %s", ifelse(length(covariates) == 0, "none", paste(covariates, collapse = ", "))), con)
  writeLines(sprintf("- Bootstrap replicates: %s", n_boot), con)
  writeLines("", con)

  writeLines("## Full model versus frozen module check", con)
  writeLines("", con)
  writeLines(knitr::kable(full_compare, format = "markdown", digits = 4), con)
  writeLines("", con)

  writeLines("## Bootstrap stability summary", con)
  writeLines("", con)
  writeLines(knitr::kable(flat, format = "markdown", digits = 4), con)
  writeLines("", con)

  writeLines("## Decision table", con)
  writeLines("", con)
  writeLines(knitr::kable(decision, format = "markdown", digits = 4), con)
  writeLines("", con)

  writeLines("## Interpretation guide", con)
  writeLines("", con)
  writeLines("- PASS for top200/top500 means the module shows acceptable recovery of frozen genes under subject-level bootstrap resampling and preserves disease-effect direction among overlapping genes.", con)
  writeLines("- INFO for top50/top100 means these are exploratory low-overlap thresholds and should not carry the main manuscript conclusions.", con)
  writeLines("- If top200/top500 are REVIEW, emphasize module-score validation, matched-random transport, graph-prioritized convergence, and external validation rather than exact single-gene membership stability.", con)
  close(con)

  msg("Finished Step12A v4")
  for (p in c(diag_out, full_stats_out, full_compare_out, stability_out, summary_out, freq_out, decision_out, md_out)) {
    msg("Wrote:", p)
  }
}

main()
