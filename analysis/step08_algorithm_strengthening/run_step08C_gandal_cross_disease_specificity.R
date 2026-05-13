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

PATH_FILE <- file.path(BASE, "neurotrace_algorithm_project/step01_resource_manifest/results/09_step01B_preferred_paths.sh")
MODULE_FILE <- file.path(BASE, "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")

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

read_path_file <- function(path_file) {
  lines <- readLines(path_file, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[!grepl("^#", lines) & grepl("=", lines)]
  out <- list()
  for (ln in lines) {
    sp <- strsplit(ln, "=", fixed = TRUE)[[1]]
    key <- sp[1]
    val <- paste(sp[-1], collapse = "=")
    val <- gsub('^"|"$', "", val)
    val <- gsub("^'|'$", "", val)
    out[[key]] <- val
  }
  out
}

clean_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

standardize_dx_multi <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))

  out[grepl("control|ctl|normal|unaffected", y)] <- "Control"
  out[grepl("^asd$|autism|autistic|idiopathic asd|idiopathic_autism", y)] <- "ASD"
  out[grepl("schiz|scz", y)] <- "SCZ"
  out[grepl("bipolar|\\bbd\\b", y)] <- "BD"
  out[grepl("dup15|dup 15|duplication", y)] <- "Dup15q"
  out[grepl("affective", y)] <- "AFF"
  out[grepl("major depression|mdd|depression", y)] <- "MDD"

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

is_numeric_matrix_like <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x))) return(FALSE)
  if (nrow(x) < 100 || ncol(x) < 20) return(FALSE)
  xx <- x[seq_len(min(20, nrow(x))), seq_len(min(20, ncol(x))), drop = FALSE]
  suppressWarnings(num <- as.matrix(apply(xx, 2, as.numeric)))
  mean(is.finite(num)) > 0.80
}

as_numeric_matrix <- function(x) {
  if (is.data.frame(x)) {
    rn <- rownames(x)
    x <- as.data.frame(x, check.names = FALSE)
    mat <- as.matrix(data.frame(lapply(x, function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))
    rownames(mat) <- rn
    colnames(mat) <- colnames(x)
    return(mat)
  }
  storage.mode(x) <- "numeric"
  x
}

pick_expression_object <- function(env) {
  nms <- ls(env)
  preferred <- c("datExpr", "expr", "expression", "dat_expression", "rsem_gene_counts")
  for (nm in preferred) {
    if (nm %in% nms && is_numeric_matrix_like(get(nm, envir = env))) {
      return(list(name = nm, object = get(nm, envir = env)))
    }
  }
  candidates <- list()
  for (nm in nms) {
    obj <- get(nm, envir = env)
    if (is_numeric_matrix_like(obj)) candidates[[nm]] <- obj
  }
  if (length(candidates) == 0) stop("No numeric expression-like object found.")
  sizes <- sapply(candidates, function(x) nrow(x) * ncol(x))
  nm <- names(which.max(sizes))
  list(name = nm, object = candidates[[nm]])
}

pick_metadata_object <- function(env, sample_names) {
  nms <- ls(env)
  preferred <- c("datMeta", "datMeta_model", "metadata", "meta", "pheno", "sample_metadata")

  score_meta <- function(x) {
    if (!is.data.frame(x)) return(-Inf)
    cn <- tolower(colnames(x))
    score <- 0
    score <- score + 10 * any(cn %in% c("sample_id", "sampleid", "sample", "id"))
    score <- score + 10 * any(cn %in% c("subject", "subjectid", "donor", "donor_id"))
    score <- score + 15 * any(cn %in% c("diagnosis", "dx", "dxreg", "condition", "group"))
    score <- score + min(20, sum(rownames(x) %in% sample_names))
    score
  }

  for (nm in preferred) {
    if (nm %in% nms && is.data.frame(get(nm, envir = env))) {
      obj <- get(nm, envir = env)
      if (score_meta(obj) > 10) return(list(name = nm, object = obj))
    }
  }

  scores <- c()
  candidates <- list()
  for (nm in nms) {
    obj <- get(nm, envir = env)
    if (is.data.frame(obj)) {
      sc <- score_meta(obj)
      if (sc > 0) {
        candidates[[nm]] <- obj
        scores[nm] <- sc
      }
    }
  }
  if (length(candidates) == 0) stop("No metadata-like data.frame found.")
  nm <- names(which.max(scores))
  list(name = nm, object = candidates[[nm]])
}

extract_ensembl <- function(x) {
  y <- as.character(x)
  hit <- regmatches(y, regexpr("ENSG[0-9]+", y))
  out <- rep(NA_character_, length(y))
  ok <- nzchar(hit)
  out[ok] <- hit[ok]
  out
}

map_ensembl_to_symbol_orgdb <- function(ens_ids) {
  ok <- requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!ok) return(data.table())

  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })

  ens_ids <- unique(na.omit(ens_ids))
  mp <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ens_ids,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  mp <- as.data.table(mp)
  setnames(mp, c("ENSEMBL", "SYMBOL"), c("ensembl_id_clean", "gene_symbol_fixed"), skip_absent = TRUE)
  mp <- mp[!is.na(ensembl_id_clean) & !is.na(gene_symbol_fixed)]
  mp <- unique(mp, by = c("ensembl_id_clean", "gene_symbol_fixed"))
  mp[, source_rank := 1L]
  setorder(mp, ensembl_id_clean, source_rank)
  mp <- mp[!duplicated(ensembl_id_clean)]
  mp
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
  if (length(common) < 5) return(list(score = rep(NA_real_, ncol(mat_z)), n_common = length(common), common = common))
  w <- weights[match(common, genes)]
  score <- as.numeric(crossprod(w, mat_z[common, , drop = FALSE]) / (sum(abs(w)) + 1e-12))
  names(score) <- colnames(mat_z)
  list(score = score, n_common = length(common), common = common)
}

make_model_matrix <- function(md, target_dx) {
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

message("[", timestamp(), "] Step08C Gandal cross-disease specificity started")

paths <- read_path_file(PATH_FILE)
rdata <- paths[["NEUROTRACE_GANDAL2022_NORMALIZED_RDATA"]]
if (is.null(rdata) || !file.exists(rdata)) stop("Missing Gandal RData path.")

if (!file.exists(MODULE_FILE)) stop("Missing NTM module file: ", MODULE_FILE)

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

env <- new.env(parent = emptyenv())
load(rdata, envir = env)

expr_pick <- pick_expression_object(env)
expr <- as_numeric_matrix(expr_pick$object)

meta_pick <- pick_metadata_object(env, colnames(expr))
meta <- as.data.table(meta_pick$object)

if (is.null(rownames(expr))) rownames(expr) <- paste0("Gene_", seq_len(nrow(expr)))
if (is.null(colnames(expr))) colnames(expr) <- paste0("Sample_", seq_len(ncol(expr)))

sample_col <- detect_col(meta, c("sample_id", "SampleID", "sample", "Sample", "id", "ID"))
if (is.na(sample_col)) {
  meta[, sample_id_detected := rownames(meta_pick$object)]
  sample_col <- "sample_id_detected"
}

subject_col <- detect_col(meta, c("subject", "Subject", "subject_id", "SubjectID", "donor", "donor_id"), regex = "subject|donor")
dx_candidates <- c("Diagnosis", "diagnosis", "Dx", "dx", "DxReg", "condition", "group", "Details_on_Diagnosis")
dx_col <- detect_col(meta, dx_candidates, regex = "diagnosis|dx|condition|group")

if (is.na(subject_col)) {
  meta[, subject_detected := get(sample_col)]
  subject_col <- "subject_detected"
}
if (is.na(dx_col)) stop("Could not detect diagnosis column.")

setnames(meta, sample_col, "sample_id", skip_absent = TRUE)
setnames(meta, subject_col, "subject", skip_absent = TRUE)
setnames(meta, dx_col, "diagnosis_raw", skip_absent = TRUE)

meta[, sample_id := as.character(sample_id)]
meta[, subject := as.character(subject)]
meta[, dx_multi := standardize_dx_multi(diagnosis_raw)]

# try to keep multi-diagnosis information if Diagnosis is too coarse
diagnosis_cols <- intersect(c("Diagnosis", "Dx", "DxReg", "Details_on_Diagnosis"), names(meta))
dx_raw_counts <- rbindlist(lapply(diagnosis_cols, function(cc) {
  data.table(column = cc, value = as.character(meta[[cc]]))[, .N, by = .(column, value)]
}), fill = TRUE)
safe_fwrite(dx_raw_counts, file.path(OUT, "55_step08C_gandal_raw_diagnosis_value_counts.tsv"))

common_samples <- intersect(colnames(expr), meta$sample_id)
expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, sample_id)]
stopifnot(all(meta$sample_id == colnames(expr)))

# map expression rownames to symbols if needed
gene_raw <- rownames(expr)
ens <- extract_ensembl(gene_raw)
if (mean(!is.na(ens)) > 0.5) {
  mp <- map_ensembl_to_symbol_orgdb(ens)
  gene_dt <- data.table(original_gene_id = gene_raw, ensembl_id_clean = ens)
  gene_dt <- merge(gene_dt, mp[, .(ensembl_id_clean, gene_symbol_fixed)], by = "ensembl_id_clean", all.x = TRUE)
  gene_dt[, gene_final := fifelse(!is.na(gene_symbol_fixed), gene_symbol_fixed, original_gene_id)]
} else {
  gene_dt <- data.table(original_gene_id = gene_raw, gene_final = gene_raw)
}
gene_dt[, gene_final := clean_gene(gene_final)]
keep <- gene_dt$gene_final != "" & !duplicated(gene_dt$gene_final)
expr <- expr[keep, , drop = FALSE]
rownames(expr) <- gene_dt$gene_final[keep]

# subject-level aggregation
meta <- meta[!is.na(dx_multi)]
expr <- expr[, meta$sample_id, drop = FALSE]

subject_ids <- unique(meta$subject)
subject_expr <- sapply(subject_ids, function(sid) {
  idx <- which(meta$subject == sid)
  if (length(idx) == 1) expr[, idx] else rowMeans(expr[, idx, drop = FALSE], na.rm = TRUE)
})
subject_expr <- as.matrix(subject_expr)
colnames(subject_expr) <- subject_ids
rownames(subject_expr) <- rownames(expr)

subject_meta <- meta[, .(
  n_samples = .N,
  dx_multi = names(sort(table(dx_multi), decreasing = TRUE))[1],
  diagnosis_raw_values = paste(unique(as.character(diagnosis_raw)), collapse = "|"),
  age_cov = suppressWarnings(mean(as.numeric(if ("Age" %in% names(.SD)) Age else if ("age" %in% names(.SD)) age else NA), na.rm = TRUE)),
  rin_cov = suppressWarnings(mean(as.numeric(if ("RIN" %in% names(.SD)) RIN else if ("rin" %in% names(.SD)) rin else NA), na.rm = TRUE)),
  pmi_cov = suppressWarnings(mean(as.numeric(if ("PMI" %in% names(.SD)) PMI else if ("pmi" %in% names(.SD)) pmi else NA), na.rm = TRUE)),
  sex_cov = if ("Sex" %in% names(.SD)) names(sort(table(as.character(Sex)), decreasing = TRUE))[1] else if ("sex" %in% names(.SD)) names(sort(table(as.character(sex)), decreasing = TRUE))[1] else NA_character_
), by = subject]

subject_meta <- subject_meta[match(colnames(subject_expr), subject)]
stopifnot(all(subject_meta$subject == colnames(subject_expr)))

safe_fwrite(subject_meta, file.path(OUT, "56_step08C_gandal_subject_multi_dx_metadata.tsv"))

dx_summary <- subject_meta[, .(
  n_subjects = .N,
  n_samples = sum(n_samples)
), by = dx_multi][order(-n_subjects)]
safe_fwrite(dx_summary, file.path(OUT, "57_step08C_gandal_multi_dx_subject_summary.tsv"))

# Score NTM modules in all subjects
if (quantile(subject_expr, 0.99, na.rm = TRUE) > 50) {
  subject_expr <- log2(subject_expr + 1)
}
finite_rate <- rowMeans(is.finite(subject_expr))
subject_expr <- subject_expr[finite_rate >= 0.95, , drop = FALSE]
mat_z <- row_zscore(subject_expr)

score_rows <- list()
match_rows <- list()

for (prog in unique(modules$program)) {
  for (N in MAIN_TOP_N) {
    msub <- modules[get("program") == prog & top_n == N]
    setorder(msub, rank_within_program)
    msub <- msub[!duplicated(gene_key)]
    msub <- msub[gene_key %in% rownames(mat_z)]

    sc <- score_weighted(mat_z, msub$gene_key, as.numeric(msub$weight))

    match_rows[[length(match_rows) + 1]] <- data.table(
      program = prog,
      top_n = N,
      n_module_genes = uniqueN(modules[get("program") == prog & top_n == N, gene_key]),
      n_common_genes = sc$n_common,
      overlap_rate = sc$n_common / uniqueN(modules[get("program") == prog & top_n == N, gene_key])
    )

    score_rows[[length(score_rows) + 1]] <- data.table(
      subject = colnames(mat_z),
      program = prog,
      top_n = N,
      ntm_score = as.numeric(sc$score)
    )
  }
}

scores <- rbindlist(score_rows, fill = TRUE)
matches <- rbindlist(match_rows, fill = TRUE)
scores <- merge(scores, subject_meta[, .(subject, dx_multi)], by = "subject", all.x = TRUE)

safe_fwrite(matches, file.path(OUT, "58_step08C_gandal_NTM_gene_match.tsv"))
safe_fwrite(scores, file.path(OUT, "59_step08C_gandal_multi_dx_NTM_scores.tsv.gz"))

# Fit disease-vs-control models
available_disease <- setdiff(unique(subject_meta$dx_multi), "Control")
available_disease <- available_disease[!is.na(available_disease)]

model_rows <- list()

for (target_dx in available_disease) {
  md_sub <- subject_meta[dx_multi %in% c("Control", target_dx)]
  if (md_sub[dx_multi == target_dx, .N] < 5 || md_sub[dx_multi == "Control", .N] < 5) next

  model <- make_model_matrix(md_sub, target_dx = target_dx)
  md2 <- model$md
  X <- model$X

  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      sc <- scores[program == prog & top_n == N]
      sc <- sc[match(md2$subject, subject)]
      ft <- fit_beta(sc$ntm_score, md2, X)

      model_rows[[length(model_rows) + 1]] <- data.table(
        target_dx = target_dx,
        program = prog,
        top_n = N,
        beta_target_vs_Control = unname(ft["beta"]),
        se = unname(ft["se"]),
        t = unname(ft["t"]),
        p_value = unname(ft["p"]),
        n = unname(ft["n"]),
        n_target = unname(ft["n_case"]),
        n_Control = unname(ft["n_control"]),
        covariates = paste(model$covars, collapse = ";")
      )
    }
  }
}

models <- rbindlist(model_rows, fill = TRUE)
if (nrow(models) > 0) {
  models[, fdr := p.adjust(p_value, method = "BH")]
  models[, direction_positive := beta_target_vs_Control > 0]
}
safe_fwrite(models, file.path(OUT, "60_step08C_gandal_cross_disease_NTM_models.tsv"))

# Specificity: compare ASD effect to other disease effects
if (nrow(models) > 0 && "ASD" %in% models$target_dx) {
  asd <- models[target_dx == "ASD", .(
    program, top_n,
    ASD_beta = beta_target_vs_Control,
    ASD_p = p_value
  )]
  other <- models[target_dx != "ASD"]
  cmp <- merge(other, asd, by = c("program", "top_n"), all.x = TRUE)
  cmp[, delta_ASD_minus_other_beta := ASD_beta - beta_target_vs_Control]
  cmp[, ASD_beta_greater_than_other := ASD_beta > beta_target_vs_Control]
  safe_fwrite(cmp, file.path(OUT, "61_step08C_ASD_vs_other_disease_effect_comparison.tsv"))

  cmp_summary <- cmp[, .(
    n_other_diseases = .N,
    n_ASD_beta_greater = sum(ASD_beta_greater_than_other, na.rm = TRUE),
    mean_delta_ASD_minus_other = mean(delta_ASD_minus_other_beta, na.rm = TRUE)
  ), by = .(program, top_n)]
  safe_fwrite(cmp_summary, file.path(OUT, "62_step08C_ASD_specificity_summary.tsv"))
} else {
  fwrite(data.table(note = "ASD or other disease comparison unavailable"), file.path(OUT, "61_step08C_ASD_vs_other_disease_effect_comparison.tsv"), sep = "\t")
  fwrite(data.table(note = "ASD or other disease comparison unavailable"), file.path(OUT, "62_step08C_ASD_specificity_summary.tsv"), sep = "\t")
}

# Plot
try({
  pdf(file.path(FIG, "17_step08C_cross_disease_NTM_scores.pdf"), width = 11, height = 7)
  par(mfrow = c(3, 2), mar = c(5, 4, 3, 1))
  for (prog in unique(scores$program)) {
    for (N in MAIN_TOP_N) {
      ss <- scores[program == prog & top_n == N]
      boxplot(ntm_score ~ dx_multi, data = ss,
              las = 2, main = paste0(prog, " top", N),
              ylab = "NTM score", xlab = "")
    }
  }
  dev.off()
}, silent = TRUE)

summary_lines <- c(
  "# NeuroTRACE Step08C Gandal cross-disease specificity summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step08C audits whether the Gandal RData contains non-ASD diagnostic groups and tests whether NeuroTRACE-native ASD modules generalize or remain specific across other adult cortical disease groups.",
  "",
  "## Input",
  paste0("- RData: `", rdata, "`"),
  paste0("- Expression object: `", expr_pick$name, "`"),
  paste0("- Metadata object: `", meta_pick$name, "`"),
  "",
  "## Multi-diagnosis subject summary",
  paste(capture.output(print(dx_summary)), collapse = "\n"),
  "",
  "## NTM model results",
  paste(capture.output(print(models)), collapse = "\n"),
  "",
  "## Interpretation",
  "If only ASD and Control are available, Step08C is an audit rather than a cross-disease test. If SCZ/BD/MDD or other groups are available, compare ASD beta with other disease beta to assess disease specificity versus pan-psychiatric generalization."
)
writeLines(summary_lines, file.path(OUT, "63_step08C_gandal_cross_disease_specificity_summary.md"))

message("[", timestamp(), "] Step08C done")
message("[", timestamp(), "] Available dx groups: ", paste(dx_summary$dx_multi, collapse = ", "))
message("[", timestamp(), "] Results: ", OUT)
