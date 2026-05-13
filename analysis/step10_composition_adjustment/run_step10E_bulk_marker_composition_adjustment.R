#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
PROJECT <- file.path(BASE, "neurotrace_algorithm_project")
STEP <- file.path(PROJECT, "step10_composition_adjustment")
OUT <- file.path(STEP, "results")
FIG <- file.path(STEP, "figures")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

MODULE_FILE <- file.path(BASE, "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")
CHOSEN_EXTERNAL <- file.path(BASE, "neurotrace_algorithm_project/step06_cross_cohort_validation/results/32_step06B_chosen_external_inputs.tsv")
GANDAL_RDATA <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"

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
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("", "NA", "NAN", "NULL", "---")] <- NA_character_
  y
}

standardize_dx <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("control|ctl|normal|unaffected", y)] <- "Control"
  out[grepl("^asd$|autism|autistic|idiopathic asd|idiopathic_autism|case", y)] <- "ASD"
  out
}

extract_ensembl <- function(x) {
  y <- as.character(x)
  hit <- regmatches(y, regexpr("ENSG[0-9]+", y))
  out <- rep(NA_character_, length(y))
  ok <- nzchar(hit)
  out[ok] <- hit[ok]
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

map_ensembl_to_symbol_orgdb <- function(ens_ids) {
  ok <- requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!ok) {
    stop("AnnotationDbi/org.Hs.eg.db unavailable, cannot map Gandal Ensembl IDs.")
  }

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
  setnames(mp, c("ENSEMBL", "SYMBOL"), c("ensembl_id_clean", "gene_symbol"), skip_absent = TRUE)
  mp <- mp[!is.na(ensembl_id_clean) & !is.na(gene_symbol)]
  mp[, gene_symbol := clean_gene(gene_symbol)]
  mp <- mp[!is.na(gene_symbol)]
  mp <- unique(mp, by = c("ensembl_id_clean", "gene_symbol"))
  setorder(mp, ensembl_id_clean, gene_symbol)
  mp <- mp[!duplicated(ensembl_id_clean)]
  mp
}

collapse_expression_by_symbol <- function(mat, gene_symbols) {
  gene_symbols <- clean_gene(gene_symbols)
  keep <- !is.na(gene_symbols) & gene_symbols != ""
  mat <- mat[keep, , drop = FALSE]
  gene_symbols <- gene_symbols[keep]

  dt <- as.data.table(mat)
  dt[, gene_symbol := gene_symbols]
  sample_cols <- setdiff(names(dt), "gene_symbol")

  for (cc in sample_cols) {
    suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
  }

  dt[, row_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = sample_cols]
  setorder(dt, gene_symbol, -row_var)
  dt <- dt[!duplicated(gene_symbol)]
  rn <- dt$gene_symbol
  out <- as.matrix(dt[, ..sample_cols])
  rownames(out) <- rn
  out
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
  preferred <- c("datExpr", "expr", "expression", "dat_expression", "rsem_gene_counts")
  nms <- ls(env)
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
  preferred <- c("datMeta", "datMeta_model", "metadata", "meta", "pheno", "sample_metadata")
  nms <- ls(env)

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

  candidates <- list()
  scores <- c()
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
  if (length(candidates) == 0) stop("No metadata-like object found.")
  nm <- names(which.max(scores))
  list(name = nm, object = candidates[[nm]])
}

read_symbol_expression <- function(path) {
  dt <- fread(path)
  gene_col <- names(dt)[1]
  setnames(dt, gene_col, "feature_id")
  dt[, gene_symbol := clean_gene(feature_id)]
  dt <- dt[!is.na(gene_symbol)]
  dt <- dt[!duplicated(gene_symbol)]
  sample_cols <- setdiff(names(dt), c("feature_id", "gene_symbol"))
  for (cc in sample_cols) {
    suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
  }
  mat <- as.matrix(dt[, ..sample_cols])
  rownames(mat) <- dt$gene_symbol
  mat
}

read_external_meta <- function(path) {
  md <- fread(path)
  sample_col <- detect_col(md, c("sample_id", "local_sample_id", "geo_accession", "sample", "id"))
  dx_col <- detect_col(md, c("diagnosis", "disease_status", "disease status", "disease status:ch1", "diagnosis:ch1", "dx", "group"), regex = "diagnosis|disease|dx|group")
  if (is.na(sample_col)) stop("Cannot detect sample_id column in ", path)
  if (is.na(dx_col)) stop("Cannot detect diagnosis column in ", path)

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

prepare_gandal <- function() {
  env <- new.env(parent = emptyenv())
  load(GANDAL_RDATA, envir = env)

  ep <- pick_expression_object(env)
  mp <- pick_metadata_object(env, colnames(as_numeric_matrix(ep$object)))

  mat0 <- as_numeric_matrix(ep$object)
  meta <- as.data.table(mp$object)

  if (is.null(rownames(mat0))) rownames(mat0) <- paste0("Gene_", seq_len(nrow(mat0)))
  if (is.null(colnames(mat0))) colnames(mat0) <- paste0("Sample_", seq_len(ncol(mat0)))

  ens <- extract_ensembl(rownames(mat0))
  if (mean(!is.na(ens)) > 0.5) {
    ens_map <- map_ensembl_to_symbol_orgdb(ens)
    gdt <- data.table(original = rownames(mat0), ensembl_id_clean = ens)
    gdt <- merge(gdt, ens_map, by = "ensembl_id_clean", all.x = TRUE)
    gene_symbols <- gdt$gene_symbol[match(rownames(mat0), gdt$original)]
  } else {
    gene_symbols <- rownames(mat0)
  }

  mat <- collapse_expression_by_symbol(mat0, gene_symbols)

  sample_col <- detect_col(meta, c("sample_id", "SampleID", "sample", "Sample", "id", "ID"))
  if (is.na(sample_col)) {
    meta[, sample_id_detected := rownames(mp$object)]
    sample_col <- "sample_id_detected"
  }
  subject_col <- detect_col(meta, c("subject", "Subject", "subject_id", "SubjectID", "donor", "donor_id"), regex = "subject|donor")
  dx_col <- detect_col(meta, c("Diagnosis", "diagnosis", "Dx", "dx", "DxReg", "condition", "group"), regex = "diagnosis|dx|condition|group")
  if (is.na(subject_col)) {
    meta[, subject_detected := get(sample_col)]
    subject_col <- "subject_detected"
  }
  if (is.na(dx_col)) stop("Gandal diagnosis column not detected.")

  setnames(meta, sample_col, "sample_id", skip_absent = TRUE)
  setnames(meta, subject_col, "subject", skip_absent = TRUE)
  setnames(meta, dx_col, "diagnosis_raw", skip_absent = TRUE)

  meta[, sample_id := as.character(sample_id)]
  meta[, subject := as.character(subject)]
  meta[, dx := standardize_dx(diagnosis_raw)]

  common <- intersect(colnames(mat), meta$sample_id)
  mat <- mat[, common, drop = FALSE]
  meta <- meta[match(common, sample_id)]
  meta <- meta[dx %in% c("ASD", "Control")]
  mat <- mat[, meta$sample_id, drop = FALSE]

  # subject-level aggregation
  subject_ids <- unique(meta$subject)
  subject_mat <- sapply(subject_ids, function(sid) {
    idx <- which(meta$subject == sid)
    if (length(idx) == 1) mat[, idx] else rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  })
  subject_mat <- as.matrix(subject_mat)
  colnames(subject_mat) <- subject_ids
  rownames(subject_mat) <- rownames(mat)

  subject_meta <- meta[, .(
    n_samples = .N,
    dx = names(sort(table(dx), decreasing = TRUE))[1],
    age_cov = suppressWarnings(mean(as.numeric(if ("Age" %in% names(.SD)) Age else if ("age" %in% names(.SD)) age else NA), na.rm = TRUE)),
    rin_cov = suppressWarnings(mean(as.numeric(if ("RIN" %in% names(.SD)) RIN else if ("rin" %in% names(.SD)) rin else NA), na.rm = TRUE)),
    pmi_cov = suppressWarnings(mean(as.numeric(if ("PMI" %in% names(.SD)) PMI else if ("pmi" %in% names(.SD)) pmi else NA), na.rm = TRUE)),
    sex_cov = if ("Sex" %in% names(.SD)) names(sort(table(as.character(Sex)), decreasing = TRUE))[1] else if ("sex" %in% names(.SD)) names(sort(table(as.character(sex)), decreasing = TRUE))[1] else NA_character_
  ), by = subject]

  subject_meta <- subject_meta[match(colnames(subject_mat), subject)]
  subject_meta[, sample_id := subject]

  list(cohort = "Gandal2022", mat = subject_mat, meta = subject_meta, unit_col = "subject")
}

prepare_external <- function(cohort, expr_path, meta_path) {
  mat <- read_symbol_expression(expr_path)
  md <- read_external_meta(meta_path)

  common <- intersect(colnames(mat), md$sample_id)
  mat <- mat[, common, drop = FALSE]
  md <- md[match(common, sample_id)]
  md <- md[dx %in% c("ASD", "Control")]
  mat <- mat[, md$sample_id, drop = FALSE]

  list(cohort = cohort, mat = mat, meta = md, unit_col = "sample_id")
}

row_zscore <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
}

marker_sets <- list(
  EXN = c("SLC17A7", "SLC17A6", "CAMK2A", "CAMK2B", "SATB2", "TBR1", "GRIN2A", "GRIN2B", "VGLUT1", "RORB", "SLC30A3", "NEUROD6"),
  INN = c("GAD1", "GAD2", "SLC6A1", "SST", "PVALB", "VIP", "RELN", "CALB1", "CALB2", "LHX6", "DLX1", "DLX2"),
  AST = c("AQP4", "GFAP", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1", "SOX9", "CLU", "Aldoc", "SPARCL1", "S100B"),
  ODC = c("MBP", "MOG", "PLP1", "MOBP", "MAG", "CNP", "CLDN11", "ERMN", "MAL", "OPALIN"),
  OPC = c("PDGFRA", "CSPG4", "VCAN", "OLIG1", "OLIG2", "SOX10", "BCAN", "NKX2-2", "GPR17"),
  MG  = c("CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "C1QA", "C1QB", "C1QC", "TYROBP", "ITGAM", "TREM2", "SPI1"),
  END = c("CLDN5", "PECAM1", "VWF", "FLT1", "KDR", "ENG", "ESAM", "RAMP2", "ABCB1", "SLC2A1")
)

score_module <- function(mat_z, genes, weights) {
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

score_marker <- function(mat_z, genes) {
  genes <- clean_gene(genes)
  common <- intersect(rownames(mat_z), genes)
  if (length(common) < 2) {
    return(list(score = rep(NA_real_, ncol(mat_z)), n_common = length(common), common = common))
  }
  score <- colMeans(mat_z[common, , drop = FALSE], na.rm = TRUE)
  names(score) <- colnames(mat_z)
  list(score = score, n_common = length(common), common = common)
}

make_model_matrix <- function(md, extra_vars = character()) {
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

  for (ev in extra_vars) {
    if (ev %in% names(md)) {
      vals <- suppressWarnings(as.numeric(md[[ev]]))
      if (sum(is.finite(vals)) >= 5 && sd(vals, na.rm = TRUE) > 0) {
        med <- median(vals, na.rm = TRUE)
        vals[!is.finite(vals)] <- med
        md[, (ev) := as.numeric(scale(vals))]
        covars <- c(covars, ev)
      }
    }
  }

  rhs <- c("dx_bin", covars)
  form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  X <- model.matrix(form, data = md)

  qrX <- qr(X)
  keep <- sort(qrX$pivot[seq_len(qrX$rank)])
  X <- X[, keep, drop = FALSE]
  if (!"dx_bin" %in% colnames(X)) stop("dx_bin missing in model matrix")

  list(md = md, X = X, covars = covars)
}

fit_score_model <- function(y, md, extra_vars = character()) {
  model <- make_model_matrix(md, extra_vars = extra_vars)
  md2 <- model$md
  X <- model$X

  y <- as.numeric(y[md2$sample_id])
  ok <- is.finite(y) & complete.cases(X)
  y <- y[ok]
  X2 <- X[ok, , drop = FALSE]
  md2 <- md2[ok]

  if (length(unique(md2$dx)) < 2) {
    return(data.table(beta = NA_real_, se = NA_real_, t = NA_real_, p_value = NA_real_, n = length(y), n_ASD = sum(md2$dx == "ASD"), n_Control = sum(md2$dx == "Control"), covariates = paste(model$covars, collapse = ";")))
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

  data.table(
    beta = as.numeric(beta),
    se = as.numeric(se),
    t = as.numeric(tval),
    p_value = as.numeric(pval),
    n = length(y),
    n_ASD = sum(md2$dx == "ASD"),
    n_Control = sum(md2$dx == "Control"),
    covariates = paste(model$covars, collapse = ";")
  )
}

message("[", timestamp(), "] Step10E bulk marker-based composition adjustment started")

if (!file.exists(MODULE_FILE)) stop("Missing module file: ", MODULE_FILE)
if (!file.exists(CHOSEN_EXTERNAL)) stop("Missing chosen external input file: ", CHOSEN_EXTERNAL)

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

chosen <- fread(CHOSEN_EXTERNAL)

cohorts <- list()
cohorts[["Gandal2022"]] <- prepare_gandal()

for (i in seq_len(nrow(chosen))) {
  cohorts[[chosen$cohort[i]]] <- prepare_external(
    cohort = chosen$cohort[i],
    expr_path = chosen$expression_path[i],
    meta_path = chosen$metadata_path[i]
  )
}

all_scores <- list()
module_match_rows <- list()
marker_match_rows <- list()
marker_assoc_rows <- list()
ntm_model_rows <- list()
corr_rows <- list()
pca_rows <- list()
pca_loading_rows <- list()

for (coh in names(cohorts)) {
  message("[", timestamp(), "] Processing cohort: ", coh)

  obj <- cohorts[[coh]]
  mat <- obj$mat
  md <- copy(obj$meta)
  md[, sample_id := as.character(sample_id)]

  if (quantile(mat, 0.99, na.rm = TRUE) > 50) {
    mat <- log2(mat + 1)
  }

  finite_rate <- rowMeans(is.finite(mat))
  mat <- mat[finite_rate >= 0.95, , drop = FALSE]
  mat_z <- row_zscore(mat)

  # score NTM modules
  score_dt <- data.table(cohort = coh, sample_id = colnames(mat_z))
  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      msub <- modules[get("program") == prog & top_n == N]
      setorder(msub, rank_within_program)
      msub <- msub[!duplicated(gene_key)]

      sc <- score_module(mat_z, msub$gene_key, as.numeric(msub$weight))
      score_name <- paste0(prog, "_top", N)
      score_dt[, (score_name) := as.numeric(sc$score)]

      module_match_rows[[length(module_match_rows) + 1]] <- data.table(
        cohort = coh,
        program = prog,
        top_n = N,
        n_module_genes = uniqueN(msub$gene_key),
        n_common_genes = sc$n_common,
        overlap_rate = sc$n_common / uniqueN(msub$gene_key),
        common_genes = paste(sc$common, collapse = ";")
      )
    }
  }

  # score marker proxies
  for (ct in names(marker_sets)) {
    sc <- score_marker(mat_z, marker_sets[[ct]])
    score_dt[, (paste0("marker_", ct)) := as.numeric(sc$score)]

    marker_match_rows[[length(marker_match_rows) + 1]] <- data.table(
      cohort = coh,
      cell_type = ct,
      n_marker_genes = length(unique(clean_gene(marker_sets[[ct]]))),
      n_common_genes = sc$n_common,
      overlap_rate = sc$n_common / length(unique(clean_gene(marker_sets[[ct]]))),
      common_genes = paste(sc$common, collapse = ";")
    )
  }

  # merge metadata
  md_keep <- md[dx %in% c("ASD", "Control")]
  score_dt <- merge(score_dt, md_keep, by = "sample_id", all.x = TRUE)
  score_dt <- score_dt[dx %in% c("ASD", "Control")]

  marker_cols <- paste0("marker_", names(marker_sets))
  marker_cols <- marker_cols[marker_cols %in% names(score_dt)]

  # marker PCA
  marker_mat <- as.matrix(score_dt[, ..marker_cols])
  storage.mode(marker_mat) <- "numeric"
  marker_ok <- complete.cases(marker_mat)
  if (sum(marker_ok) >= 10 && ncol(marker_mat) >= 2) {
    marker_scaled <- scale(marker_mat[marker_ok, , drop = FALSE])
    pc <- prcomp(marker_scaled, center = FALSE, scale. = FALSE)

    pc_scores <- matrix(NA_real_, nrow = nrow(score_dt), ncol = min(3, ncol(pc$x)))
    pc_scores[marker_ok, ] <- pc$x[, seq_len(ncol(pc_scores)), drop = FALSE]
    for (j in seq_len(ncol(pc_scores))) {
      score_dt[, paste0("marker_PC", j) := pc_scores[, j]]
    }

    var_expl <- pc$sdev^2 / sum(pc$sdev^2)
    for (j in seq_len(min(5, length(var_expl)))) {
      pca_rows[[length(pca_rows) + 1]] <- data.table(
        cohort = coh,
        PC = paste0("marker_PC", j),
        variance_explained = var_expl[j],
        cumulative_variance_explained = sum(var_expl[seq_len(j)])
      )
    }

    load <- as.data.table(pc$rotation, keep.rownames = "marker")
    load[, cohort := coh]
    pca_loading_rows[[length(pca_loading_rows) + 1]] <- load
  }

  all_scores[[length(all_scores) + 1]] <- score_dt

  # marker-diagnosis association
  for (mk in marker_cols) {
    fit <- fit_score_model(
      y = setNames(score_dt[[mk]], score_dt$sample_id),
      md = score_dt,
      extra_vars = character()
    )
    marker_assoc_rows[[length(marker_assoc_rows) + 1]] <- cbind(
      data.table(cohort = coh, marker = mk, model_type = "marker_dx_base"),
      fit
    )
  }

  # NTM-marker correlations
  ntm_cols <- grep("^NTM", names(score_dt), value = TRUE)
  for (ntm in ntm_cols) {
    for (mk in marker_cols) {
      ok <- is.finite(score_dt[[ntm]]) & is.finite(score_dt[[mk]])
      if (sum(ok) >= 5) {
        ct <- suppressWarnings(cor.test(score_dt[[ntm]][ok], score_dt[[mk]][ok], method = "spearman"))
        corr_rows[[length(corr_rows) + 1]] <- data.table(
          cohort = coh,
          ntm_score = ntm,
          marker = mk,
          n = sum(ok),
          spearman_rho = as.numeric(ct$estimate),
          p_value = as.numeric(ct$p.value)
        )
      }
    }
  }

  # NTM models before/after composition adjustment
  extra_pc <- intersect(c("marker_PC1", "marker_PC2"), names(score_dt))
  extra_all_markers <- marker_cols

  for (ntm in ntm_cols) {
    y <- setNames(score_dt[[ntm]], score_dt$sample_id)

    models <- list(
      M0_base_covariates = character(),
      M1_marker_PC1_PC2 = extra_pc,
      M2_all_marker_scores = extra_all_markers
    )

    for (mt in names(models)) {
      fit <- fit_score_model(y = y, md = score_dt, extra_vars = models[[mt]])
      ntm_model_rows[[length(ntm_model_rows) + 1]] <- cbind(
        data.table(cohort = coh, ntm_score = ntm, model_type = mt),
        fit
      )
    }
  }
}

scores_all <- rbindlist(all_scores, fill = TRUE)
module_match <- rbindlist(module_match_rows, fill = TRUE)
marker_match <- rbindlist(marker_match_rows, fill = TRUE)
marker_assoc <- rbindlist(marker_assoc_rows, fill = TRUE)
ntm_models <- rbindlist(ntm_model_rows, fill = TRUE)
corr_dt <- rbindlist(corr_rows, fill = TRUE)
pca_dt <- rbindlist(pca_rows, fill = TRUE)
pca_loading <- rbindlist(pca_loading_rows, fill = TRUE)

ntm_models[, fdr := p.adjust(p_value, method = "BH"), by = .(model_type)]
marker_assoc[, fdr := p.adjust(p_value, method = "BH")]
corr_dt[, fdr := p.adjust(p_value, method = "BH")]

# Attenuation relative to base model
base <- ntm_models[model_type == "M0_base_covariates", .(
  cohort,
  ntm_score,
  beta_base = beta,
  se_base = se,
  p_base = p_value,
  fdr_base = fdr
)]

adj <- ntm_models[model_type != "M0_base_covariates"]
atten <- merge(adj, base, by = c("cohort", "ntm_score"), all.x = TRUE)
atten[, beta_adjusted := beta]
atten[, attenuation_fraction := fifelse(
  is.finite(beta_base) & abs(beta_base) > 1e-12,
  1 - beta_adjusted / beta_base,
  NA_real_
)]
atten[, direction_preserved := sign(beta_base) == sign(beta_adjusted)]
atten[, abs_beta_ratio_adjusted_vs_base := abs(beta_adjusted) / abs(beta_base)]

safe_fwrite(scores_all, file.path(OUT, "01_step10E_bulk_NTM_marker_scores.tsv.gz"))
safe_fwrite(module_match, file.path(OUT, "02_step10E_NTM_gene_match.tsv"))
safe_fwrite(marker_match, file.path(OUT, "03_step10E_marker_gene_match.tsv"))
safe_fwrite(marker_assoc, file.path(OUT, "04_step10E_marker_dx_association.tsv"))
safe_fwrite(corr_dt, file.path(OUT, "05_step10E_NTM_marker_spearman_correlations.tsv"))
safe_fwrite(pca_dt, file.path(OUT, "06_step10E_marker_PCA_variance.tsv"))
safe_fwrite(pca_loading, file.path(OUT, "07_step10E_marker_PCA_loadings.tsv"))
safe_fwrite(ntm_models, file.path(OUT, "08_step10E_NTM_models_before_after_composition.tsv"))
safe_fwrite(atten, file.path(OUT, "09_step10E_NTM_composition_attenuation.tsv"))

# summary
summary <- atten[, .(
  n_tests = .N,
  n_direction_preserved = sum(direction_preserved, na.rm = TRUE),
  median_abs_beta_ratio = median(abs_beta_ratio_adjusted_vs_base, na.rm = TRUE),
  median_attenuation_fraction = median(attenuation_fraction, na.rm = TRUE),
  n_adjusted_p_lt_0.05 = sum(p_value < 0.05, na.rm = TRUE),
  min_adjusted_p = min(p_value, na.rm = TRUE)
), by = .(model_type)]
safe_fwrite(summary, file.path(OUT, "10_step10E_composition_adjustment_summary.tsv"))

# plot
try({
  pdf(file.path(FIG, "01_step10E_NTM_beta_before_after_composition.pdf"), width = 11, height = 7)
  par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
  for (mt in unique(atten$model_type)) {
    ss <- atten[model_type == mt]
    plot(ss$beta_base, ss$beta_adjusted,
         pch = 16,
         xlab = "Base beta",
         ylab = "Composition-adjusted beta",
         main = mt)
    abline(0, 1, lty = 2, col = "grey")
    abline(h = 0, lty = 3, col = "grey")
    abline(v = 0, lty = 3, col = "grey")
  }
  dev.off()
}, silent = TRUE)

summary_lines <- c(
  "# NeuroTRACE Step10E bulk marker-based composition adjustment summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Assess whether NTM module associations in adult bulk cortex are explained by broad cell-type composition proxies.",
  "",
  "## Cohorts",
  paste(names(cohorts), collapse = ", "),
  "",
  "## Cell-type marker proxies",
  paste(names(marker_sets), collapse = ", "),
  "",
  "## Model types",
  "- M0_base_covariates: NTM ~ diagnosis + available covariates",
  "- M1_marker_PC1_PC2: NTM ~ diagnosis + covariates + marker_PC1 + marker_PC2",
  "- M2_all_marker_scores: NTM ~ diagnosis + covariates + EXN + INN + AST + ODC + OPC + MG + END marker proxies",
  "",
  "## Adjustment summary",
  paste(capture.output(print(summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "If diagnosis betas remain direction-preserved after marker-PC and all-marker adjustment, NTM associations are not fully explained by broad marker-based composition proxies. If betas attenuate strongly, composition coupling should be acknowledged."
)
writeLines(summary_lines, file.path(OUT, "11_step10E_bulk_marker_composition_adjustment_summary.md"))

message("[", timestamp(), "] Step10E done")
message("[", timestamp(), "] Results: ", OUT)
