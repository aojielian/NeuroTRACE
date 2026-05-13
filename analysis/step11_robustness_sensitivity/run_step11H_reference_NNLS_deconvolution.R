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

clean_gene <- function(x) {
  x <- as.character(x)
  x <- sub("^gene:", "", x)
  x <- sub("^GENE:", "", x)
  x <- trimws(x)
  x
}


harmonize_gene_symbols <- function(x, label = "gene_ids") {
  x0 <- as.character(x)
  x1 <- clean_gene(x0)

  # Critical fix:
  # Gandal rownames can look like ENSG00000000003.14_1.
  # Extract the core Ensembl ID from any string before mapping.
  ens_core <- rep(NA_character_, length(x1))
  has_ens <- grepl("ENSG[0-9]+", x1)
  if (any(has_ens)) {
    m <- regexpr("ENSG[0-9]+", x1[has_ens])
    ens_core[has_ens] <- regmatches(x1[has_ens], m)
  }

  if (any(!is.na(ens_core))) {
    if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)) {

      keytypes_available <- tryCatch(
        AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db),
        error = function(e) character(0)
      )

      if (!"ENSEMBL" %in% keytypes_available) {
        cat(sprintf("[%s] WARNING: ENSEMBL keytype unavailable for %s; keeping original IDs\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), label))
        return(x1)
      }

      valid_keys <- tryCatch(
        AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL"),
        error = function(e) character(0)
      )

      query_keys <- unique(ens_core[!is.na(ens_core)])
      valid_query <- intersect(query_keys, valid_keys)

      if (length(valid_query) == 0) {
        cat(sprintf("[%s] WARNING: %s has %d Ensembl-like IDs but none are valid ENSEMBL keys in current org.Hs.eg.db; keeping original IDs\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    label,
                    length(query_keys)))
        return(x1)
      }

      mapped <- tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          keys = valid_query,
          keytype = "ENSEMBL",
          column = "SYMBOL",
          multiVals = "first"
        ),
        error = function(e) {
          cat(sprintf("[%s] WARNING: mapIds failed for %s: %s\n",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      label,
                      conditionMessage(e)))
          return(NULL)
        }
      )

      if (!is.null(mapped)) {
        pos <- which(!is.na(ens_core) & ens_core %in% names(mapped))
        vals <- as.character(mapped[ens_core[pos]])
        ok <- !is.na(vals) & vals != ""
        if (any(ok)) {
          x1[pos[ok]] <- vals[ok]
        }

        cat(sprintf("[%s] %s: mapped %d/%d Ensembl-like IDs to symbols using extracted ENSG core\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    label,
                    sum(ok),
                    sum(!is.na(ens_core))))
      }
    } else {
      cat(sprintf("[%s] WARNING: %s has Ensembl-like IDs but AnnotationDbi/org.Hs.eg.db is unavailable; keeping original IDs\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  label))
    }
  }

  x1 <- clean_gene(x1)
  x1
}

collapse_duplicate_gene_rows <- function(mat, label = "matrix") {
  genes <- rownames(mat)
  keep <- !is.na(genes) & genes != "" & genes != "NA" & genes != "nan"
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]

  if (!any(duplicated(genes))) return(mat)

  cat(sprintf("[%s] Collapsing duplicated gene symbols in %s: %d duplicated rows\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              label,
              sum(duplicated(genes))))

  idx <- split(seq_len(nrow(mat)), genes)
  collapsed <- do.call(rbind, lapply(idx, function(ii) {
    Matrix::colMeans(mat[ii, , drop = FALSE])
  }))
  collapsed <- Matrix::Matrix(collapsed, sparse = TRUE)
  rownames(collapsed) <- names(idx)
  colnames(collapsed) <- colnames(mat)
  collapsed
}


zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
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

map_broad_celltype <- function(x) {
  y <- as.character(x)
  low <- tolower(y)

  out <- rep("Other", length(y))

  out[grepl("exn|excit|glut|pyram|proj|l2|l3|l4|l5|l6|it|et|ct", low)] <- "EXN"
  out[grepl("inn|inhib|gaba|sst|pvalb|pv|vip|lamp5|sncg", low)] <- "INN"
  out[grepl("astro|ast", low)] <- "AST"
  out[grepl("\\bopc\\b|precursor", low)] <- "OPC"
  out[grepl("oligo|odc|oligodend", low)] <- "ODC"
  out[grepl("micro|\\bmg\\b|myeloid", low)] <- "MG"
  out[grepl("endo|endothelial|vascular|pericyte|vlmc|smooth", low)] <- "END"

  out
}

extract_seurat_matrix <- function(obj) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed in this R environment.")
  }

  assay <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) "RNA")

  mat <- tryCatch({
    Seurat::GetAssayData(obj, assay = assay, slot = "data")
  }, error = function(e1) {
    tryCatch({
      Seurat::GetAssayData(obj, assay = assay, layer = "data")
    }, error = function(e2) {
      tryCatch({
        Seurat::GetAssayData(obj, assay = assay, slot = "counts")
      }, error = function(e3) {
        Seurat::GetAssayData(obj, assay = assay, layer = "counts")
      })
    })
  })

  mat
}

build_reference <- function(ref_rds, outdir, max_cells_per_type = 5000, seed = 20260509) {
  ref_cache <- file.path(outdir, "01_step11H_PsychENCODE_broad_reference_average.tsv.gz")
  marker_cache <- file.path(outdir, "02_step11H_PsychENCODE_reference_marker_genes.tsv")

  if (file.exists(ref_cache) && file.exists(marker_cache)) {
    msg("Reference cache found. Loading:", ref_cache)
    ref_avg <- read.table(gzfile(ref_cache), header = TRUE, sep = "\t", check.names = FALSE)
    rownames(ref_avg) <- ref_avg$gene
    ref_avg$gene <- NULL
    markers <- read.table(marker_cache, header = TRUE, sep = "\t", check.names = FALSE)
    return(list(ref_avg = as.matrix(ref_avg), markers = markers))
  }

  msg("Loading PsychENCODE reference object:", ref_rds)
  obj <- readRDS(ref_rds)

  meta <- obj@meta.data
  mat <- extract_seurat_matrix(obj)

  if (is.null(colnames(mat))) stop("Reference expression matrix has no cell barcodes.")
  if (is.null(rownames(mat))) stop("Reference expression matrix has no gene names.")

  anno_col <- find_col(meta, c(
    "broad_cell_type", "broad_celltype", "cell_type_broad", "annotation",
    "cell_type", "celltype", "class", "subclass", "cluster", "cluster_name"
  ))

  if (is.na(anno_col)) {
    stop("Cannot detect cell-type annotation column in PsychENCODE metadata.")
  }

  msg("Using annotation column:", anno_col)

  meta$cell_barcode_tmp <- rownames(meta)
  meta$broad_cell_type_NNLS <- map_broad_celltype(meta[[anno_col]])

  keep <- meta$broad_cell_type_NNLS %in% c("EXN", "INN", "AST", "ODC", "OPC", "MG", "END")
  meta2 <- meta[keep, , drop = FALSE]
  meta2 <- meta2[meta2$cell_barcode_tmp %in% colnames(mat), , drop = FALSE]

  set.seed(seed)

  sampled_cells <- unlist(lapply(split(meta2$cell_barcode_tmp, meta2$broad_cell_type_NNLS), function(v) {
    if (length(v) > max_cells_per_type) sample(v, max_cells_per_type) else v
  }), use.names = FALSE)

  meta2 <- meta2[match(sampled_cells, meta2$cell_barcode_tmp), , drop = FALSE]
  mat2 <- mat[, sampled_cells, drop = FALSE]

  msg("Reference cells after broad-type filtering/sampling:", ncol(mat2))
  msg("Broad cell-type counts:")
  print(table(meta2$broad_cell_type_NNLS))

  genes <- harmonize_gene_symbols(rownames(mat2), label = "PsychENCODE reference")
  rownames(mat2) <- genes
  mat2 <- collapse_duplicate_gene_rows(mat2, label = "PsychENCODE reference")

  broad_levels <- c("EXN", "INN", "AST", "ODC", "OPC", "MG", "END")
  broad_levels <- broad_levels[broad_levels %in% unique(meta2$broad_cell_type_NNLS)]

  ref_avg <- matrix(NA_real_, nrow = nrow(mat2), ncol = length(broad_levels))
  rownames(ref_avg) <- rownames(mat2)
  colnames(ref_avg) <- broad_levels

  for (ct in broad_levels) {
    cells <- meta2$cell_barcode_tmp[meta2$broad_cell_type_NNLS == ct]
    cells <- intersect(cells, colnames(mat2))
    msg("Averaging reference for", ct, "cells:", length(cells))
    ref_avg[, ct] <- Matrix::rowMeans(mat2[, cells, drop = FALSE])
  }

  # Remove genes with no information.
  keep_gene <- rowSums(is.finite(ref_avg)) == ncol(ref_avg)
  keep_gene <- keep_gene & apply(ref_avg, 1, sd, na.rm = TRUE) > 0
  ref_avg <- ref_avg[keep_gene, , drop = FALSE]

  # Marker selection: top specificity per broad cell type.
  marker_list <- list()
  for (ct in colnames(ref_avg)) {
    others <- setdiff(colnames(ref_avg), ct)
    spec <- ref_avg[, ct] - apply(ref_avg[, others, drop = FALSE], 1, max, na.rm = TRUE)
    ord <- order(spec, decreasing = TRUE)
    top <- head(ord[is.finite(spec[ord]) & spec[ord] > 0], 250)
    marker_list[[ct]] <- data.frame(
      cell_type = ct,
      gene = rownames(ref_avg)[top],
      specificity = spec[top],
      mean_expr = ref_avg[top, ct],
      stringsAsFactors = FALSE
    )
  }

  markers <- do.call(rbind, marker_list)
  markers <- markers[!duplicated(paste(markers$cell_type, markers$gene)), , drop = FALSE]

  ref_out <- data.frame(gene = rownames(ref_avg), ref_avg, check.names = FALSE)
  safe_write_tsv(ref_out, ref_cache)
  safe_write_tsv(markers, marker_cache)

  msg("Wrote reference average:", ref_cache)
  msg("Wrote reference markers:", marker_cache)

  list(ref_avg = ref_avg, markers = markers)
}

load_gandal_bulk <- function(gandal_rdata) {
  msg("Loading Gandal RData:", gandal_rdata)
  env <- new.env()
  load(gandal_rdata, envir = env)

  obj_names <- ls(env)
  msg("Objects in Gandal RData:", paste(obj_names, collapse = ", "))

  # Find expression matrix/data.frame.
  expr_candidates <- list()
  for (nm in obj_names) {
    x <- get(nm, env)
    if ((is.matrix(x) || is.data.frame(x)) && nrow(x) > 1000 && ncol(x) > 20) {
      numeric_frac <- mean(sapply(as.data.frame(x[seq_len(min(20, nrow(x))), seq_len(min(20, ncol(x))), drop = FALSE]), is.numeric))
      expr_candidates[[nm]] <- list(object = x, numeric_frac = numeric_frac, dim = dim(x))
    }
  }

  if (length(expr_candidates) == 0) stop("Cannot detect expression matrix in Gandal RData.")

  # Prefer datExpr if available.
  expr_name <- if ("datExpr" %in% names(expr_candidates)) "datExpr" else names(expr_candidates)[1]
  expr <- as.matrix(expr_candidates[[expr_name]]$object)
  msg("Using expression object:", expr_name, "dim:", paste(dim(expr), collapse = " x "))

  # Find metadata.
  meta_candidates <- list()
  for (nm in obj_names) {
    x <- get(nm, env)
    if (is.data.frame(x) && nrow(x) > 20 && ncol(x) >= 2) {
      meta_candidates[[nm]] <- x
    }
  }

  meta_name <- NA_character_
  for (nm in names(meta_candidates)) {
    x <- meta_candidates[[nm]]
    if (nrow(x) == nrow(expr) || nrow(x) == ncol(expr)) {
      meta_name <- nm
      break
    }
  }
  if (is.na(meta_name) && "datMeta_model" %in% names(meta_candidates)) meta_name <- "datMeta_model"
  if (is.na(meta_name)) stop("Cannot detect metadata object in Gandal RData.")

  meta <- meta_candidates[[meta_name]]
  msg("Using metadata object:", meta_name, "dim:", paste(dim(meta), collapse = " x "))

  # Orient expression to genes x samples.
  # If rownames match metadata rows, expr is samples x genes, so transpose.
  rn <- rownames(expr)
  cn <- colnames(expr)
  meta_ids <- rownames(meta)

  if (!is.null(rn) && !is.null(meta_ids) && length(intersect(rn, meta_ids)) > 10) {
    expr <- t(expr)
  } else if (!is.null(cn) && !is.null(meta_ids) && length(intersect(cn, meta_ids)) > 10) {
    # already genes x samples
    expr <- expr
  } else {
    # Heuristic: more genes than samples -> genes x samples; otherwise transpose.
    if (nrow(expr) < ncol(expr)) expr <- t(expr)
  }

  rownames(expr) <- harmonize_gene_symbols(rownames(expr), label = "Gandal bulk")
  expr <- collapse_duplicate_gene_rows(expr, label = "Gandal bulk")

  # Align metadata.
  common_samples <- intersect(colnames(expr), rownames(meta))
  if (length(common_samples) < 20) {
    # Try sample ID columns.
    sid_col <- find_col(meta, c("sample", "sample_id", "SampleID", "subject", "individual", "id"))
    if (!is.na(sid_col)) {
      rownames(meta) <- as.character(meta[[sid_col]])
      common_samples <- intersect(colnames(expr), rownames(meta))
    }
  }

  if (length(common_samples) < 20) {
    stop("Cannot align Gandal expression samples with metadata.")
  }

  expr <- expr[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]

  list(dataset = "Gandal2022", expr = expr, meta = meta)
}

prepare_dx <- function(meta) {
  dx_col <- find_col(meta, c("Diagnosis", "diagnosis", "dx", "Dx", "group", "Group", "diagnosis_raw", "case_control"))
  if (is.na(dx_col)) stop("Cannot detect diagnosis column.")

  raw <- as.character(meta[[dx_col]])
  low <- tolower(raw)

  dx <- rep(NA_character_, length(raw))
  dx[grepl("asd|autism", low)] <- "ASD"
  dx[grepl("control|ctl|normal|unaffected", low)] <- "Control"

  # Numeric fallback sometimes 2=ASD, 1=Control, but use only if labels not found.
  if (all(is.na(dx))) {
    dx[raw %in% c("2", "case", "Case")] <- "ASD"
    dx[raw %in% c("1", "0", "control", "Control")] <- "Control"
  }

  dx
}

get_covariates <- function(meta) {
  cov_candidates <- list(
    age = c("age", "Age", "age_cov"),
    sex = c("sex", "Sex", "gender", "Gender", "sex_cov"),
    rin = c("RIN", "rin", "RIN_cov", "rin_cov"),
    pmi = c("PMI", "pmi", "PMI_cov", "pmi_cov", "postmortem_interval")
  )

  covs <- c()
  for (nm in names(cov_candidates)) {
    cc <- find_col(meta, cov_candidates[[nm]])
    if (!is.na(cc)) covs <- c(covs, cc)
  }
  unique(covs)
}

softmax_prop <- function(theta) {
  z <- theta - max(theta)
  p <- exp(z)
  p / sum(p)
}

deconv_one_sample <- function(y, A) {
  # A: genes x celltypes; y: genes
  k <- ncol(A)
  obj <- function(theta) {
    p <- softmax_prop(theta)
    pred <- as.vector(A %*% p)
    sum((pred - y)^2, na.rm = TRUE)
  }
  fit <- tryCatch(
    optim(rep(0, k), obj, method = "BFGS", control = list(maxit = 500)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(rep(NA_real_, k))
  }
  softmax_prop(fit$par)
}

recompute_markers_from_refavg <- function(ref_avg, top_n_per_type = 300) {
  ref_avg <- as.matrix(ref_avg)
  out <- list()
  ct_names <- colnames(ref_avg)

  for (ct in ct_names) {
    others <- setdiff(ct_names, ct)
    if (length(others) == 0) next

    spec <- ref_avg[, ct] - apply(ref_avg[, others, drop = FALSE], 1, max, na.rm = TRUE)
    spec[!is.finite(spec)] <- NA_real_

    ord <- order(spec, decreasing = TRUE, na.last = NA)
    ord <- ord[spec[ord] > 0]

    if (length(ord) == 0) next

    use <- head(ord, top_n_per_type)
    out[[ct]] <- data.frame(
      cell_type = ct,
      gene = rownames(ref_avg)[use],
      specificity = spec[use],
      mean_expr = ref_avg[use, ct],
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) {
    return(data.frame(cell_type = character(), gene = character(), specificity = numeric(), mean_expr = numeric()))
  }

  ans <- do.call(rbind, out)
  ans <- ans[!duplicated(paste(ans$cell_type, ans$gene)), , drop = FALSE]
  rownames(ans) <- NULL
  ans
}

run_deconvolution <- function(expr, ref_avg, markers, dataset, outdir) {
  rownames(expr) <- harmonize_gene_symbols(rownames(expr), label = paste0(dataset, " bulk before NNLS"))
  rownames(ref_avg) <- harmonize_gene_symbols(rownames(ref_avg), label = "reference average before NNLS")

  expr <- collapse_duplicate_gene_rows(expr, label = paste0(dataset, " bulk before NNLS"))
  ref_avg <- collapse_duplicate_gene_rows(ref_avg, label = "reference average before NNLS")

  common_ref_bulk <- intersect(rownames(expr), rownames(ref_avg))

  # Recompute markers after harmonizing reference rownames.
  markers_recomputed <- recompute_markers_from_refavg(ref_avg, top_n_per_type = 300)
  markers_recomputed$gene <- clean_gene(markers_recomputed$gene)
  marker_genes_recomputed <- unique(markers_recomputed$gene)

  # Also keep original marker file as diagnostic only.
  markers$gene <- harmonize_gene_symbols(markers$gene, label = "original reference marker genes before NNLS")
  marker_genes_original <- unique(clean_gene(markers$gene))

  common_recomputed_marker <- intersect(common_ref_bulk, marker_genes_recomputed)
  common_original_marker <- intersect(common_ref_bulk, marker_genes_original)

  gene_selection_mode <- "recomputed_reference_markers"
  common <- common_recomputed_marker

  # If recomputed markers are still too sparse but ref-bulk overlap is sufficient,
  # fall back to all shared genes with high reference variance. This is a conservative
  # reference-based deconvolution sensitivity, not a marker-only method.
  if (length(common) < 100 && length(common_ref_bulk) >= 500) {
    gene_selection_mode <- "fallback_all_shared_high_reference_variance"
    ref_shared <- as.matrix(ref_avg[common_ref_bulk, , drop = FALSE])
    rv <- apply(ref_shared, 1, sd, na.rm = TRUE)
    rv[!is.finite(rv)] <- 0
    ord <- order(rv, decreasing = TRUE)
    common <- common_ref_bulk[ord]
    common <- head(common, min(length(common), 3000))
  }

  diag_path <- file.path(outdir, "00_step11H_gene_id_overlap_diagnostics.tsv")
  diag_df <- data.frame(
    dataset = dataset,
    n_bulk_genes = length(unique(rownames(expr))),
    n_reference_genes = length(unique(rownames(ref_avg))),
    n_ref_bulk_overlap = length(common_ref_bulk),
    n_original_marker_genes = length(unique(marker_genes_original)),
    n_recomputed_marker_genes = length(unique(marker_genes_recomputed)),
    n_ref_bulk_original_marker_overlap = length(common_original_marker),
    n_ref_bulk_recomputed_marker_overlap = length(common_recomputed_marker),
    n_final_deconvolution_genes = length(common),
    gene_selection_mode = gene_selection_mode,
    bulk_first20 = paste(head(rownames(expr), 20), collapse = ","),
    reference_first20 = paste(head(rownames(ref_avg), 20), collapse = ","),
    original_marker_first20 = paste(head(marker_genes_original, 20), collapse = ","),
    recomputed_marker_first20 = paste(head(marker_genes_recomputed, 20), collapse = ","),
    final_gene_first20 = paste(head(common, 20), collapse = ","),
    stringsAsFactors = FALSE
  )
  safe_write_tsv(diag_df, diag_path)

  marker_recomputed_out <- file.path(outdir, "02b_step11H_recomputed_reference_marker_genes.tsv")
  safe_write_tsv(markers_recomputed, marker_recomputed_out)

  msg(dataset, "gene-id overlap diagnostics written:", diag_path)
  msg(dataset, "bulk-reference overlap:", length(common_ref_bulk))
  msg(dataset, "original marker overlap:", length(common_original_marker))
  msg(dataset, "recomputed marker overlap:", length(common_recomputed_marker))
  msg(dataset, "final deconvolution genes:", length(common))
  msg(dataset, "gene selection mode:", gene_selection_mode)

  if (length(common_ref_bulk) < 500) {
    stop("Too few genes shared between bulk and reference after gene-ID harmonization: ",
         length(common_ref_bulk), ". See diagnostics: ", diag_path)
  }

  if (length(common) < 100) {
    stop("Too few final genes for reference-based deconvolution: ",
         length(common), ". See diagnostics: ", diag_path)
  }

  ref_use <- ref_avg[common, , drop = FALSE]
  expr_use <- expr[common, , drop = FALSE]

  # Joint gene-wise z-score for compatibility between single-nucleus reference and bulk.
  combined <- cbind(ref_use, expr_use)
  combined_z <- zscore_rows(combined)

  A <- combined_z[, colnames(ref_use), drop = FALSE]
  Y <- combined_z[, colnames(expr_use), drop = FALSE]

  props <- matrix(NA_real_, nrow = ncol(Y), ncol = ncol(A))
  rownames(props) <- colnames(Y)
  colnames(props) <- paste0("NNLS_", colnames(A))

  for (i in seq_len(ncol(Y))) {
    if (i %% 20 == 0) msg(dataset, "deconvolving sample", i, "/", ncol(Y))
    props[i, ] <- deconv_one_sample(Y[, i], A)
  }

  props_df <- data.frame(dataset = dataset, sample_id = rownames(props), props, check.names = FALSE)
  props_df
}

load_module_weights <- function(module_weight_file) {
  mw <- read.table(module_weight_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
  gene_col <- find_col(mw, c("gene_symbol_fixed", "gene_symbol", "gene", "symbol"))
  program_col <- find_col(mw, c("program", "module_family", "module"))
  topn_col <- find_col(mw, c("top_n", "topn"))
  weight_col <- find_col(mw, c("weight", "module_weight", "signed_weight", "score"))

  if (any(is.na(c(gene_col, program_col, topn_col, weight_col)))) {
    stop("Cannot detect required columns in module weight file.")
  }

  out <- data.frame(
    gene = clean_gene(mw[[gene_col]]),
    program = as.character(mw[[program_col]]),
    top_n = as.numeric(mw[[topn_col]]),
    weight = as.numeric(mw[[weight_col]]),
    stringsAsFactors = FALSE
  )

  out <- out[out$program %in% c("NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed") &
               out$top_n %in% c(200, 500) &
               is.finite(out$weight), , drop = FALSE]
  out
}

score_modules <- function(expr, module_weights, dataset) {
  expr_z <- zscore_rows(expr)
  records <- list()

  idx <- 1
  for (key in split(module_weights, paste(module_weights$program, module_weights$top_n, sep = "|"))) {
    genes <- intersect(key$gene, rownames(expr_z))
    if (length(genes) < 20) next

    w <- key$weight[match(genes, key$gene)]
    names(w) <- genes

    mat <- expr_z[genes, , drop = FALSE]
    score <- as.numeric(crossprod(w, mat)) / sum(abs(w))
    records[[idx]] <- data.frame(
      dataset = dataset,
      sample_id = colnames(expr_z),
      program = key$program[1],
      top_n = key$top_n[1],
      n_common_genes = length(genes),
      ntm_score = score,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }

  do.call(rbind, records)
}

fit_one_model <- function(df, outcome, model_type, covariates, prop_cols) {
  dat <- df
  dat$dx_binary <- ifelse(dat$dx == "ASD", 1, ifelse(dat$dx == "Control", 0, NA))
  dat <- dat[is.finite(dat$dx_binary) & is.finite(dat[[outcome]]), , drop = FALSE]

  covs <- covariates[covariates %in% colnames(dat)]
  # Drop covariates with no variation or excessive missingness.
  keep_covs <- c()
  for (cc in covs) {
    v <- dat[[cc]]
    miss <- mean(is.na(v))
    if (miss > 0.3) next
    if (length(unique(v[!is.na(v)])) <= 1) next
    keep_covs <- c(keep_covs, cc)
  }
  covs <- keep_covs

  if (model_type == "base") {
    rhs <- c("dx_binary", covs)
  } else if (model_type == "NNLS_PC12") {
    pc_cols <- c("NNLS_PC1", "NNLS_PC2")
    rhs <- c("dx_binary", pc_cols[pc_cols %in% colnames(dat)], covs)
  } else if (model_type == "NNLS_props") {
    # Drop one proportion to reduce collinearity.
    pp <- prop_cols[prop_cols %in% colnames(dat)]
    if (length(pp) > 1) pp <- pp[-length(pp)]
    rhs <- c("dx_binary", pp, covs)
  } else {
    stop("Unknown model type: ", model_type)
  }

  rhs <- unique(rhs)
  f <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))

  fit <- tryCatch(lm(f, data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(
      model_type = model_type, beta_ASD = NA_real_, se_ASD = NA_real_,
      p_ASD = NA_real_, n = nrow(dat), n_ASD = sum(dat$dx == "ASD"),
      n_Control = sum(dat$dx == "Control"), formula = paste(deparse(f), collapse = " ")
    ))
  }

  sm <- summary(fit)$coefficients
  if (!"dx_binary" %in% rownames(sm)) {
    beta <- se <- p <- NA_real_
  } else {
    beta <- sm["dx_binary", "Estimate"]
    se <- sm["dx_binary", "Std. Error"]
    p <- sm["dx_binary", "Pr(>|t|)"]
  }

  data.frame(
    model_type = model_type,
    beta_ASD = beta,
    se_ASD = se,
    p_ASD = p,
    n = nrow(model.frame(fit)),
    n_ASD = sum(model.frame(fit)$dx_binary == 1),
    n_Control = sum(model.frame(fit)$dx_binary == 0),
    formula = paste(deparse(f), collapse = " "),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  ref_rds <- get_arg("--ref_rds", "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object.rds")
  gandal_rdata <- get_arg("--gandal_rdata", "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData")
  module_weight_file <- get_arg("--module_weights", "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")
  outdir <- get_arg("--outdir", "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step11_robustness_sensitivity/results")
  max_cells_per_type <- as.integer(get_arg("--max_cells_per_type", "5000"))
  seed <- as.integer(get_arg("--seed", "20260509"))

  safe_dir_create(outdir)

  msg("Step11H started")
  msg("Reference RDS:", ref_rds)
  msg("Gandal RData:", gandal_rdata)
  msg("Module weights:", module_weight_file)
  msg("Outdir:", outdir)
  msg("Max cells per type:", max_cells_per_type)

  ref <- build_reference(ref_rds, outdir, max_cells_per_type = max_cells_per_type, seed = seed)
  ref_avg <- ref$ref_avg
  markers <- ref$markers

  bulk <- load_gandal_bulk(gandal_rdata)
  expr <- bulk$expr
  meta <- bulk$meta
  dataset <- bulk$dataset

  dx <- prepare_dx(meta)
  meta$dx_NNLS <- dx
  keep <- meta$dx_NNLS %in% c("ASD", "Control")
  meta <- meta[keep, , drop = FALSE]
  expr <- expr[, rownames(meta), drop = FALSE]

  msg("Gandal ASD/Control samples:", ncol(expr))
  print(table(meta$dx_NNLS))

  props <- run_deconvolution(expr, ref_avg, markers, dataset, outdir)

  # Add deconvolution PCs.
  prop_cols <- grep("^NNLS_", colnames(props), value = TRUE)
  pc <- prcomp(props[, prop_cols, drop = FALSE], center = TRUE, scale. = TRUE)
  props$NNLS_PC1 <- pc$x[, 1]
  props$NNLS_PC2 <- if (ncol(pc$x) >= 2) pc$x[, 2] else 0

  prop_out <- file.path(outdir, "03_step11H_Gandal_NNLS_deconvolution_proportions.tsv")
  safe_write_tsv(props, prop_out)

  # NTM scoring.
  mw <- load_module_weights(module_weight_file)
  scores <- score_modules(expr, mw, dataset)
  score_out <- file.path(outdir, "04_step11H_Gandal_NTM_scores_from_bulk.tsv")
  safe_write_tsv(scores, score_out)

  # Merge modeling data.
  meta$sample_id <- rownames(meta)
  meta$dx <- meta$dx_NNLS

  dat <- merge(scores, meta, by = "sample_id", all.x = TRUE)
  dat <- merge(dat, props, by = c("dataset", "sample_id"), all.x = TRUE)

  model_input_out <- file.path(outdir, "04b_step11H_Gandal_model_input.tsv.gz")
  safe_write_tsv(dat, model_input_out)

  covariates <- get_covariates(meta)
  msg("Detected covariates:", paste(covariates, collapse = ", "))

  model_records <- list()
  idx <- 1
  for (program in unique(dat$program)) {
    for (top_n in sort(unique(dat$top_n))) {
      sub <- dat[dat$program == program & dat$top_n == top_n, , drop = FALSE]
      if (nrow(sub) < 20) next

      for (mt in c("base", "NNLS_PC12", "NNLS_props")) {
        res <- fit_one_model(sub, "ntm_score", mt, covariates, prop_cols)
        res$dataset <- dataset
        res$program <- program
        res$top_n <- top_n
        res$n_common_genes <- unique(sub$n_common_genes)[1]
        model_records[[idx]] <- res
        idx <- idx + 1
      }
    }
  }

  models <- do.call(rbind, model_records)
  models <- models[, c("dataset", "program", "top_n", "model_type", "beta_ASD", "se_ASD", "p_ASD",
                       "n", "n_ASD", "n_Control", "n_common_genes", "formula")]

  model_out <- file.path(outdir, "05_step11H_Gandal_NNLS_adjusted_models.tsv")
  safe_write_tsv(models, model_out)

  # Model comparison summary.
  base <- models[models$model_type == "base", c("dataset", "program", "top_n", "beta_ASD", "p_ASD")]
  colnames(base)[4:5] <- c("base_beta", "base_p")

  pcmod <- models[models$model_type == "NNLS_PC12", c("dataset", "program", "top_n", "beta_ASD", "p_ASD")]
  colnames(pcmod)[4:5] <- c("NNLS_PC12_beta", "NNLS_PC12_p")

  propmod <- models[models$model_type == "NNLS_props", c("dataset", "program", "top_n", "beta_ASD", "p_ASD")]
  colnames(propmod)[4:5] <- c("NNLS_props_beta", "NNLS_props_p")

  comp <- merge(base, pcmod, by = c("dataset", "program", "top_n"), all = TRUE)
  comp <- merge(comp, propmod, by = c("dataset", "program", "top_n"), all = TRUE)

  comp$PC12_direction_preserved <- sign(comp$base_beta) == sign(comp$NNLS_PC12_beta)
  comp$props_direction_preserved <- sign(comp$base_beta) == sign(comp$NNLS_props_beta)

  comp$PC12_beta_ratio <- comp$NNLS_PC12_beta / comp$base_beta
  comp$props_beta_ratio <- comp$NNLS_props_beta / comp$base_beta

  comp_out <- file.path(outdir, "06_step11H_Gandal_NNLS_model_comparison_summary.tsv")
  safe_write_tsv(comp, comp_out)

  # Deconvolution diagnosis associations.
  prop_diag_records <- list()
  j <- 1
  prop_model_dat <- merge(meta[, c("sample_id", "dx", covariates), drop = FALSE], props, by = "sample_id")
  prop_model_dat$dx_binary <- ifelse(prop_model_dat$dx == "ASD", 1, 0)

  for (pcell in prop_cols) {
    f <- as.formula(paste(pcell, "~ dx_binary"))
    fit <- lm(f, data = prop_model_dat)
    sm <- summary(fit)$coefficients
    prop_diag_records[[j]] <- data.frame(
      dataset = dataset,
      proportion = pcell,
      beta_ASD = sm["dx_binary", "Estimate"],
      se_ASD = sm["dx_binary", "Std. Error"],
      p_ASD = sm["dx_binary", "Pr(>|t|)"],
      mean_ASD = mean(prop_model_dat[prop_model_dat$dx == "ASD", pcell], na.rm = TRUE),
      mean_Control = mean(prop_model_dat[prop_model_dat$dx == "Control", pcell], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    j <- j + 1
  }

  prop_diag <- do.call(rbind, prop_diag_records)
  prop_diag$p_fdr <- p.adjust(prop_diag$p_ASD, method = "BH")
  prop_diag_out <- file.path(outdir, "07_step11H_Gandal_NNLS_proportion_diagnosis_models.tsv")
  safe_write_tsv(prop_diag, prop_diag_out)

  # Markdown summary.
  md_out <- file.path(outdir, "08_step11H_overall_summary.md")
  con <- file(md_out, open = "wt")
  writeLines("# Step11H reference-based NNLS deconvolution adjustment", con)
  writeLines("", con)

  writeLines("## Purpose", con)
  writeLines("", con)
  writeLines("This analysis uses PsychENCODE broad cell-type pseudobulk reference profiles to estimate broad cell-type proportions in Gandal 2022 bulk cortex samples by constrained NNLS-style deconvolution, then re-tests NTM score diagnosis effects before and after deconvolution adjustment.", con)
  writeLines("", con)

  writeLines("## Reference construction", con)
  writeLines("", con)
  writeLines(sprintf("- Reference object: `%s`", ref_rds), con)
  writeLines(sprintf("- Maximum cells sampled per broad cell type: %s", max_cells_per_type), con)
  writeLines(sprintf("- Broad cell types in reference: %s", paste(colnames(ref_avg), collapse = ", ")), con)
  writeLines(sprintf("- Reference genes after filtering: %s", nrow(ref_avg)), con)
  writeLines(sprintf("- Marker genes used for NNLS: %s unique genes", length(unique(markers$gene))), con)
  writeLines("", con)

  writeLines("## Bulk cohort", con)
  writeLines("", con)
  writeLines(sprintf("- Dataset: %s", dataset), con)
  writeLines(sprintf("- ASD samples: %s", sum(meta$dx == "ASD")), con)
  writeLines(sprintf("- Control samples: %s", sum(meta$dx == "Control")), con)
  writeLines(sprintf("- Covariates detected: %s", ifelse(length(covariates) == 0, "none", paste(covariates, collapse = ", "))), con)
  writeLines("", con)

  writeLines("## NTM diagnosis models before and after NNLS adjustment", con)
  writeLines("", con)
  writeLines(knitr::kable(comp, format = "markdown", digits = 4), con)
  writeLines("", con)

  writeLines("## Estimated cell-type proportion diagnosis associations", con)
  writeLines("", con)
  writeLines(knitr::kable(prop_diag, format = "markdown", digits = 4), con)
  writeLines("", con)

  writeLines("## Interpretation guide", con)
  writeLines("", con)
  writeLines("- If NNLS_PC12_direction_preserved and/or props_direction_preserved remain TRUE for most NTM programs, this supports that NTM effects are not fully explained by reference-estimated broad cell-type proportions.", con)
  writeLines("- If beta ratios are attenuated, the correct interpretation is composition-coupled but not fully composition-explained, consistent with Step10E.", con)
  writeLines("- If some modules reverse after all-proportion adjustment, report this as a sensitivity boundary rather than a failed primary result, because broad-proportion adjustment is highly collinear and conservative.", con)
  close(con)

  msg("Finished Step11H")
  msg("Wrote:", prop_out)
  msg("Wrote:", score_out)
  msg("Wrote:", model_out)
  msg("Wrote:", comp_out)
  msg("Wrote:", prop_diag_out)
  msg("Wrote:", md_out)
}

main()
