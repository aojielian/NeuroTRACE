#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
STEP <- file.path(BASE, "neurotrace_algorithm_project/step10_single_cell_validation")
OUT <- file.path(STEP, "results")
FIG <- file.path(STEP, "figures")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

MODULE_FILE <- file.path(BASE, "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")

OBJECTS <- data.table(
  dataset = c("PsychENCODE", "Velmeshev"),
  rds_path = c(
    "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object.rds",
    "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/Velmeshev_Object.rds"
  ),
  diagnosis_col = c("Diagnosis", "diagnosis"),
  donor_col = c("individual_ID", "individual"),
  celltype_col = c("annotation", "cluster"),
  age_col = c("Age", "age"),
  sex_col = c("Sex_Chromosome", "sex"),
  pmi_col = c("PMI", "post.mortem.interval..hours."),
  rin_col = c("RIN", "RNA.Integrity.Number")
)
RUN_DATASET <- Sys.getenv("NEUROTRACE_SC_DATASET", "ALL")
if (RUN_DATASET != "ALL") {
  OBJECTS <- OBJECTS[dataset == RUN_DATASET]
  if (nrow(OBJECTS) == 0) stop("No matching dataset in OBJECTS for NEUROTRACE_SC_DATASET=", RUN_DATASET)
}

out_file <- function(filename) {
  if (RUN_DATASET == "ALL") {
    file.path(OUT, filename)
  } else {
    file.path(OUT, paste0(RUN_DATASET, "_", filename))
  }
}

fig_file <- function(filename) {
  if (RUN_DATASET == "ALL") {
    file.path(FIG, filename)
  } else {
    file.path(FIG, paste0(RUN_DATASET, "_", filename))
  }
}

MAIN_TOP_N <- c(200, 500)
MIN_CELLS_PER_DONOR_CELLTYPE <- 20

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
  out[grepl("ctl|control|normal", y)] <- "Control"
  out[grepl("asd|autism", y)] <- "ASD"
  out
}

map_velmeshev_broad <- function(x) {
  y <- as.character(x)
  low <- tolower(y)
  out <- rep("OTHER", length(y))
  out[grepl("oligodendro", low)] <- "ODC"
  out[grepl("^opc$|opc", low)] <- "OPC"
  out[grepl("astro", low)] <- "AST"
  out[grepl("micro|mg", low)] <- "MG"
  out[grepl("endo|endothel", low)] <- "END"
  out[grepl("l2|l3|l4|l5|l6|neu-nrgn|excit|projection", low)] <- "EXN"
  out[grepl("interneuron|inh|gad|sst|pvalb|vip", low)] <- "INN"
  out
}

map_psychencode_broad <- function(x) {
  y <- as.character(x)
  valid <- c("EXN", "INN", "AST", "ODC", "OPC", "MG", "END")
  out <- toupper(y)
  out[!out %in% valid] <- "OTHER"
  out
}

get_counts_or_data <- function(obj) {
  assay <- tryCatch(DefaultAssay(obj), error = function(e) "RNA")
  assay_obj <- obj[[assay]]

  mat <- NULL

  # Seurat v5 layer-compatible
  mat <- tryCatch(GetAssayData(obj, assay = assay, layer = "data"), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0) {
    mat <- tryCatch(GetAssayData(obj, assay = assay, slot = "data"), error = function(e) NULL)
  }
  if (is.null(mat) || nrow(mat) == 0) {
    mat <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  }
  if (is.null(mat) || nrow(mat) == 0) {
    mat <- tryCatch(GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  }
  if (is.null(mat) || nrow(mat) == 0) {
    stop("Could not extract expression matrix from Seurat object.")
  }

  mat
}


row_zscore_matrix <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mu <- rowMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  out <- sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
  out[!is.finite(out)] <- 0
  out
}

score_module_on_matrix <- function(mat, module_sub) {
  genes <- clean_gene(module_sub$gene_symbol_fixed)
  weights <- as.numeric(module_sub$weight)
  names(weights) <- genes
  common <- intersect(clean_gene(rownames(mat)), genes)

  if (length(common) < 5) {
    return(list(score = rep(NA_real_, ncol(mat)), n_common = length(common), common = common))
  }

  # map rownames
  rn_clean <- clean_gene(rownames(mat))
  idx <- match(common, rn_clean)
  w <- weights[common]

  score <- as.numeric(crossprod(w, mat[idx, , drop = FALSE]) / (sum(abs(w)) + 1e-12))
  names(score) <- colnames(mat)
  list(score = score, n_common = length(common), common = common)
}

fit_donor_model <- function(dt, score_col) {
  d <- copy(dt)
  d <- d[dx %in% c("ASD", "Control")]
  d[, dx_bin := ifelse(dx == "ASD", 1, 0)]

  covars <- c()
  for (cc in c("age_cov", "pmi_cov", "rin_cov")) {
    if (cc %in% names(d)) {
      vals <- suppressWarnings(as.numeric(d[[cc]]))
      if (sum(is.finite(vals)) >= 5 && sd(vals, na.rm = TRUE) > 0) {
        med <- median(vals, na.rm = TRUE)
        vals[!is.finite(vals)] <- med
        d[, (cc) := as.numeric(scale(vals))]
        covars <- c(covars, cc)
      }
    }
  }

  if ("sex_cov" %in% names(d)) {
    d[is.na(sex_cov) | sex_cov == "", sex_cov := "Unknown"]
    if (length(unique(d$sex_cov)) >= 2) covars <- c(covars, "sex_cov")
  }

  rhs <- c("dx_bin", covars)
  form <- as.formula(paste(score_col, "~", paste(rhs, collapse = " + ")))

  if (length(unique(d$dx)) < 2 || nrow(d) < length(rhs) + 3) {
    return(data.table(beta = NA_real_, se = NA_real_, t = NA_real_, p_value = NA_real_, n_donors = nrow(d), n_ASD = sum(d$dx == "ASD"), n_Control = sum(d$dx == "Control"), covariates = paste(covars, collapse = ";")))
  }

  fit <- tryCatch(lm(form, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(beta = NA_real_, se = NA_real_, t = NA_real_, p_value = NA_real_, n_donors = nrow(d), n_ASD = sum(d$dx == "ASD"), n_Control = sum(d$dx == "Control"), covariates = paste(covars, collapse = ";")))
  }

  sm <- summary(fit)$coefficients
  if (!"dx_bin" %in% rownames(sm)) {
    return(data.table(beta = NA_real_, se = NA_real_, t = NA_real_, p_value = NA_real_, n_donors = nrow(d), n_ASD = sum(d$dx == "ASD"), n_Control = sum(d$dx == "Control"), covariates = paste(covars, collapse = ";")))
  }

  data.table(
    beta = sm["dx_bin", "Estimate"],
    se = sm["dx_bin", "Std. Error"],
    t = sm["dx_bin", "t value"],
    p_value = sm["dx_bin", "Pr(>|t|)"],
    n_donors = nrow(d),
    n_ASD = sum(d$dx == "ASD"),
    n_Control = sum(d$dx == "Control"),
    covariates = paste(covars, collapse = ";")
  )
}

aggregate_by_group <- function(mat, groups) {
  groups <- as.factor(groups)
  mm <- sparse.model.matrix(~ 0 + groups)
  colnames(mm) <- levels(groups)
  sums <- mat %*% mm
  counts <- as.numeric(table(groups))
  avg <- t(t(sums) / counts)
  list(avg = avg, counts = counts, group_levels = levels(groups))
}

message("[", timestamp(), "] Step10F single-cell NTM localization started")

if (!file.exists(MODULE_FILE)) stop("Missing module file: ", MODULE_FILE)

modules <- fread(MODULE_FILE)
modules <- modules[top_n %in% MAIN_TOP_N]
modules[, gene_key := clean_gene(gene_symbol_fixed)]

object_summary_rows <- list()
group_meta_rows <- list()
score_rows <- list()
gene_match_rows <- list()
model_rows <- list()
celltype_summary_rows <- list()

for (i in seq_len(nrow(OBJECTS))) {
  ds <- OBJECTS$dataset[i]
  rds <- OBJECTS$rds_path[i]
  message("[", timestamp(), "] Loading object: ", ds)

  if (!file.exists(rds)) stop("Missing RDS: ", rds)

  obj <- readRDS(rds)
  message("[", timestamp(), "] Loaded object: ", ds,
          " cells=", ncol(obj), " genes=", nrow(obj))
  flush.console()

  md <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  message("[", timestamp(), "] Metadata loaded: ", ds,
          " rows=", nrow(md), " cols=", ncol(md))
  flush.console()

  dx_col <- OBJECTS$diagnosis_col[i]
  donor_col <- OBJECTS$donor_col[i]
  cell_col <- OBJECTS$celltype_col[i]
  age_col <- OBJECTS$age_col[i]
  sex_col <- OBJECTS$sex_col[i]
  pmi_col <- OBJECTS$pmi_col[i]
  rin_col <- OBJECTS$rin_col[i]

  md[, dx := standardize_dx(get(dx_col))]
  md[, donor_id := as.character(get(donor_col))]
  md[, raw_celltype := as.character(get(cell_col))]

  if (ds == "PsychENCODE") {
    md[, broad_celltype := map_psychencode_broad(raw_celltype)]
  } else {
    md[, broad_celltype := map_velmeshev_broad(raw_celltype)]
  }

  if (age_col %in% names(md)) suppressWarnings(md[, age_cov := as.numeric(get(age_col))])
  if (pmi_col %in% names(md)) suppressWarnings(md[, pmi_cov := as.numeric(get(pmi_col))])
  if (rin_col %in% names(md)) suppressWarnings(md[, rin_cov := as.numeric(get(rin_col))])
  if (sex_col %in% names(md)) md[, sex_cov := as.factor(as.character(get(sex_col)))]

  keep_cells <- md[dx %in% c("ASD", "Control") & broad_celltype %in% c("EXN", "INN", "AST", "ODC", "OPC", "MG", "END") & !is.na(donor_id), cell_id]

  message("[", timestamp(), "] Cells retained after dx/celltype filtering: ", ds,
          " n=", length(keep_cells))
  flush.console()

  obj <- subset(obj, cells = keep_cells)
  md <- md[cell_id %in% keep_cells]
  md <- md[match(colnames(obj), cell_id)]

  message("[", timestamp(), "] Extracting expression matrix: ", ds)
  flush.console()
  mat <- get_counts_or_data(obj)
  message("[", timestamp(), "] Expression matrix extracted: ", ds,
          " genes=", nrow(mat), " cells=", ncol(mat))
  flush.console()

  # Ensure gene symbols uppercase
  rownames(mat) <- clean_gene(rownames(mat))
  keep_gene <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep_gene, , drop = FALSE]
  message("[", timestamp(), "] Gene filtering complete: ", ds,
          " genes=", nrow(mat), " cells=", ncol(mat))
  flush.console()

  group_id <- paste(md$donor_id, md$broad_celltype, sep = "||")
  message("[", timestamp(), "] Aggregating donor x celltype pseudobulk: ", ds,
          " groups=", length(unique(group_id)))
  flush.console()
  agg <- aggregate_by_group(mat, group_id)
  avg <- agg$avg
  message("[", timestamp(), "] Pseudobulk aggregation complete: ", ds,
          " genes=", nrow(avg), " groups=", ncol(avg))
  flush.console()

  message("[", timestamp(), "] Applying gene-wise z-score to pseudobulk matrix: ", ds)
  flush.console()
  avg_z <- row_zscore_matrix(avg)
  message("[", timestamp(), "] Gene-wise z-score complete: ", ds)
  flush.console()

  group_dt <- data.table(group_id = colnames(avg))
  group_dt[, donor_id := sub("\\|\\|.*$", "", group_id)]
  group_dt[, broad_celltype := sub("^.*\\|\\|", "", group_id)]
  group_dt[, n_cells := as.integer(agg$counts[match(group_id, agg$group_levels)])]

  # donor covariates by majority/mean
  donor_meta <- md[, .(
    dx = names(sort(table(dx), decreasing = TRUE))[1],
    age_cov = suppressWarnings(mean(as.numeric(age_cov), na.rm = TRUE)),
    pmi_cov = suppressWarnings(mean(as.numeric(pmi_cov), na.rm = TRUE)),
    rin_cov = suppressWarnings(mean(as.numeric(rin_cov), na.rm = TRUE)),
    sex_cov = if ("sex_cov" %in% names(.SD)) names(sort(table(as.character(sex_cov)), decreasing = TRUE))[1] else NA_character_
  ), by = donor_id]

  group_dt <- merge(group_dt, donor_meta, by = "donor_id", all.x = TRUE)

  object_summary_rows[[length(object_summary_rows) + 1]] <- data.table(
    dataset = ds,
    n_cells_used = ncol(mat),
    n_genes = nrow(mat),
    n_donors = uniqueN(md$donor_id),
    n_groups = ncol(avg),
    n_ASD_cells = sum(md$dx == "ASD"),
    n_Control_cells = sum(md$dx == "Control")
  )

  celltype_summary_tmp <- md[, .(
    n_cells = .N,
    n_donors = uniqueN(donor_id),
    n_ASD_donors = uniqueN(donor_id[dx == "ASD"]),
    n_Control_donors = uniqueN(donor_id[dx == "Control"])
  ), by = broad_celltype]
  celltype_summary_tmp[, dataset := ds]
  setcolorder(celltype_summary_tmp, c("dataset", "broad_celltype", "n_cells", "n_donors", "n_ASD_donors", "n_Control_donors"))
  celltype_summary_rows[[length(celltype_summary_rows) + 1]] <- celltype_summary_tmp

  # Score modules on pseudobulk group matrix
  message("[", timestamp(), "] Starting NTM scoring on pseudobulk matrix: ", ds)
  flush.console()

  for (prog in unique(modules$program)) {
    for (N in MAIN_TOP_N) {
      message("[", timestamp(), "]   Scoring ", ds, " ", prog, " top", N)
      flush.console()
      msub <- modules[get("program") == prog & top_n == N]
      setorder(msub, rank_within_program)
      msub <- msub[!duplicated(gene_key)]

      sc <- score_module_on_matrix(avg_z, msub)

      score_col <- paste0(prog, "_top", N)
      group_dt[, (score_col) := as.numeric(sc$score[match(group_id, names(sc$score))])]

      gene_match_rows[[length(gene_match_rows) + 1]] <- data.table(
        dataset = ds,
        program = prog,
        top_n = N,
        n_module_genes = uniqueN(msub$gene_key),
        n_common_genes = sc$n_common,
        overlap_rate = sc$n_common / uniqueN(msub$gene_key),
        common_genes = paste(sc$common, collapse = ";")
      )
    }
  }

  group_meta_rows[[length(group_meta_rows) + 1]] <- group_dt

  score_long <- melt(
    group_dt,
    id.vars = c("donor_id", "broad_celltype", "n_cells", "dx", "age_cov", "pmi_cov", "rin_cov", "sex_cov"),
    measure.vars = grep("^NTM", names(group_dt), value = TRUE),
    variable.name = "ntm_score_name",
    value.name = "ntm_score"
  )
  score_long[, dataset := ds]
  score_rows[[length(score_rows) + 1]] <- score_long

  # donor-level ASD vs Control within each cell type and NTM score
  message("[", timestamp(), "] Starting donor-level models: ", ds)
  flush.console()

  for (ct in sort(unique(group_dt$broad_celltype))) {
    message("[", timestamp(), "]   Modeling cell type: ", ds, " ", ct)
    flush.console()
    dt_ct <- group_dt[broad_celltype == ct & n_cells >= MIN_CELLS_PER_DONOR_CELLTYPE]
    if (nrow(dt_ct) == 0) next

    for (score_col in grep("^NTM", names(dt_ct), value = TRUE)) {
      fit <- fit_donor_model(dt_ct, score_col)

      # parse program/top_n
      prog <- sub("_top[0-9]+$", "", score_col)
      top_n <- as.integer(sub("^.*_top", "", score_col))

      model_rows[[length(model_rows) + 1]] <- cbind(
        data.table(dataset = ds, broad_celltype = ct, program = prog, top_n = top_n, ntm_score_name = score_col),
        fit
      )
    }
  }

  message("[", timestamp(), "] Finished dataset: ", ds)
  flush.console()
}

object_summary <- rbindlist(object_summary_rows, fill = TRUE)
celltype_summary <- rbindlist(celltype_summary_rows, fill = TRUE)
group_meta <- rbindlist(group_meta_rows, fill = TRUE)
score_long <- rbindlist(score_rows, fill = TRUE)
gene_match <- rbindlist(gene_match_rows, fill = TRUE)
models <- rbindlist(model_rows, fill = TRUE)

models[, fdr := p.adjust(p_value, method = "BH")]
models[, direction_positive := beta > 0]

# localization by mean score across diagnosis groups and cell types
localization <- score_long[, .(
  mean_score = mean(ntm_score, na.rm = TRUE),
  median_score = median(ntm_score, na.rm = TRUE),
  n_groups = .N,
  n_donors = uniqueN(donor_id),
  n_ASD_donors = uniqueN(donor_id[dx == "ASD"]),
  n_Control_donors = uniqueN(donor_id[dx == "Control"])
), by = .(dataset, broad_celltype, ntm_score_name)]

# model summary
model_summary <- models[, .(
  n_tests = .N,
  n_positive = sum(direction_positive, na.rm = TRUE),
  n_p_lt_0.05 = sum(p_value < 0.05, na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  min_fdr = min(fdr, na.rm = TRUE),
  median_beta = median(beta, na.rm = TRUE)
), by = .(dataset, program, top_n)]

safe_fwrite(object_summary, out_file("07_step10F_sc_object_used_summary.tsv"))
safe_fwrite(celltype_summary, out_file("08_step10F_sc_celltype_donor_summary.tsv"))
safe_fwrite(group_meta, out_file("09_step10F_sc_donor_celltype_NTM_scores_wide.tsv.gz"))
safe_fwrite(score_long, out_file("10_step10F_sc_donor_celltype_NTM_scores_long.tsv.gz"))
safe_fwrite(gene_match, out_file("11_step10F_sc_NTM_gene_match.tsv"))
safe_fwrite(localization, out_file("12_step10F_sc_NTM_celltype_localization_summary.tsv"))
safe_fwrite(models, out_file("13_step10F_sc_donor_level_ASD_models.tsv"))
safe_fwrite(model_summary, out_file("14_step10F_sc_donor_level_model_summary.tsv"))

# plots
try({
  pdf(fig_file("02_step10F_sc_NTM_celltype_localization_heatmap.pdf"), width = 11, height = 7)
  par(mfrow = c(2, 3), mar = c(8, 4, 3, 1))
  for (ds in unique(localization$dataset)) {
    for (score_name in unique(localization$ntm_score_name)) {
      ss <- localization[dataset == ds & ntm_score_name == score_name]
      ss <- ss[order(broad_celltype)]
      barplot(ss$mean_score, names.arg = ss$broad_celltype, las = 2,
              main = paste(ds, score_name), ylab = "Mean donor-celltype score")
    }
  }
  dev.off()
}, silent = TRUE)

summary_lines <- c(
  "# NeuroTRACE Step10F single-cell NTM localization summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Localize NeuroTRACE-native NTM modules across adult ASD single-cell/single-nucleus cell types using gene-wise standardized donor-celltype pseudobulk matrices and test donor-level ASD versus Control effects within broad cell types.",
  "",
  paste0("## Run dataset filter\n", RUN_DATASET),
  "",
  "## Objects used",
  paste(capture.output(print(object_summary)), collapse = "\n"),
  "",
  "## Cell-type donor summary",
  paste(capture.output(print(celltype_summary)), collapse = "\n"),
  "",
  "## Gene matching",
  paste(capture.output(print(gene_match[, .(n_common_genes = min(n_common_genes), median_overlap = median(overlap_rate)), by = .(dataset, program, top_n)])), collapse = "\n"),
  "",
  "## Donor-level model summary",
  paste(capture.output(print(model_summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "Use localization summaries to identify the cell-type context of NTM modules. Donor-level models avoid treating cells as independent observations. These analyses are supporting validation and cell-type contextualization, not a new discovery branch."
)
writeLines(summary_lines, out_file("15_step10F_sc_NTM_localization_summary.md"))

message("[", timestamp(), "] Step10F done")
message("[", timestamp(), "] Results: ", OUT)
