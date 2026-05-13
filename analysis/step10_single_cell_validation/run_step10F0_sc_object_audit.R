#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
STEP <- file.path(BASE, "neurotrace_algorithm_project/step10_single_cell_validation")
OUT <- file.path(STEP, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

OBJECTS <- data.table(
  dataset = c("PsychENCODE", "Velmeshev"),
  rds_path = c(
    "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object.rds",
    "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/Velmeshev_Object.rds"
  )
)

safe_fwrite <- function(x, file) {
  fwrite(as.data.table(x), file, sep = "\t", quote = FALSE, na = "NA")
}

detect_cols <- function(md) {
  cn <- names(md)
  low <- tolower(cn)

  list(
    diagnosis_candidates = cn[grepl("diagnosis|dx|condition|group|disease|asd", low)],
    donor_candidates = cn[grepl("donor|subject|individual|sample|case|patient|brain", low)],
    celltype_candidates = cn[grepl("cell.?type|celltype|annotation|annot|cluster|class|subclass|lineage|major|minor|broad", low)],
    region_candidates = cn[grepl("region|area|cortex|brain", low)],
    batch_candidates = cn[grepl("batch|seq|library|run|platform", low)]
  )
}

summarize_col <- function(md, col, dataset, max_levels = 40) {
  if (!col %in% names(md)) return(data.table())
  tab <- as.data.table(table(as.character(md[[col]]), useNA = "ifany"))
  setnames(tab, c("value", "n"))
  tab[, dataset := dataset]
  tab[, column := col]
  setcolorder(tab, c("dataset", "column", "value", "n"))
  tab[order(-n)][seq_len(min(.N, max_levels))]
}

message("[", timestamp(), "] Step10F0 single-cell object audit started")

object_rows <- list()
meta_col_rows <- list()
candidate_rows <- list()
preview_rows <- list()
level_rows <- list()

for (i in seq_len(nrow(OBJECTS))) {
  dataset <- OBJECTS$dataset[i]
  rds_path <- OBJECTS$rds_path[i]

  message("[", timestamp(), "] Loading ", dataset, ": ", rds_path)
  if (!file.exists(rds_path)) stop("Missing RDS: ", rds_path)

  obj <- readRDS(rds_path)

  cls <- paste(class(obj), collapse = ";")
  assays <- NA_character_
  reductions <- NA_character_
  n_cells <- NA_integer_
  n_features_default <- NA_integer_
  default_assay <- NA_character_

  if ("Seurat" %in% class(obj)) {
    assays <- paste(names(obj@assays), collapse = ";")
    reductions <- paste(names(obj@reductions), collapse = ";")
    n_cells <- ncol(obj)
    default_assay <- tryCatch(DefaultAssay(obj), error = function(e) NA_character_)
    n_features_default <- tryCatch(nrow(obj[[default_assay]]), error = function(e) NA_integer_)
    md <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  } else {
    md <- tryCatch(as.data.table(obj@meta.data, keep.rownames = "cell_id"), error = function(e) data.table())
    if (nrow(md) == 0) {
      md <- data.table(note = "Object is not Seurat or lacks meta.data")
    }
  }

  object_rows[[length(object_rows) + 1]] <- data.table(
    dataset = dataset,
    rds_path = rds_path,
    class = cls,
    assays = assays,
    default_assay = default_assay,
    reductions = reductions,
    n_cells = n_cells,
    n_features_default_assay = n_features_default,
    n_metadata_cols = ncol(md)
  )

  meta_col_rows[[length(meta_col_rows) + 1]] <- data.table(
    dataset = dataset,
    column = names(md),
    class = sapply(md, function(x) paste(class(x), collapse = ";")),
    n_unique = sapply(md, function(x) uniqueN(as.character(x), na.rm = FALSE)),
    n_missing = sapply(md, function(x) sum(is.na(x)))
  )

  cand <- detect_cols(md)
  for (nm in names(cand)) {
    candidate_rows[[length(candidate_rows) + 1]] <- data.table(
      dataset = dataset,
      candidate_type = nm,
      columns = paste(cand[[nm]], collapse = ";")
    )
  }

  # preview first 20 metadata columns
  pcols <- names(md)[seq_len(min(20, ncol(md)))]
  prev <- md[seq_len(min(5, nrow(md))), ..pcols]
  prev[, dataset := dataset]
  preview_rows[[length(preview_rows) + 1]] <- prev

  # level summaries for candidate columns
  candidate_cols <- unique(unlist(cand))
  candidate_cols <- candidate_cols[candidate_cols %in% names(md)]
  for (cc in candidate_cols) {
    level_rows[[length(level_rows) + 1]] <- summarize_col(md, cc, dataset)
  }
}

object_summary <- rbindlist(object_rows, fill = TRUE)
metadata_columns <- rbindlist(meta_col_rows, fill = TRUE)
candidate_summary <- rbindlist(candidate_rows, fill = TRUE)
metadata_preview <- rbindlist(preview_rows, fill = TRUE)
level_summary <- rbindlist(level_rows, fill = TRUE)

safe_fwrite(object_summary, file.path(OUT, "01_step10F0_sc_object_summary.tsv"))
safe_fwrite(metadata_columns, file.path(OUT, "02_step10F0_sc_metadata_columns.tsv"))
safe_fwrite(candidate_summary, file.path(OUT, "03_step10F0_sc_candidate_columns.tsv"))
safe_fwrite(metadata_preview, file.path(OUT, "04_step10F0_sc_metadata_preview.tsv"))
safe_fwrite(level_summary, file.path(OUT, "05_step10F0_sc_candidate_level_summary.tsv"))

summary_lines <- c(
  "# NeuroTRACE Step10F0 single-cell object audit summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Audit PsychENCODE and Velmeshev single-cell/single-nucleus objects before NTM scoring and donor-level validation.",
  "",
  "## Object summary",
  paste(capture.output(print(object_summary)), collapse = "\n"),
  "",
  "## Candidate columns",
  paste(capture.output(print(candidate_summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "Use this audit to choose diagnosis, donor and cell-type columns for Step10F single-cell NTM localization."
)
writeLines(summary_lines, file.path(OUT, "06_step10F0_sc_object_audit_summary.md"))

message("[", timestamp(), "] Step10F0 done")
message("[", timestamp(), "] Results: ", OUT)
