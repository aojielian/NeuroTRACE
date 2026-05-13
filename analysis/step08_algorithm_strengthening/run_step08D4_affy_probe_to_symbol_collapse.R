#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
PROJECT <- file.path(BASE, "neurotrace_algorithm_project")
STEP <- file.path(PROJECT, "step08_algorithm_strengthening")
OUT <- file.path(STEP, "results")
PARSED <- file.path(PROJECT, "cross_disease_GEO/parsed")
PLATFORM_DIR <- file.path(PROJECT, "cross_disease_GEO/platforms")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DATASETS <- data.table(
  dataset = c("GSE53987", "GSE12649", "GSE21138"),
  platform = c("GPL570", "GPL96", "GPL570"),
  expression_probe_path = file.path(PARSED, c(
    "GSE53987/GSE53987.expression_probe.tsv.gz",
    "GSE12649/GSE12649.expression_probe.tsv.gz",
    "GSE21138/GSE21138.expression_probe.tsv.gz"
  )),
  sample_path = file.path(PARSED, c(
    "GSE53987/GSE53987.samples.tsv",
    "GSE12649/GSE12649.samples.tsv",
    "GSE21138/GSE21138.samples.tsv"
  ))
)

safe_fwrite <- function(x, file) {
  x <- as.data.table(x)
  if (grepl(".gz", file, fixed = TRUE) && endsWith(file, ".gz")) {
    tmp <- sub(".gz$", "", file)
    fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
    system(sprintf("gzip -f %s", shQuote(tmp)))
  } else {
    fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
  }
}

clean_gene <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %in% c("", "NA", "---", "NULL", "NAN")] <- NA_character_
  y
}

ends_any <- function(x, suffixes) {
  out <- rep(FALSE, length(x))
  for (suf in suffixes) out <- out | endsWith(x, suf)
  out
}

locate_platform_file <- function(platform) {
  if (!dir.exists(PLATFORM_DIR)) {
    stop("Platform directory does not exist: ", PLATFORM_DIR)
  }

  files <- list.files(PLATFORM_DIR, full.names = TRUE)
  if (length(files) == 0) {
    stop("No files found in platform directory: ", PLATFORM_DIR)
  }

  bn <- basename(files)

  keep <- startsWith(bn, platform) & ends_any(
    bn,
    c(".txt", ".txt.gz", ".annot", ".annot.gz", ".soft", ".soft.gz")
  )

  cand <- files[keep]

  if (length(cand) == 0) {
    stop(
      "No local platform annotation file found for ", platform, ".\n",
      "Expected files such as ", platform, "-55999.txt or ", platform, "-57554.txt in: ",
      PLATFORM_DIR, "\nCurrent files:\n", paste(bn, collapse = "\n")
    )
  }

  # Prefer manually downloaded GEO table txt files such as GPL570-55999.txt
  cbn <- basename(cand)
  priority <- rep(9L, length(cand))
  priority[endsWith(cbn, ".txt") & grepl("-", cbn, fixed = TRUE)] <- 1L
  priority[endsWith(cbn, ".txt.gz") & grepl("-", cbn, fixed = TRUE)] <- 2L
  priority[endsWith(cbn, ".annot")] <- 3L
  priority[endsWith(cbn, ".annot.gz")] <- 4L
  priority[endsWith(cbn, ".soft")] <- 5L
  priority[endsWith(cbn, ".soft.gz")] <- 6L

  cand <- cand[order(priority, cbn)]
  cand[1]
}

read_platform_lines <- function(file) {
  if (endsWith(file, ".gz")) {
    con <- gzfile(file, "rt")
  } else {
    con <- file(file, "rt")
  }
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

read_platform_table <- function(platform) {
  f <- locate_platform_file(platform)
  message("[", timestamp(), "] Using platform file for ", platform, ": ", f)

  lines <- read_platform_lines(f)

  # GEO platform txt usually has comment lines beginning with #, then a real header beginning with ID
  header_i <- which(startsWith(lines, "ID\t"))
  if (length(header_i) == 0) {
    # SOFT platform table case
    b <- which(startsWith(lines, "!platform_table_begin"))
    e <- which(startsWith(lines, "!platform_table_end"))
    if (length(b) > 0 && length(e) > 0 && e[1] > b[1] + 1) {
      header_i <- b[1] + 1
      end_i <- e[1] - 1
      tab_lines <- lines[header_i:end_i]
    } else {
      stop("Could not find platform table header beginning with ID in file: ", f)
    }
  } else {
    header_i <- header_i[1]
    tab_lines <- lines[header_i:length(lines)]
  }

  dt <- fread(
    text = paste(tab_lines, collapse = "\n"),
    sep = "\t",
    quote = "",
    fill = TRUE,
    showProgress = FALSE
  )

  if (nrow(dt) == 0 || ncol(dt) < 2) {
    stop("Platform table parsed but appears empty: ", f)
  }

  cn <- names(dt)
  low <- tolower(cn)

  probe_col <- "ID"
  if (!probe_col %in% cn) probe_col <- cn[1]

  sym_idx <- which(low == "gene symbol")
  if (length(sym_idx) == 0) sym_idx <- which(low == "gene_symbol")
  if (length(sym_idx) == 0) sym_idx <- which(low == "symbol")
  if (length(sym_idx) == 0) sym_idx <- grep("symbol", low, fixed = TRUE)

  if (length(sym_idx) == 0) {
    stop("Could not detect gene symbol column in platform file. Columns: ", paste(cn, collapse = " | "))
  }

  symbol_col <- cn[sym_idx[1]]

  entrez_idx <- which(low == "entrez_gene_id")
  if (length(entrez_idx) == 0) entrez_idx <- which(low == "entrez gene")
  if (length(entrez_idx) == 0) entrez_idx <- grep("entrez", low, fixed = TRUE)
  entrez_col <- if (length(entrez_idx) > 0) cn[entrez_idx[1]] else NA_character_

  mp <- data.table(
    probe_id = as.character(dt[[probe_col]]),
    gene_symbol_raw = as.character(dt[[symbol_col]])
  )

  if (!is.na(entrez_col)) {
    mp[, entrez_id := as.character(dt[[entrez_col]])]
  } else {
    mp[, entrez_id := NA_character_]
  }

  mp <- mp[
    !is.na(probe_id) &
      probe_id != "" &
      !is.na(gene_symbol_raw) &
      gene_symbol_raw != ""
  ]

  # Split multi-symbol annotations
  mp[, gene_symbol_raw := gsub(" // ", " /// ", gene_symbol_raw, fixed = TRUE)]
  mp[, gene_symbol_raw := gsub("///", " /// ", gene_symbol_raw, fixed = TRUE)]

  mp_long <- mp[, .(
    gene_symbol = unlist(strsplit(gene_symbol_raw, " /// ", fixed = TRUE))
  ), by = .(probe_id, entrez_id)]

  mp_long[, gene_symbol := clean_gene(gene_symbol)]
  mp_long <- mp_long[!is.na(gene_symbol)]
  mp_long[, platform := platform]
  mp_long <- unique(mp_long, by = c("probe_id", "gene_symbol"))

  if (nrow(mp_long) == 0) {
    stop("No valid probe-to-symbol rows after parsing platform file: ", f)
  }

  mp_long
}

standardize_dx_region <- function(dataset, md) {
  md <- as.data.table(md)

  dx <- rep(NA_character_, nrow(md))

  # GSE53987 has disease_state
  if ("disease_state" %in% names(md)) {
    y <- tolower(as.character(md$disease_state))
    dx[grepl("control", y, fixed = TRUE)] <- "Control"
    dx[grepl("schizophrenia", y, fixed = TRUE)] <- "SCZ"
    dx[grepl("bipolar", y, fixed = TRUE)] <- "BD"
    dx[grepl("depressive", y, fixed = TRUE)] <- "MDD"
    dx[grepl("depression", y, fixed = TRUE)] <- "MDD"
  }

  # GSE12649 uses source_name_ch1 / characteristics_ch1_2
  if (all(is.na(dx)) && "source_name_ch1" %in% names(md)) {
    y <- tolower(as.character(md$source_name_ch1))
    dx[grepl("control", y, fixed = TRUE)] <- "Control"
    dx[grepl("schizophrenia", y, fixed = TRUE)] <- "SCZ"
    dx[grepl("bipolar", y, fixed = TRUE)] <- "BD"
    dx[grepl("depressive", y, fixed = TRUE)] <- "MDD"
    dx[grepl("depression", y, fixed = TRUE)] <- "MDD"
  }

  if (all(is.na(dx)) && "characteristics_ch1_2" %in% names(md)) {
    y <- tolower(as.character(md$characteristics_ch1_2))
    dx[grepl("control", y, fixed = TRUE)] <- "Control"
    dx[grepl("schizophrenia", y, fixed = TRUE)] <- "SCZ"
    dx[grepl("bipolar", y, fixed = TRUE)] <- "BD"
    dx[grepl("depressive", y, fixed = TRUE)] <- "MDD"
    dx[grepl("depression", y, fixed = TRUE)] <- "MDD"
  }

  # GSE21138 title has Control / illness labels
  if (all(is.na(dx)) && "title" %in% names(md)) {
    y <- tolower(as.character(md$title))
    dx[grepl("control", y, fixed = TRUE)] <- "Control"
    dx[grepl("schiz", y, fixed = TRUE)] <- "SCZ"
  }

  # GSE21138 stage_of_illness column may contain control or patient stages
  illness_cols <- names(md)[grepl("stage_of_illness", names(md), fixed = TRUE)]
  if (any(is.na(dx)) && length(illness_cols) > 0) {
    y <- tolower(as.character(md[[illness_cols[1]]]))
    dx[is.na(dx) & grepl("control", y, fixed = TRUE)] <- "Control"
    dx[is.na(dx) & !grepl("control", y, fixed = TRUE) & y != "" & y != "nan"] <- "SCZ"
  }

  md[, dx_multi := dx]

  region <- rep(NA_character_, nrow(md))

  if ("tissue" %in% names(md)) {
    region <- as.character(md$tissue)
  }

  if ("source_name_ch1" %in% names(md)) {
    y <- tolower(as.character(md$source_name_ch1))
    region[grepl("pre-frontal", y, fixed = TRUE)] <- "PFC_BA46"
    region[grepl("prefrontal", y, fixed = TRUE)] <- "PFC_BA46"
    region[grepl("ba46", y, fixed = TRUE)] <- "PFC_BA46"
    region[grepl("hippocampus", y, fixed = TRUE)] <- "hippocampus"
    region[grepl("striatum", y, fixed = TRUE)] <- "striatum"
  }

  if ("brain_region" %in% names(md)) {
    y <- tolower(as.character(md$brain_region))
    region[grepl("ba46", y, fixed = TRUE)] <- "PFC_BA46"
    region[grepl("prefrontal", y, fixed = TRUE)] <- "PFC_BA46"
  }

  if ("characteristics_ch1" %in% names(md)) {
    y <- tolower(as.character(md$characteristics_ch1))
    region[grepl("ba46", y, fixed = TRUE)] <- "PFC_BA46"
    region[grepl("prefrontal", y, fixed = TRUE)] <- "PFC_BA46"
  }

  y <- tolower(as.character(region))
  region[grepl("pre-frontal", y, fixed = TRUE)] <- "PFC_BA46"
  region[grepl("prefrontal", y, fixed = TRUE)] <- "PFC_BA46"
  region[grepl("ba46", y, fixed = TRUE)] <- "PFC_BA46"
  region[grepl("hippocampus", y, fixed = TRUE)] <- "hippocampus"
  region[grepl("striatum", y, fixed = TRUE)] <- "striatum"

  md[, region_standard := region]

  # Numeric / categorical covariates if present
  low_names <- tolower(names(md))

  age_hit <- which(low_names %in% c("age", "age_ch1"))
  if (length(age_hit) > 0) suppressWarnings(md[, age_cov := as.numeric(get(names(md)[age_hit[1]]))])

  rin_hit <- which(low_names %in% c("rin", "rin_ch1"))
  if (length(rin_hit) > 0) suppressWarnings(md[, rin_cov := as.numeric(get(names(md)[rin_hit[1]]))])

  pmi_hit <- which(low_names %in% c("pmi", "pmi_ch1", "postmortem_interval", "post_mortem_interval"))
  if (length(pmi_hit) > 0) suppressWarnings(md[, pmi_cov := as.numeric(get(names(md)[pmi_hit[1]]))])

  sex_hit <- which(low_names %in% c("gender", "sex", "gender_ch1", "sex_ch1"))
  if (length(sex_hit) > 0) md[, sex_cov := as.factor(as.character(get(names(md)[sex_hit[1]])))]

  md
}

collapse_probes_to_symbol <- function(expr, mp) {
  dt <- as.data.table(expr)
  setnames(dt, names(dt)[1], "probe_id")

  dt <- merge(
    mp[, .(probe_id, gene_symbol)],
    dt,
    by = "probe_id",
    allow.cartesian = TRUE
  )

  sample_cols <- setdiff(names(dt), c("probe_id", "gene_symbol"))

  for (cc in sample_cols) {
    suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
  }

  dt[, probe_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = sample_cols]
  dt[, probe_var := apply(.SD, 1, var, na.rm = TRUE), .SDcols = sample_cols]

  setorder(dt, gene_symbol, -probe_var, -probe_mean)
  best <- dt[!duplicated(gene_symbol)]

  out <- best[, c("gene_symbol", sample_cols), with = FALSE]
  setnames(out, "gene_symbol", "feature_id")
  out
}

message("[", timestamp(), "] Step08D4 Affymetrix probe-to-symbol collapse started")
message("[", timestamp(), "] Platform directory: ", PLATFORM_DIR)

mapping_rows <- list()
manifest_rows <- list()
sample_summary_rows <- list()

for (i in seq_len(nrow(DATASETS))) {
  ds <- DATASETS$dataset[i]
  platform <- DATASETS$platform[i]
  expr_path <- DATASETS$expression_probe_path[i]
  sample_path <- DATASETS$sample_path[i]

  if (!file.exists(expr_path)) stop("Missing expression file: ", expr_path)
  if (!file.exists(sample_path)) stop("Missing sample file: ", sample_path)

  message("[", timestamp(), "] Processing ", ds, " platform=", platform)

  mp <- read_platform_table(platform)

  mapping_rows[[length(mapping_rows) + 1]] <- data.table(
    dataset = ds,
    platform = platform,
    n_probe_mapped = uniqueN(mp$probe_id),
    n_gene_symbols = uniqueN(mp$gene_symbol)
  )

  expr <- fread(expr_path, showProgress = FALSE)
  md <- fread(sample_path, showProgress = FALSE)

  md <- standardize_dx_region(ds, md)
  collapsed <- collapse_probes_to_symbol(expr, mp)

  out_dir <- dirname(expr_path)
  out_expr <- file.path(out_dir, paste0(ds, ".expression_symbol.tsv.gz"))
  out_sample <- file.path(out_dir, paste0(ds, ".samples.standardized.tsv"))

  safe_fwrite(collapsed, out_expr)
  safe_fwrite(md, out_sample)

  dx_tab <- md[, .N, by = dx_multi][order(-N)]
  region_tab <- md[, .N, by = region_standard][order(-N)]

  safe_fwrite(dx_tab, file.path(out_dir, paste0(ds, ".dx_summary.tsv")))
  safe_fwrite(region_tab, file.path(out_dir, paste0(ds, ".region_summary.tsv")))

  manifest_rows[[length(manifest_rows) + 1]] <- data.table(
    dataset = ds,
    platform = platform,
    input_expression_probe = expr_path,
    input_samples = sample_path,
    output_expression_symbol = out_expr,
    output_samples_standardized = out_sample,
    n_probe_features = nrow(expr),
    n_symbol_features = nrow(collapsed),
    n_samples = ncol(collapsed) - 1,
    n_dx_groups = uniqueN(md$dx_multi, na.rm = TRUE),
    n_regions = uniqueN(md$region_standard, na.rm = TRUE)
  )

  sample_summary_rows[[length(sample_summary_rows) + 1]] <- data.table(
    dataset = ds,
    dx_summary = paste(dx_tab$dx_multi, dx_tab$N, sep = ":", collapse = ";"),
    region_summary = paste(region_tab$region_standard, region_tab$N, sep = ":", collapse = ";")
  )
}

mapping_summary <- rbindlist(mapping_rows, fill = TRUE)
manifest <- rbindlist(manifest_rows, fill = TRUE)
sample_summary <- rbindlist(sample_summary_rows, fill = TRUE)

safe_fwrite(mapping_summary, file.path(OUT, "75_step08D4_affy_mapping_summary.tsv"))
safe_fwrite(manifest, file.path(OUT, "76_step08D4_symbol_expression_manifest.tsv"))
safe_fwrite(sample_summary, file.path(OUT, "77_step08D4_standardized_sample_summary.tsv"))

summary_lines <- c(
  "# NeuroTRACE Step08D4 Affymetrix probe-to-symbol collapse summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Map Affymetrix probe IDs from GSE53987/GSE12649/GSE21138 to gene symbols and collapse probes to a single gene-symbol expression matrix.",
  "",
  "## Mapping summary",
  paste(capture.output(print(mapping_summary)), collapse = "\n"),
  "",
  "## Manifest",
  paste(capture.output(print(manifest)), collapse = "\n"),
  "",
  "## Sample summary",
  paste(capture.output(print(sample_summary)), collapse = "\n"),
  "",
  "## Interpretation",
  "Use expression_symbol.tsv.gz and samples.standardized.tsv for Step08E cross-disease NTM validation."
)

writeLines(summary_lines, file.path(OUT, "78_step08D4_affy_probe_to_symbol_summary.md"))

message("[", timestamp(), "] Step08D4 done")
message("[", timestamp(), "] Results: ", OUT)
