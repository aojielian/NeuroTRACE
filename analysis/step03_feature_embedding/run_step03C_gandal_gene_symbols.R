#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

BASE <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"
PROJECT <- file.path(BASE, "neurotrace_algorithm_project")
STEP <- file.path(PROJECT, "step03_feature_embedding")
OUT <- file.path(STEP, "results")
FIG <- file.path(STEP, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

DE_IN <- file.path(OUT, "16_step03C_gandal_native_DE_results.tsv.gz")
MODULE_IN <- file.path(OUT, "17_step03C_neurotrace_native_module_weights.tsv")
BRAINSPAN_EMBED <- file.path(OUT, "01_step03A_brainspan_gene_embedding_pca.tsv.gz")
SFARI_FIXED <- file.path(OUT, "13_step03B_sfari_prior_on_embedding_universe_fixed.tsv")

safe_fwrite <- function(x, file) {
  x <- as.data.table(x)
  if (grepl("\\.gz$", file)) {
    tmp <- sub("\\.gz$", "", file)
    data.table::fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
    system(sprintf("gzip -f %s", shQuote(tmp)))
  } else {
    data.table::fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
  }
}

extract_ensembl <- function(x) {
  y <- as.character(x)
  hit <- regmatches(y, regexpr("ENSG[0-9]+", y))
  out <- rep(NA_character_, length(y))
  ok <- nzchar(hit)
  out[ok] <- hit[ok]
  out
}

is_symbol_like <- function(x) {
  y <- as.character(x)
  !grepl("^ENSG[0-9]+", y) & grepl("^[A-Za-z][A-Za-z0-9_.-]*$", y)
}

parse_gtf_attributes <- function(attr, key) {
  # key example: gene_id, gene_name
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regexec(pattern, attr)
  r <- regmatches(attr, m)
  out <- vapply(r, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  out
}

find_gtf_candidates <- function() {
  roots <- c(
    "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD",
    "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA",
    "/gpfs/hpc/home/lijc/lianaoj/SPARK_WGS",
    "/gpfs/hpc/home/lijc/lianaoj/AAAY_data",
    "/gpfs/hpc/home/lijc/lianaoj"
  )

  cmd <- paste(
    "find",
    paste(shQuote(roots[file.exists(roots)]), collapse = " "),
    "\\( -iname '*gencode*.gtf' -o -iname '*gencode*.gtf.gz' -o -iname '*.gtf' -o -iname '*.gtf.gz' \\)",
    "2>/dev/null | head -200"
  )

  res <- suppressWarnings(system(cmd, intern = TRUE))
  res[file.exists(res)]
}

build_map_from_orgdb <- function(ens_ids) {
  message("[", timestamp(), "] Trying AnnotationDbi/org.Hs.eg.db mapping")
  ok <- requireNamespace("AnnotationDbi", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!ok) {
    message("[", timestamp(), "] org.Hs.eg.db not available")
    return(data.table())
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
    columns = c("ENSEMBL", "SYMBOL", "GENENAME")
  )
  mp <- as.data.table(mp)
  setnames(mp, c("ENSEMBL", "SYMBOL", "GENENAME"), c("ensembl_id_clean", "gene_symbol_fixed", "gene_name"), skip_absent = TRUE)
  mp <- mp[!is.na(ensembl_id_clean) & !is.na(gene_symbol_fixed)]
  mp <- unique(mp, by = c("ensembl_id_clean", "gene_symbol_fixed"))
  mp[, mapping_source := "org.Hs.eg.db"]
  mp
}

build_map_from_gtf <- function(ens_ids) {
  message("[", timestamp(), "] Trying local GTF mapping")
  gtfs <- find_gtf_candidates()

  if (length(gtfs) == 0) {
    message("[", timestamp(), "] No local GTF candidates found")
    return(data.table())
  }

  message("[", timestamp(), "] GTF candidates:")
  message(paste(gtfs, collapse = "\n"))

  ens_set <- unique(na.omit(ens_ids))
  rows <- list()

  for (gtf in gtfs) {
    message("[", timestamp(), "] Reading GTF candidate: ", gtf)

    cmd <- if (grepl("\\.gz$", gtf)) {
      paste("zcat", shQuote(gtf))
    } else {
      paste("cat", shQuote(gtf))
    }

    # Use fread cmd; keep only gene rows where possible
    dt <- tryCatch(
      fread(cmd = paste(cmd, "| awk '$3==\"gene\" {print}'"), sep = "\t", header = FALSE, quote = ""),
      error = function(e) data.table()
    )

    if (nrow(dt) == 0) next
    if (ncol(dt) < 9) next

    attr <- dt[[9]]
    gene_id <- parse_gtf_attributes(attr, "gene_id")
    gene_name <- parse_gtf_attributes(attr, "gene_name")

    gene_id_clean <- sub("\\..*$", "", gene_id)

    mp <- data.table(
      ensembl_id_clean = gene_id_clean,
      gene_symbol_fixed = gene_name,
      gene_name = gene_name,
      mapping_source = paste0("GTF:", gtf)
    )

    mp <- mp[
      !is.na(ensembl_id_clean) &
      ensembl_id_clean %in% ens_set &
      !is.na(gene_symbol_fixed) &
      gene_symbol_fixed != ""
    ]

    if (nrow(mp) > 0) {
      rows[[length(rows) + 1]] <- unique(mp, by = c("ensembl_id_clean", "gene_symbol_fixed"))
    }

    # 如果已经覆盖大部分，就不用继续扫
    if (length(rows) > 0) {
      tmp <- rbindlist(rows, fill = TRUE)
      cov <- uniqueN(tmp$ensembl_id_clean) / length(ens_set)
      message("[", timestamp(), "] Current GTF mapping coverage: ", round(cov, 4))
      if (cov > 0.80) break
    }
  }

  if (length(rows) == 0) return(data.table())
  rbindlist(rows, fill = TRUE)
}

collapse_duplicate_symbols <- function(de) {
  de <- as.data.table(de)

  # 对同一 gene_symbol_fixed 多个 transcript/duplicate probe，保留 p 最小、效应绝对值较大的一条
  # 注意：data.table::setorder 不能直接使用 -abs(beta_ASD_vs_Control)，所以先显式生成排序列
  de[, abs_beta_for_collapse := abs(beta_ASD_vs_Control)]
  de[, p_value_for_collapse := as.numeric(p_value)]

  setorder(de, gene_symbol_fixed, p_value_for_collapse, -abs_beta_for_collapse)

  de2 <- de[!duplicated(gene_symbol_fixed)]
  de2[, c("abs_beta_for_collapse", "p_value_for_collapse") := NULL]
  de2
}

derive_symbol_modules <- function(de_symbol, top_ns = c(50, 100, 200, 500)) {
  de <- copy(de_symbol)
  de[, signed_score := sign(beta_ASD_vs_Control) * (-log10(pmax(p_value, 1e-300)))]
  de[, abs_score := abs(signed_score)]

  rows <- list()

  for (N in top_ns) {
    up <- de[beta_ASD_vs_Control > 0][order(-signed_score)][seq_len(min(N, .N))]
    if (nrow(up) > 0) {
      up[, `:=`(
        program = "NTM1_ASD_up",
        top_n = N,
        direction = "ASD_up",
        weight = signed_score
      )]
      rows[[length(rows) + 1]] <- up
    }

    down <- de[beta_ASD_vs_Control < 0][order(signed_score)][seq_len(min(N, .N))]
    if (nrow(down) > 0) {
      down[, `:=`(
        program = "NTM2_ASD_down",
        top_n = N,
        direction = "ASD_down",
        weight = signed_score
      )]
      rows[[length(rows) + 1]] <- down
    }

    signed <- de[order(-abs_score)][seq_len(min(N, .N))]
    if (nrow(signed) > 0) {
      signed[, `:=`(
        program = "NTM3_ASD_signed",
        top_n = N,
        direction = ifelse(beta_ASD_vs_Control >= 0, "ASD_up", "ASD_down"),
        weight = signed_score
      )]
      rows[[length(rows) + 1]] <- signed
    }
  }

  mod <- rbindlist(rows, fill = TRUE)
  mod[, weight_z := as.numeric(scale(weight)), by = .(program, top_n)]
  mod[, rank_within_program := frank(-abs(weight), ties.method = "first"), by = .(program, top_n)]

  keep_cols <- c(
    "program", "top_n", "direction", "rank_within_program",
    "gene_symbol_fixed", "ensembl_id_clean", "gene_symbol",
    "weight", "weight_z",
    "beta_ASD_vs_Control", "se", "t", "p_value", "fdr",
    "mean_expression", "sd_expression", "mapping_source"
  )

  mod[, ..keep_cols]
}

message("[", timestamp(), "] Step03C2 fixing Gandal gene symbols started")

for (f in c(DE_IN, MODULE_IN, BRAINSPAN_EMBED, SFARI_FIXED)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

de <- fread(DE_IN)
modules_old <- fread(MODULE_IN)
brain <- fread(BRAINSPAN_EMBED, nThread = 1)
sfari <- fread(SFARI_FIXED)

de[, original_gene_id := gene_symbol]
de[, ensembl_id_clean := extract_ensembl(original_gene_id)]
de[, original_is_symbol_like := is_symbol_like(original_gene_id)]

ens_ids <- unique(na.omit(de$ensembl_id_clean))

message("[", timestamp(), "] DE rows: ", nrow(de))
message("[", timestamp(), "] Unique clean Ensembl IDs: ", length(ens_ids))

map1 <- build_map_from_orgdb(ens_ids)
map2 <- data.table()

if (nrow(map1) == 0 || uniqueN(map1$ensembl_id_clean) / length(ens_ids) < 0.50) {
  map2 <- build_map_from_gtf(ens_ids)
}

mapping <- rbindlist(list(map1, map2), fill = TRUE)
if (nrow(mapping) > 0) {
  # Prefer orgdb, then first GTF
  mapping[, source_rank := fifelse(mapping_source == "org.Hs.eg.db", 1L, 2L)]
  setorder(mapping, ensembl_id_clean, source_rank)
  mapping <- mapping[!duplicated(ensembl_id_clean)]
} else {
  mapping <- data.table(
    ensembl_id_clean = character(),
    gene_symbol_fixed = character(),
    gene_name = character(),
    mapping_source = character()
  )
}

safe_fwrite(mapping, file.path(OUT, "21_step03C2_ensembl_to_symbol_mapping.tsv"))

de2 <- merge(de, mapping, by = "ensembl_id_clean", all.x = TRUE)

# 如果本来就是 symbol-like，则保留原始 symbol；否则使用 mapping
de2[, gene_symbol_final := fifelse(
  !is.na(gene_symbol_fixed) & gene_symbol_fixed != "",
  gene_symbol_fixed,
  fifelse(original_is_symbol_like, original_gene_id, NA_character_)
)]

de2[, mapping_status := fifelse(
  !is.na(gene_symbol_fixed) & gene_symbol_fixed != "",
  "mapped_ensembl_to_symbol",
  fifelse(original_is_symbol_like, "original_symbol_like", "unmapped")
)]

safe_fwrite(de2, file.path(OUT, "22_step03C2_native_DE_with_symbol_mapping.tsv.gz"))

# 只保留能映射到 symbol 的基因
de_mapped <- de2[!is.na(gene_symbol_final) & gene_symbol_final != ""]
setnames(de_mapped, "gene_symbol_final", "gene_symbol_fixed", skip_absent = TRUE)

de_symbol <- collapse_duplicate_symbols(de_mapped)
safe_fwrite(de_symbol, file.path(OUT, "23_step03C2_native_DE_symbol_level.tsv.gz"))

symbol_modules <- derive_symbol_modules(de_symbol, top_ns = c(50, 100, 200, 500))
safe_fwrite(symbol_modules, file.path(OUT, "24_step03C2_neurotrace_native_module_weights_symbol.tsv"))

module_summary <- symbol_modules[, .(
  n_genes = .N,
  n_ASD_up = sum(direction == "ASD_up"),
  n_ASD_down = sum(direction == "ASD_down"),
  median_abs_weight = median(abs(weight), na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  min_fdr = min(fdr, na.rm = TRUE)
), by = .(program, top_n)]

safe_fwrite(module_summary, file.path(OUT, "25_step03C2_native_symbol_module_summary.tsv"))

# BrainSpan overlap
brain[, gene_key_upper := toupper(gene_key)]
symbol_modules[, gene_key_upper := toupper(gene_symbol_fixed)]

brain_overlap <- symbol_modules[, .(
  n_module_genes = .N,
  n_unique_module_genes = uniqueN(gene_symbol_fixed),
  n_overlap_brainspan = uniqueN(gene_symbol_fixed[gene_key_upper %in% brain$gene_key_upper]),
  overlap_rate_brainspan = uniqueN(gene_symbol_fixed[gene_key_upper %in% brain$gene_key_upper]) / uniqueN(gene_symbol_fixed)
), by = .(program, top_n)]

safe_fwrite(brain_overlap, file.path(OUT, "26_step03C2_native_module_brainspan_overlap.tsv"))

# SFARI overlap
sfari_cols <- names(sfari)
sfari_gene_col <- if ("gene_key" %in% sfari_cols) "gene_key" else NA_character_
if (!is.na(sfari_gene_col)) {
  sfari[, gene_key_upper := toupper(gene_key)]
  sfari_sets <- c(
    "is_SFARI_all_fixed",
    "is_SFARI_score_le1_fixed",
    "is_SFARI_score_le2_fixed",
    "is_SFARI_syndromic_fixed",
    "is_SFARI_S_plus_1_fixed",
    "is_SFARI_S_plus_1_2_fixed"
  )
  sfari_sets <- intersect(sfari_sets, names(sfari))

  sfari_long <- rbindlist(lapply(sfari_sets, function(ss) {
    data.table(
      risk_set = ss,
      gene_key_upper = sfari[get(ss) == 1, gene_key_upper]
    )
  }))

  sfari_overlap <- symbol_modules[, {
    genes <- unique(gene_key_upper)
    rbindlist(lapply(unique(sfari_long$risk_set), function(rs) {
      risk_genes <- unique(sfari_long[risk_set == rs, gene_key_upper])
      data.table(
        risk_set = rs,
        n_module_genes = length(genes),
        overlap_n = sum(genes %in% risk_genes),
        overlap_rate = sum(genes %in% risk_genes) / length(genes)
      )
    }))
  }, by = .(program, top_n)]

  safe_fwrite(sfari_overlap, file.path(OUT, "27_step03C2_native_module_sfari_overlap.tsv"))
}

# audit
audit <- data.table(
  item = c(
    "DE_input_rows",
    "unique_clean_ensembl_ids",
    "mapped_ensembl_ids",
    "mapping_coverage_unique_ensembl",
    "DE_rows_mapped_or_symbol",
    "symbol_level_DE_rows_after_duplicate_collapse",
    "symbol_module_rows",
    "BrainSpan_embedding_genes"
  ),
  value = c(
    nrow(de),
    length(ens_ids),
    uniqueN(mapping$ensembl_id_clean),
    round(uniqueN(mapping$ensembl_id_clean) / max(length(ens_ids), 1), 4),
    nrow(de_mapped),
    nrow(de_symbol),
    nrow(symbol_modules),
    nrow(brain)
  )
)
safe_fwrite(audit, file.path(OUT, "28_step03C2_mapping_audit.tsv"))

top_symbol_de <- rbind(
  de_symbol[beta_ASD_vs_Control > 0][order(p_value)][1:min(.N, 100)][, direction := "ASD_up"],
  de_symbol[beta_ASD_vs_Control < 0][order(p_value)][1:min(.N, 100)][, direction := "ASD_down"]
)
safe_fwrite(top_symbol_de, file.path(OUT, "29_step03C2_top_native_symbol_DE_genes.tsv"))

summary_lines <- c(
  "# NeuroTRACE Step03C2 Gandal gene-symbol correction summary",
  "",
  paste0("Generated: ", timestamp()),
  "",
  "## Purpose",
  "Step03C produced NeuroTRACE-native modules from Gandal 2022 without using DevMap/DevBridge seeds, but the gene identifiers were Ensembl-style IDs. Step03C2 maps these IDs to gene symbols and rebuilds symbol-level NTM modules for BrainSpan/SFARI/graph integration.",
  "",
  "## Mapping audit",
  paste0("- DE input rows: ", nrow(de)),
  paste0("- Unique clean Ensembl IDs: ", length(ens_ids)),
  paste0("- Mapped Ensembl IDs: ", uniqueN(mapping$ensembl_id_clean)),
  paste0("- Mapping coverage: ", round(uniqueN(mapping$ensembl_id_clean) / max(length(ens_ids), 1), 4)),
  paste0("- Symbol-level DE rows after duplicate collapse: ", nrow(de_symbol)),
  paste0("- Symbol module rows: ", nrow(symbol_modules)),
  "",
  "## Main outputs",
  "- `21_step03C2_ensembl_to_symbol_mapping.tsv`",
  "- `23_step03C2_native_DE_symbol_level.tsv.gz`",
  "- `24_step03C2_neurotrace_native_module_weights_symbol.tsv`",
  "- `25_step03C2_native_symbol_module_summary.tsv`",
  "- `26_step03C2_native_module_brainspan_overlap.tsv`",
  "- `27_step03C2_native_module_sfari_overlap.tsv`",
  "- `29_step03C2_top_native_symbol_DE_genes.tsv`",
  "",
  "## Interpretation",
  "Use Step03C2 symbol-level NTM modules for all downstream NeuroTRACE real-data analyses. Do not use the Ensembl-style Step03C module file directly for BrainSpan alignment."
)
writeLines(summary_lines, file.path(OUT, "30_step03C2_gene_symbol_fix_summary.md"))

message("[", timestamp(), "] Step03C2 done")
message("[", timestamp(), "] Mapping coverage: ", round(uniqueN(mapping$ensembl_id_clean) / max(length(ens_ids), 1), 4))
message("[", timestamp(), "] Symbol-level DE rows: ", nrow(de_symbol))
message("[", timestamp(), "] Results: ", OUT)
