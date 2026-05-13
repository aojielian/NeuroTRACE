suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Supplementary_FigureS3")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ZHOU_COMBINED_FILE <- file.path(ROOT, "step10_zhou2022_genesets_fixed/results/12_step10G4_combined_Zhou_enrichment.tsv")
ZHOU_TOP_FILE <- file.path(ROOT, "step10_zhou2022_genesets_fixed/results/13_step10G4_top_Zhou_enrichment_hits.tsv")
SATTER_COMBINED_FILE <- file.path(ROOT, "step10_independent_genetic_risk/results/04_step10G2_combined_Satterstrom102_enrichment.tsv")
DEGREE_NULL_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12B_degree_matched_graph_null/07b_step12B_decision_table.corrected.tsv")

PDF_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S3_v5.pdf")
PNG_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S3_v5.png")
DIAG_OUT <- file.path(OUTDIR, "Supplementary_Figure_S3_v5_diagnostics.txt")
PLOTDIR  <- file.path(OUTDIR, "plotdata_v5")
dir.create(PLOTDIR, recursive = TRUE, showWarnings = FALSE)

COL_TEXT   <- "#263238"
COL_MUTED  <- "#6B7280"
COL_LIGHT  <- "#D7DEE6"
COL_GREEN  <- "#4F9A73"
COL_BLUE   <- "#4777A8"
COL_ORANGE <- "#D96C4A"
COL_PURPLE <- "#7A5AA6"
COL_RED    <- "#B64A4A"
COL_GREY   <- "#8A96A3"

font_family <- "Helvetica"

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.30, colour = COL_TEXT),
      axis.text = element_text(colour = COL_TEXT, size = base_size),
      axis.title = element_text(colour = COL_TEXT, size = base_size + 0.8),
      strip.background = element_rect(fill = "white", colour = COL_LIGHT, linewidth = 0.35),
      strip.text = element_text(face = "bold", colour = COL_TEXT, size = base_size),
      legend.title = element_text(colour = COL_TEXT, size = base_size - 0.2),
      legend.text = element_text(colour = COL_TEXT, size = base_size - 0.7),
      legend.key.size = unit(0.34, "cm"),
      plot.tag = element_text(face = "bold", size = 16, colour = "black"),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(6, 7, 6, 7)
    )
}

read_dt <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
  fread(path)
}

pick_col <- function(dt, patterns, numeric_only = FALSE, exclude = character()) {
  nm <- names(dt)
  nm_low <- tolower(nm)
  ok <- rep(TRUE, length(nm))
  if (numeric_only) ok <- vapply(dt, is.numeric, logical(1))
  if (length(exclude) > 0) ok <- ok & !(nm_low %in% tolower(exclude))
  for (pat in patterns) {
    hit <- which(grepl(pat, nm_low, perl = TRUE) & ok)
    if (length(hit) > 0) return(nm[hit[1]])
  }
  NA_character_
}

pretty_module <- function(x) {
  y <- as.character(x)
  y <- gsub("^module:", "", y)
  y <- gsub("\\|.*$", "", y)
  y <- gsub("_top[0-9]+$", "", y)
  out <- y
  out[grepl("NTM1", y, ignore.case = TRUE)] <- "NTM1\nASD-up"
  out[grepl("NTM2", y, ignore.case = TRUE)] <- "NTM2\nASD-down"
  out[grepl("NTM3", y, ignore.case = TRUE)] <- "NTM3\nASD-signed"
  out
}

module_levels <- c("NTM1\nASD-up", "NTM2\nASD-down", "NTM3\nASD-signed")

clean_resource <- function(x) {
  y <- as.character(x)
  y <- gsub("^risk_set:", "", y)
  y <- gsub("^gene_set:", "", y)
  y <- gsub("_fixed$", "", y)
  y <- gsub("_", " ", y)
  y <- gsub("TARGET GENES", "", y, ignore.case = TRUE)
  y <- gsub("GENES", "", y, ignore.case = TRUE)
  y <- gsub("\\s+", " ", y)
  trimws(y)
}

short_resource <- function(x, max_chars = 36) {
  y <- clean_resource(x)
  y <- gsub("Known ASD NDD", "Known ASD/NDD", y, ignore.case = TRUE)
  y <- gsub("Archetype", "Arch.", y, ignore.case = TRUE)
  y <- gsub("neurotransmission", "neurotrans.", y, ignore.case = TRUE)
  y <- gsub("vesicle mediated transport", "vesicle transport", y, ignore.case = TRUE)
  y <- gsub("constrained background", "constrained bg", y, ignore.case = TRUE)
  ifelse(nchar(y) > max_chars, paste0(substr(y, 1, max_chars - 3), "..."), y)
}

resource_category <- function(x) {
  y <- clean_resource(x)
  out <- rep("Other", length(y))
  out[grepl("Archetype", y, ignore.case = TRUE)] <- "Zhou archetypes"
  out[grepl("Known|ASD/NDD|Meta", y, ignore.case = TRUE)] <- "Known ASD/NDD"
  out[grepl("TDT|overtransmitted|inherited", y, ignore.case = TRUE)] <- "TDT / inherited"
  out[grepl("constrain|LOEUF|pLI", y, ignore.case = TRUE)] <- "Constraint"
  out[grepl("deNovo|de novo|dnLoF|LoF", y, ignore.case = TRUE)] <- "de novo / LoF"
  out
}

standardize_enrichment <- function(dt) {
  module_col <- pick_col(dt, c("^module$", "module_family", "^program$", "^ntm$", "query", "signature"))
  topn_col <- pick_col(dt, c("^top_n$", "^topn$", "module_top_n", "cutoff", "threshold"))
  resource_col <- pick_col(dt, c("gene_set", "geneset", "set_name", "term", "resource", "category", "risk_set", "class", "label", "pathway"))
  source_col <- pick_col(dt, c("source_type", "gene_set_source", "^source$", "collection", "database"))
  or_col <- pick_col(dt, c("median_or", "odds_ratio", "^or$", "enrichment", "enrichment_ratio", "effect", "estimate"), numeric_only = TRUE)
  p_col <- pick_col(dt, c("min_p", "^p$", "p_value", "pval", "p\\.value"), numeric_only = TRUE)
  fdr_col <- pick_col(dt, c("min_fdr", "fdr", "padj", "adj.*p", "q_value", "qval", "global_fdr"), numeric_only = TRUE)
  overlap_col <- pick_col(dt, c("n_overlap", "overlap_count", "overlap", "^k$", "^a$", "n_hits", "n_genes"), numeric_only = TRUE)
  total_col <- pick_col(dt, c("n_total", "n_tested", "n_gene_sets", "n_risk_sets", "n_constraint_sets", "set_size", "n_set"), numeric_only = TRUE)
  status_col <- pick_col(dt, c("status", "decision", "call"))

  out <- data.table(
    module = if (!is.na(module_col)) pretty_module(dt[[module_col]]) else NA_character_,
    top_n = if (!is.na(topn_col)) suppressWarnings(as.numeric(dt[[topn_col]])) else NA_real_,
    resource = if (!is.na(resource_col)) as.character(dt[[resource_col]]) else NA_character_,
    resource_source = if (!is.na(source_col)) as.character(dt[[source_col]]) else NA_character_,
    odds_ratio = if (!is.na(or_col)) as.numeric(dt[[or_col]]) else NA_real_,
    p_value = if (!is.na(p_col)) as.numeric(dt[[p_col]]) else NA_real_,
    fdr = if (!is.na(fdr_col)) as.numeric(dt[[fdr_col]]) else NA_real_,
    overlap_n = if (!is.na(overlap_col)) as.numeric(dt[[overlap_col]]) else NA_real_,
    total_n = if (!is.na(total_col)) as.numeric(dt[[total_col]]) else NA_real_,
    status = if (!is.na(status_col)) as.character(dt[[status_col]]) else NA_character_
  )

  out[, resource := fifelse(is.na(resource) | resource == "", resource_source, resource)]
  out[, resource := clean_resource(resource)]
  out <- out[!is.na(resource) & resource != ""]
  out
}

## ============================================================
## Read and standardize inputs
## ============================================================

zhou_combined_raw <- read_dt(ZHOU_COMBINED_FILE)
zhou_top_raw <- read_dt(ZHOU_TOP_FILE)
satter_raw <- read_dt(SATTER_COMBINED_FILE)
degree_dt <- read_dt(DEGREE_NULL_FILE)

zhou_dt <- standardize_enrichment(zhou_combined_raw)
zhou_top_dt <- standardize_enrichment(zhou_top_raw)
satter_dt <- standardize_enrichment(satter_raw)

## keep top500 when available
zhou_main <- zhou_dt[!is.na(module) & module %in% module_levels]
if (500 %in% zhou_main$top_n) zhou_main <- zhou_main[top_n == 500]

zhou_main[, score := -log10(pmax(fdr, 1e-300))]
zhou_main[!is.finite(score), score := NA_real_]
zhou_main[, category := resource_category(resource)]

## ============================================================
## Panel A: Zhou category-level convergence heatmap
## ============================================================

a_dt <- zhou_main[
  ,
  .(
    score = max(score, na.rm = TRUE),
    max_or = suppressWarnings(max(odds_ratio, na.rm = TRUE)),
    max_overlap = suppressWarnings(max(overlap_n, na.rm = TRUE))
  ),
  by = .(module, category)
]

a_dt[!is.finite(score), score := NA_real_]
a_dt[!is.finite(max_or), max_or := NA_real_]
a_dt[!is.finite(max_overlap), max_overlap := NA_real_]

cat_levels <- c("Zhou archetypes", "Constraint", "Known ASD/NDD", "TDT / inherited", "de novo / LoF", "Other")
cat_levels <- cat_levels[cat_levels %in% unique(a_dt$category)]

a_dt[, module := factor(module, levels = module_levels)]
a_dt[, category := factor(category, levels = rev(cat_levels))]

fwrite(a_dt, file.path(PLOTDIR, "S3A_zhou_category_convergence.tsv"), sep = "\t")

panelA <- ggplot(a_dt, aes(x = module, y = category, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.40, width = 0.90, height = 0.82) +
  scale_fill_gradient(low = "#EEF2F6", high = COL_GREEN, name = expression(-log[10](FDR))) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(size = 7.2),
    axis.text.y = element_text(size = 7.2),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

## ============================================================
## Panel B: Zhou 2022 top gene-set convergence
## ============================================================

b_dt <- copy(zhou_top_dt)
b_dt <- b_dt[!is.na(module) & module %in% module_levels]

if (500 %in% b_dt$top_n) {
  b_dt <- b_dt[top_n == 500]
}

b_dt[, score := -log10(pmax(fdr, 1e-300))]
b_dt[!is.finite(score), score := NA_real_]

b_best <- b_dt[
  ,
  .(
    score = max(score, na.rm = TRUE),
    odds_ratio = suppressWarnings(max(odds_ratio, na.rm = TRUE)),
    overlap_n = suppressWarnings(max(overlap_n, na.rm = TRUE))
  ),
  by = .(module, resource)
]

b_best[!is.finite(score), score := NA_real_]
b_best[!is.finite(odds_ratio), odds_ratio := NA_real_]
b_best[!is.finite(overlap_n), overlap_n := NA_real_]

top_resources <- b_best[
  ,
  .(max_score = max(score, na.rm = TRUE)),
  by = resource
][order(-max_score)]$resource
top_resources <- top_resources[seq_len(min(12, length(top_resources)))]

b_plot <- b_best[resource %in% top_resources]
b_plot[, module := factor(module, levels = module_levels)]
b_plot[, resource_label := short_resource(resource, 36)]
b_plot[, resource_label := factor(resource_label, levels = rev(unique(b_plot[order(score)]$resource_label)))]

fwrite(b_plot, file.path(PLOTDIR, "S3B_zhou2022_top_convergence.tsv"), sep = "\t")

panelB <- ggplot(b_plot, aes(x = module, y = resource_label)) +
  geom_point(aes(size = pmax(overlap_n, 1), colour = score), alpha = 0.92) +
  scale_colour_gradient(low = COL_BLUE, high = COL_RED, name = expression(-log[10](FDR))) +
  scale_size_continuous(range = c(1.8, 6.0), name = "Overlap") +
  labs(x = NULL, y = NULL, tag = "B") +
  theme_pub(7.2) +
  theme(
    axis.text.x = element_text(size = 7.0),
    axis.text.y = element_text(size = 6.7),
    legend.position = "right"
  )

## ============================================================
## Panel C: Satterstrom ASD exome-risk enrichment
## ============================================================

c_dt <- copy(satter_dt)
c_dt <- c_dt[!is.na(module) & module %in% module_levels]
c_dt <- c_dt[top_n %in% c(200, 500)]

c_dt[, score := -log10(pmax(fdr, 1e-300))]
c_dt[!is.finite(score), score := NA_real_]

c_best <- c_dt[
  ,
  .(
    odds_ratio = suppressWarnings(max(odds_ratio, na.rm = TRUE)),
    score = suppressWarnings(max(score, na.rm = TRUE)),
    overlap_n = suppressWarnings(max(overlap_n, na.rm = TRUE)),
    total_n = suppressWarnings(max(total_n, na.rm = TRUE))
  ),
  by = .(module, top_n)
]

c_best[!is.finite(odds_ratio), odds_ratio := NA_real_]
c_best[!is.finite(score), score := NA_real_]
c_best[!is.finite(overlap_n), overlap_n := NA_real_]
c_best[!is.finite(total_n), total_n := NA_real_]

c_best[, module := factor(module, levels = module_levels)]
c_best[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]
c_best[, overlap_label := fifelse(
  !is.na(overlap_n) & !is.na(total_n),
  paste0(overlap_n, "/", total_n),
  fifelse(!is.na(overlap_n), as.character(overlap_n), "")
)]

fwrite(c_best, file.path(PLOTDIR, "S3C_satterstrom102_enrichment.tsv"), sep = "\t")

panelC <- ggplot(c_best, aes(x = cutoff, y = odds_ratio, fill = module)) +
  geom_hline(yintercept = 1, colour = COL_LIGHT, linewidth = 0.35, linetype = "dashed") +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "white", linewidth = 0.20) +
  geom_text(
    aes(label = overlap_label),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 2.35,
    colour = COL_TEXT,
    family = font_family
  ) +
  scale_fill_manual(values = c(
    "NTM1\nASD-up" = COL_ORANGE,
    "NTM2\nASD-down" = COL_BLUE,
    "NTM3\nASD-signed" = COL_PURPLE
  ), drop = FALSE, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Satterstrom enrichment OR", tag = "C") +
  theme_pub(7.8) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 7.2)
  )

## ============================================================
## Panel D: degree-matched graph null summary
## ============================================================

d_dt <- copy(degree_dt)

d_dt[, module := pretty_module(module_family)]
d_dt[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]

d_dt[, analysis_label := fifelse(
  grepl("SFARI", analysis, ignore.case = TRUE),
  "NTM2 SFARI degree null",
  "Transport degree null"
)]

d_dt[, neglog_fdr := -log10(pmax(primary_fdr, 1e-300))]
d_dt[!is.finite(neglog_fdr), neglog_fdr := NA_real_]

d_dt[, analysis_label := factor(
  analysis_label,
  levels = c("Transport degree null", "NTM2 SFARI degree null")
)]

fwrite(d_dt, file.path(PLOTDIR, "S3D_degree_matched_null_summary.tsv"), sep = "\t")

shape_vals <- c(
  "Transport degree null" = 21,
  "NTM2 SFARI degree null" = 24
)

panelD <- ggplot(d_dt, aes(x = cutoff, y = analysis_label)) +
  geom_point(
    aes(fill = neglog_fdr, shape = analysis_label),
    size = 4.0,
    colour = "white",
    stroke = 0.35,
    alpha = 0.95
  ) +
  scale_fill_gradient(
    low = COL_BLUE,
    high = COL_RED,
    name = expression(-log[10](FDR))
  ) +
  scale_shape_manual(values = shape_vals, name = NULL) +
  labs(x = NULL, y = NULL, tag = "D") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(size = 7.2),
    axis.text.y = element_text(size = 7.2),
    legend.position = "right"
  )


## ============================================================
## Diagnostics
## ============================================================

diag <- c(
  "Supplementary Figure S3 v5 diagnostics",
  paste("Generated:", as.character(Sys.time())),
  "",
  paste("ZHOU_COMBINED_FILE:", ZHOU_COMBINED_FILE),
  paste("ZHOU_TOP_FILE:", ZHOU_TOP_FILE),
  paste("SATTER_COMBINED_FILE:", SATTER_COMBINED_FILE),
  paste("DEGREE_NULL_FILE:", DEGREE_NULL_FILE),
  "",
  paste("zhou_combined rows:", nrow(zhou_combined_raw)),
  paste("zhou_top rows:", nrow(zhou_top_raw)),
  paste("satter rows:", nrow(satter_raw)),
  paste("degree rows:", nrow(degree_dt)),
  "",
  "zhou_combined columns:",
  paste(names(zhou_combined_raw), collapse = ", "),
  "",
  "zhou_top columns:",
  paste(names(zhou_top_raw), collapse = ", "),
  "",
  "satter columns:",
  paste(names(satter_raw), collapse = ", "),
  "",
  "degree columns:",
  paste(names(degree_dt), collapse = ", "),
  "",
  "plotdata files:",
  paste(list.files(PLOTDIR, full.names = FALSE), collapse = ", ")
)

writeLines(diag, DIAG_OUT)

## ============================================================
## Export
## ============================================================

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 2,
    widths = unit(c(1.02, 1.12), "null"),
    heights = unit(c(1.00, 1.00), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 10.8, height = 8.2, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 3240, height = 2460, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Supplementary Figure S3 v5 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plotdata:\n")
cat(PLOTDIR, "\n")
cat("Diagnostics:\n")
cat(DIAG_OUT, "\n")
