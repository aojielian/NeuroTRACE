suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Supplementary_FigureS5")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

STEP12A_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12A_v4_full/04_step12A_decision_table.tsv")
STEP12B_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12B_degree_matched_graph_null/07b_step12B_decision_table.corrected.tsv")
STEP12C_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12C_stage_alignment_weight_sensitivity_v3/05_step12C_v3_decision_table.tsv")
STEP12D_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12D_gnomad_constraint/05_step12D1_decision_table.tsv")

DISEASE_SUMMARY_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/02_step11E_v2_disease_level_summary.tsv")
MODULE_SUMMARY_FILE  <- file.path(ROOT, "step11_robustness_sensitivity/results/03_step11E_v2_module_level_summary.tsv")

PDF_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S5_v3.pdf")
PNG_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S5_v3.png")
DIAG_OUT <- file.path(OUTDIR, "Supplementary_Figure_S5_v3_diagnostics.txt")
PLOTDIR  <- file.path(OUTDIR, "plotdata_v3")
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
    warning("Missing input file: ", path)
    return(NULL)
  }
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("zcat", shQuote(path)))
  } else {
    fread(path)
  }
}

pick_col <- function(dt, patterns, numeric_only = FALSE, exclude = character()) {
  if (is.null(dt)) return(NA_character_)
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
  y <- gsub("_top[0-9]+$", "", y)
  out <- y
  out[grepl("NTM1", y, ignore.case = TRUE)] <- "NTM1\nASD-up"
  out[grepl("NTM2", y, ignore.case = TRUE)] <- "NTM2\nASD-down"
  out[grepl("NTM3", y, ignore.case = TRUE)] <- "NTM3\nASD-signed"
  out
}

module_levels <- c(
  "NTM1\nASD-up",
  "NTM2\nASD-down",
  "NTM3\nASD-signed"
)

module_cols <- c(
  "NTM1\nASD-up" = COL_ORANGE,
  "NTM2\nASD-down" = COL_BLUE,
  "NTM3\nASD-signed" = COL_PURPLE
)

make_placeholder <- function(tag, label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.55,
      label = label,
      size = 3.0,
      colour = COL_MUTED,
      family = font_family,
      lineheight = 0.95
    ) +
    xlim(0, 1) + ylim(0, 1) +
    labs(tag = tag) +
    theme_void(base_family = font_family) +
    theme(
      plot.tag = element_text(face = "bold", size = 16, colour = "black"),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(6, 7, 6, 7)
    )
}

## ============================================================
## Read inputs
## ============================================================

step12a <- read_dt(STEP12A_FILE)
step12b <- read_dt(STEP12B_FILE)
step12c <- read_dt(STEP12C_FILE)
step12d <- read_dt(STEP12D_FILE)
disease_summary <- read_dt(DISEASE_SUMMARY_FILE)
module_summary <- read_dt(MODULE_SUMMARY_FILE)

## ============================================================
## Panel A: bootstrap module stability
## ============================================================

panelA <- NULL

if (!is.null(step12a) && nrow(step12a) > 0) {
  a_dt <- copy(step12a)

  a_dt[, module := pretty_module(module)]
  a_dt <- a_dt[module %in% module_levels]

  ## Keep interpretable median-level bootstrap stability metrics.
  ## Jaccard is omitted from the figure because it is sensitive to boundary
  ## changes in ranked gene sets and is better left in Supplementary Table S10.
  metric_map <- c(
    "frozen_recall_median" = "Frozen recall",
    "bootstrap_precision_median" = "Bootstrap precision",
    "sign_concordance_overlap_median" = "Sign concordance",
    "full_vs_boot_signed_score_spearman_median" = "Signed-score Spearman"
  )

  keep_metrics <- names(metric_map)[names(metric_map) %in% names(a_dt)]

  if (length(keep_metrics) == 0) {
    stop("No expected bootstrap stability metric columns were found in STEP12A file.")
  }

  a_long <- melt(
    a_dt,
    id.vars = "module",
    measure.vars = keep_metrics,
    variable.name = "metric",
    value.name = "value"
  )

  a_long[, metric_label := unname(metric_map[as.character(metric)])]
  a_long[, module := factor(module, levels = module_levels)]
  a_long[, metric_label := factor(
    metric_label,
    levels = rev(c(
      "Frozen recall",
      "Bootstrap precision",
      "Sign concordance",
      "Signed-score Spearman"
    ))
  )]

  fwrite(a_long, file.path(PLOTDIR, "S5A_bootstrap_module_stability.tsv"), sep = "\t")

  panelA <- ggplot(a_long, aes(x = module, y = metric_label, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.42, width = 0.88, height = 0.78) +
    scale_fill_gradient(
      low = "#EEF2F6",
      high = COL_GREEN,
      limits = c(0, 1),
      name = "Stability"
    ) +
    labs(x = NULL, y = NULL, tag = "A") +
    theme_pub(7.8) +
    theme(
      axis.text.x = element_text(size = 7.1),
      axis.text.y = element_text(size = 7.2),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

if (is.null(panelA)) {
  panelA <- make_placeholder("A", "Bootstrap module stability\ncould not be generated.\nSee diagnostics.")
}


## ============================================================
## Panel B: stage-alignment stability across top200/top500
## ============================================================

panelB <- NULL

if (!is.null(step12c) && nrow(step12c) > 0) {
  b_dt <- copy(step12c)

  b_dt[, module := pretty_module(module_family)]
  b_dt <- b_dt[module %in% module_levels]
  b_dt <- b_dt[top_n %in% c(200, 500)]

  b_dt[, stability_fraction := fifelse(
    module == "NTM2\nASD-down",
    postnatal_top_fraction,
    late_prenatal_top_fraction
  )]

  b_dt[, stability_axis := fifelse(
    module == "NTM2\nASD-down",
    "Postnatal/adult-like",
    "Late-prenatal"
  )]

  b_dt[, module := factor(module, levels = module_levels)]
  b_dt[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]

  fwrite(b_dt, file.path(PLOTDIR, "S5B_stage_alignment_stability.tsv"), sep = "\t")

  panelB <- ggplot(b_dt, aes(x = cutoff, y = stability_fraction, group = module, colour = module)) +
    geom_hline(yintercept = 0.5, colour = COL_LIGHT, linewidth = 0.35, linetype = "dashed") +
    geom_line(linewidth = 0.55, alpha = 0.9) +
    geom_point(size = 2.6) +
    geom_text(
      aes(label = sprintf("%.2f", stability_fraction)),
      vjust = -0.70,
      size = 2.35,
      colour = COL_TEXT,
      family = font_family,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = module_cols, name = NULL) +
    scale_y_continuous(limits = c(0, 1.10), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Stage-call stability fraction", tag = "B") +
    theme_pub(7.8) +
    theme(
      axis.text.x = element_text(size = 7.2),
      legend.position = "bottom"
    )
}

if (is.null(panelB)) {
  panelB <- make_placeholder("B", "Stage-alignment stability\ncould not be generated.\nSee diagnostics.")
}

## ============================================================
## Panel C: topology-aware degree-matched null
## ============================================================

panelC <- NULL

if (!is.null(step12b) && nrow(step12b) > 0) {
  c_dt <- copy(step12b)

  c_dt[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]
  c_dt[, analysis_label := fifelse(
    grepl("SFARI", analysis, ignore.case = TRUE),
    "NTM2 SFARI degree null",
    "Transport degree null"
  )]
  c_dt[, analysis_label := factor(
    analysis_label,
    levels = c("Transport degree null", "NTM2 SFARI degree null")
  )]
  c_dt[, neglog_fdr := -log10(pmax(primary_fdr, 1e-300))]

  fwrite(c_dt, file.path(PLOTDIR, "S5C_degree_matched_null.tsv"), sep = "\t")

  panelC <- ggplot(c_dt, aes(x = cutoff, y = analysis_label)) +
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
    scale_shape_manual(
      values = c(
        "Transport degree null" = 21,
        "NTM2 SFARI degree null" = 24
      ),
      name = NULL
    ) +
    labs(x = NULL, y = NULL, tag = "C") +
    theme_pub(7.8) +
    theme(
      axis.text.x = element_text(size = 7.2),
      axis.text.y = element_text(size = 7.2),
      legend.position = "right"
    )
}

if (is.null(panelC)) {
  panelC <- make_placeholder("C", "Degree-matched topology null\ncould not be generated.\nSee diagnostics.")
}

## ============================================================
## Panel D: cross-disease boundary concordance
## ============================================================

panelD <- NULL

if (!is.null(disease_summary) && nrow(disease_summary) > 0) {
  d_dt <- copy(disease_summary)

  disease_col <- pick_col(d_dt, c("^disease$", "^dx$", "diagnosis", "disorder"))
  if (is.na(disease_col)) {
    d_dt[, disease := NA_character_]
  } else {
    setnames(d_dt, disease_col, "disease")
  }

  d_dt <- d_dt[disease %in% c("SCZ", "BD", "MDD")]

  if ("positive_rate" %in% names(d_dt)) {
    d_dt[, concordance := positive_rate]
  } else {
    pos_col <- pick_col(d_dt, c("n_positive", "positive"), numeric_only = TRUE)
    n_col <- pick_col(d_dt, c("n_tests", "n_total", "total"), numeric_only = TRUE)
    d_dt[, concordance := as.numeric(get(pos_col)) / as.numeric(get(n_col))]
  }

  if (!"stouffer_fdr_across_diseases" %in% names(d_dt)) {
    fdr_col <- pick_col(d_dt, c("stouffer.*fdr", "fdr", "q_value", "padj"), numeric_only = TRUE)
    if (!is.na(fdr_col)) setnames(d_dt, fdr_col, "stouffer_fdr_across_diseases")
  }

  if (!"n_positive" %in% names(d_dt)) {
    pos_col <- pick_col(d_dt, c("n_positive", "positive"), numeric_only = TRUE)
    if (!is.na(pos_col)) setnames(d_dt, pos_col, "n_positive")
  }

  if (!"n_tests" %in% names(d_dt)) {
    n_col <- pick_col(d_dt, c("n_tests", "n_total", "total"), numeric_only = TRUE)
    if (!is.na(n_col)) setnames(d_dt, n_col, "n_tests")
  }

  d_dt[, disease := factor(disease, levels = c("SCZ", "BD", "MDD"))]
  d_dt[, neglog_fdr := -log10(pmax(stouffer_fdr_across_diseases, 1e-300))]
  d_dt[, label := paste0(n_positive, "/", n_tests)]
  d_dt[, boundary := disease == "MDD"]

  fwrite(d_dt, file.path(PLOTDIR, "S5D_cross_disease_boundary_concordance.tsv"), sep = "\t")

  panelD <- ggplot(d_dt, aes(x = disease, y = concordance)) +
    geom_hline(yintercept = 0.5, colour = COL_LIGHT, linewidth = 0.35, linetype = "dashed") +
    geom_segment(
      aes(xend = disease, y = 0.5, yend = concordance),
      colour = COL_LIGHT,
      linewidth = 0.65
    ) +
    geom_point(
      aes(size = neglog_fdr, fill = boundary),
      shape = 21,
      colour = "white",
      stroke = 0.35
    ) +
    geom_text(
      aes(label = label),
      vjust = -1.1,
      size = 2.55,
      colour = COL_TEXT,
      family = font_family
    ) +
    scale_fill_manual(values = c("FALSE" = COL_GREEN, "TRUE" = COL_GREY), guide = "none") +
    scale_size_continuous(range = c(2.6, 6.4), name = expression(-log[10](FDR))) +
    scale_y_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Directional concordance", tag = "D") +
    theme_pub(7.8) +
    theme(
      axis.text.x = element_text(size = 7.3),
      legend.position = "right"
    )
}

if (is.null(panelD)) {
  panelD <- make_placeholder("D", "Cross-disease boundary concordance\ncould not be generated.\nSee diagnostics.")
}

## ============================================================
## Diagnostics
## ============================================================

diag <- c(
  "Supplementary Figure S5 v3 diagnostics",
  paste("Generated:", as.character(Sys.time())),
  "",
  paste("STEP12A_FILE:", STEP12A_FILE),
  paste("STEP12B_FILE:", STEP12B_FILE),
  paste("STEP12C_FILE:", STEP12C_FILE),
  paste("STEP12D_FILE:", STEP12D_FILE),
  paste("DISEASE_SUMMARY_FILE:", DISEASE_SUMMARY_FILE),
  paste("MODULE_SUMMARY_FILE:", MODULE_SUMMARY_FILE),
  "",
  paste("step12a rows:", if (is.null(step12a)) "NA" else nrow(step12a)),
  paste("step12b rows:", if (is.null(step12b)) "NA" else nrow(step12b)),
  paste("step12c rows:", if (is.null(step12c)) "NA" else nrow(step12c)),
  paste("step12d rows:", if (is.null(step12d)) "NA" else nrow(step12d)),
  paste("disease_summary rows:", if (is.null(disease_summary)) "NA" else nrow(disease_summary)),
  paste("module_summary rows:", if (is.null(module_summary)) "NA" else nrow(module_summary)),
  "",
  "step12a columns:",
  if (is.null(step12a)) "NULL" else paste(names(step12a), collapse = ", "),
  "",
  "step12b columns:",
  if (is.null(step12b)) "NULL" else paste(names(step12b), collapse = ", "),
  "",
  "step12c columns:",
  if (is.null(step12c)) "NULL" else paste(names(step12c), collapse = ", "),
  "",
  "step12d columns:",
  if (is.null(step12d)) "NULL" else paste(names(step12d), collapse = ", "),
  "",
  "disease_summary columns:",
  if (is.null(disease_summary)) "NULL" else paste(names(disease_summary), collapse = ", "),
  "",
  "module_summary columns:",
  if (is.null(module_summary)) "NULL" else paste(names(module_summary), collapse = ", "),
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
    widths = unit(c(1.05, 1.05), "null"),
    heights = unit(c(1.00, 1.00), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 10.6, height = 8.0, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 3180, height = 2400, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Supplementary Figure S5 v3 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plotdata:\n")
cat(PLOTDIR, "\n")
cat("Diagnostics:\n")
cat(DIAG_OUT, "\n")
