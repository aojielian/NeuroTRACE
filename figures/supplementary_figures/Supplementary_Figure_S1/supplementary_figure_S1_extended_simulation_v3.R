suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Supplementary_FigureS1")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

SIM_SUMMARY_FILE <- file.path(ROOT, "step02_simulation_benchmark/results/13_step02B_summary_by_condition_method.tsv")
ABLATION_FILE <- file.path(ROOT, "step02_simulation_benchmark/results/22_step02C_ablation_summary_by_component.tsv")
SCORE_MODEL_FILE <- file.path(ROOT, "step08_algorithm_strengthening/results/48_step08B_realdata_baseline_score_models.tsv")

PDF_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S1_v3.pdf")
PNG_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S1_v3.png")
DIAG_OUT <- file.path(OUTDIR, "Supplementary_Figure_S1_v3_diagnostics.txt")
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

pretty_condition <- function(x) {
  y <- as.character(x)
  out <- y
  out[y == "weak_signal"] <- "Weak signal"
  out[y == "low_overlap"] <- "Low overlap"
  out[y == "composition_confounded"] <- "Composition confounded"
  out[y == "high_dropout"] <- "High dropout"
  out[y == "false_prior_stress"] <- "False-prior stress"
  out[y == "combined_hard"] <- "Combined hard"

  ## Fallback for already-human-readable labels
  out <- gsub("_", " ", out)
  out <- gsub("false prior stress", "False-prior stress", out, ignore.case = TRUE)
  out <- gsub("weak signal", "Weak signal", out, ignore.case = TRUE)
  out <- gsub("low overlap", "Low overlap", out, ignore.case = TRUE)
  out <- gsub("composition confounded", "Composition confounded", out, ignore.case = TRUE)
  out <- gsub("high dropout", "High dropout", out, ignore.case = TRUE)
  out <- gsub("combined hard", "Combined hard", out, ignore.case = TRUE)
  out
}

pretty_method <- function(x) {
  y <- as.character(x)
  y[y == "neurotrace_full"] <- "NeuroTRACE"
  y[y == "neurotrace_lite"] <- "NeuroTRACE"
  y[y == "graph_ot"] <- "Graph OT"
  y[y == "simple_ot"] <- "Simple OT"
  y[y == "embedding_nn"] <- "Embedding NN"
  y[y == "simple_embedding_nn"] <- "Embedding NN"
  y[y == "mean_signature"] <- "Mean signature"
  y[y == "random"] <- "Random"
  y[y == "neurotrace_no_embedding"] <- "No embedding"
  y[y == "neurotrace_no_graph"] <- "No graph"
  y[y == "neurotrace_no_risk"] <- "No risk"
  y[y == "neurotrace_no_invariance"] <- "No invariance"
  y
}

pretty_score <- function(x) {
  y <- as.character(x)
  y[y == "neurotrace_weighted_signed"] <- "NeuroTRACE\nweighted signed"
  y[y == "sign_only"] <- "Sign only"
  y[y == "unsigned_mean"] <- "Unsigned\nmean"
  y[y == "abs_weight_unsigned"] <- "Abs-weight\nunsigned"
  y
}

pretty_program <- function(x) {
  y <- as.character(x)
  y <- gsub("_top[0-9]+$", "", y)
  out <- y
  out[grepl("NTM1", y)] <- "NTM1\nASD-up"
  out[grepl("NTM2", y)] <- "NTM2\nASD-down"
  out[grepl("NTM3", y)] <- "NTM3\nASD-signed"
  out
}

## ============================================================
## Read inputs
## ============================================================

sim <- fread(SIM_SUMMARY_FILE)
abl <- fread(ABLATION_FILE)
score <- fread(SCORE_MODEL_FILE)

## ============================================================
## Panel A: condition-specific benchmark heatmap
## ============================================================

selected_methods_raw <- c(
  "neurotrace_full",
  "neurotrace_lite",
  "graph_ot",
  "simple_ot",
  "embedding_nn",
  "simple_embedding_nn",
  "mean_signature",
  "random"
)

metric_cols <- c(
  "alignment_top1_mean",
  "true_transport_mass_mean",
  "gene_auprc_mean",
  "precision_at_50_mean"
)

metric_cols <- metric_cols[metric_cols %in% names(sim)]

if (length(metric_cols) == 0) {
  stop("No expected simulation metric columns found in: ", SIM_SUMMARY_FILE)
}

sim_long <- melt(
  sim,
  id.vars = c("condition", "method"),
  measure.vars = metric_cols,
  variable.name = "metric",
  value.name = "value"
)

sim_long <- sim_long[method %in% selected_methods_raw]

sim_long[, method_label := pretty_method(method)]
sim_long[, condition_label := pretty_condition(condition)]

sim_long[, metric_label := fifelse(
  metric == "alignment_top1_mean", "Top-1 alignment",
  fifelse(
    metric == "true_transport_mass_mean", "True transport mass",
    fifelse(
      metric == "gene_auprc_mean", "Gene AUPRC",
      "Precision@50"
    )
  )
)]

method_levels <- c("NeuroTRACE", "Graph OT", "Simple OT", "Embedding NN", "Mean signature", "Random")
method_levels <- method_levels[method_levels %in% unique(sim_long$method_label)]

condition_levels <- c(
  "Weak signal",
  "Low overlap",
  "Composition confounded",
  "High dropout",
  "False-prior stress",
  "Combined hard"
)
condition_levels <- condition_levels[condition_levels %in% unique(sim_long$condition_label)]

metric_levels <- c("Top-1 alignment", "True transport mass", "Gene AUPRC", "Precision@50")

sim_long[, method_label := factor(method_label, levels = rev(method_levels))]
sim_long[, condition_label := factor(condition_label, levels = condition_levels)]
sim_long[, metric_label := factor(metric_label, levels = metric_levels)]

fwrite(sim_long, file.path(PLOTDIR, "S1A_condition_specific_benchmark_heatmap.tsv"), sep = "\t")

panelA <- ggplot(sim_long, aes(x = condition_label, y = method_label, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_fill_gradient(
    low = "#EEF2F6",
    high = COL_GREEN,
    limits = c(0, 1),
    name = "Mean\nperformance"
  ) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub(7.6) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.6),
    axis.text.y = element_text(size = 7.2),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.spacing = unit(0.7, "lines")
  )

## ============================================================
## Panel B: method-rank summary
## ============================================================

rank_dt <- copy(sim_long)
rank_dt[, rank := frank(-value, ties.method = "average"), by = .(condition_label, metric_label)]

rank_sum <- rank_dt[
  ,
  .(
    mean_rank = mean(rank, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE)
  ),
  by = .(method_label, metric_label)
]

rank_sum[, method_label := factor(method_label, levels = rev(method_levels))]
rank_sum[, metric_label := factor(metric_label, levels = metric_levels)]

fwrite(rank_sum, file.path(PLOTDIR, "S1B_method_rank_summary.tsv"), sep = "\t")

panelB <- ggplot(rank_sum, aes(x = mean_rank, y = method_label)) +
  geom_vline(xintercept = 1, colour = COL_LIGHT, linewidth = 0.35, linetype = "dashed") +
  geom_point(aes(size = mean_value, colour = mean_rank)) +
  facet_wrap(~ metric_label, ncol = 2) +
  scale_colour_gradient(low = COL_GREEN, high = COL_RED, trans = "reverse", name = "Mean\nrank") +
  scale_size_continuous(range = c(2.2, 6.5), name = "Mean\nscore") +
  scale_x_reverse() +
  labs(x = "Average rank across hard conditions", y = NULL, tag = "B") +
  theme_pub(7.8) +
  theme(
    axis.text.y = element_text(size = 7.2),
    legend.position = "right",
    panel.spacing = unit(0.7, "lines")
  )

## ============================================================
## Panel C: component ablation details
## ============================================================

abl[, component_label := pretty_method(method)]

abl_long <- melt(
  abl,
  id.vars = c("method", "component_label"),
  measure.vars = c("mean_loss_alignment", "mean_loss_transport", "mean_decoy_gain_when_removed"),
  variable.name = "quantity",
  value.name = "value"
)

abl_long[, quantity_label := fifelse(
  quantity == "mean_loss_alignment", "Alignment loss",
  fifelse(
    quantity == "mean_loss_transport", "Transport loss",
    "Decoy gain"
  )
)]

component_levels <- c("No graph", "No embedding", "No risk", "No invariance")
component_levels <- component_levels[component_levels %in% unique(abl_long$component_label)]

abl_long[, component_label := factor(component_label, levels = rev(component_levels))]
abl_long[, quantity_label := factor(quantity_label, levels = c("Alignment loss", "Transport loss", "Decoy gain"))]

fwrite(abl_long, file.path(PLOTDIR, "S1C_component_ablation_details.tsv"), sep = "\t")

panelC <- ggplot(abl_long, aes(x = value, y = component_label, fill = quantity_label)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Alignment loss" = COL_GREEN,
      "Transport loss" = COL_BLUE,
      "Decoy gain" = COL_RED
    ),
    name = NULL
  ) +
  labs(x = "Change after component removal", y = NULL, tag = "C") +
  theme_pub(8.0) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 7.4)
  )

## ============================================================
## Panel D: all-module signed vs unsigned scoring comparison
## ============================================================

score <- score[top_n == max(top_n, na.rm = TRUE)]

score[, program_label := pretty_program(program)]
score[, method_label := pretty_score(scoring_method)]
score[, cohort := factor(cohort, levels = c("GSE102741", "GSE64018"))]

score[, program_label := factor(
  program_label,
  levels = c("NTM1\nASD-up", "NTM2\nASD-down", "NTM3\nASD-signed")
)]

score[, method_label := factor(
  method_label,
  levels = c(
    "NeuroTRACE\nweighted signed",
    "Sign only",
    "Unsigned\nmean",
    "Abs-weight\nunsigned"
  )
)]

score[, ci_low := beta_ASD_vs_Control - 1.96 * se]
score[, ci_high := beta_ASD_vs_Control + 1.96 * se]

fwrite(score, file.path(PLOTDIR, "S1D_signed_unsigned_scoring_all_modules.tsv"), sep = "\t")

score_cols <- c(
  "NeuroTRACE\nweighted signed" = COL_GREEN,
  "Sign only" = COL_BLUE,
  "Unsigned\nmean" = COL_ORANGE,
  "Abs-weight\nunsigned" = COL_PURPLE
)

panelD <- ggplot(score, aes(x = method_label, y = beta_ASD_vs_Control, colour = method_label)) +
  geom_hline(yintercept = 0, colour = COL_LIGHT, linewidth = 0.35) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.16,
    linewidth = 0.28,
    alpha = 0.75
  ) +
  geom_point(size = 1.9) +
  facet_grid(program_label ~ cohort) +
  scale_colour_manual(values = score_cols, drop = FALSE) +
  labs(x = NULL, y = "ASD effect", tag = "D") +
  theme_pub(7.4) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.4),
    axis.text.y = element_text(size = 6.8),
    legend.position = "none",
    panel.spacing = unit(0.55, "lines")
  )

## ============================================================
## Diagnostics
## ============================================================

diag <- c(
  "Supplementary Figure S1 v3 diagnostics",
  paste("Generated:", as.character(Sys.time())),
  "",
  paste("SIM_SUMMARY_FILE:", SIM_SUMMARY_FILE),
  paste("ABLATION_FILE:", ABLATION_FILE),
  paste("SCORE_MODEL_FILE:", SCORE_MODEL_FILE),
  "",
  paste("sim rows:", nrow(sim)),
  paste("abl rows:", nrow(abl)),
  paste("score rows:", nrow(score)),
  "",
  "sim columns:",
  paste(names(sim), collapse = ", "),
  "",
  "abl columns:",
  paste(names(abl), collapse = ", "),
  "",
  "score columns:",
  paste(names(score), collapse = ", "),
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
    heights = unit(c(1.00, 1.05), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 10.4, height = 8.1, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 3120, height = 2430, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Supplementary Figure S1 v3 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plotdata:\n")
cat(PLOTDIR, "\n")
cat("Diagnostics:\n")
cat(DIAG_OUT, "\n")
