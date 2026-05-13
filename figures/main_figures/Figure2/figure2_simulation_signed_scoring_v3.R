suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Figure2")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

CONFIG_FILE <- file.path(ROOT, "step02_simulation_benchmark/results/11_step02B_hard_simulation_config.tsv")
SIM_FILE    <- file.path(ROOT, "step02_simulation_benchmark/results/13_step02B_summary_by_condition_method.tsv")
ABL_FILE    <- file.path(ROOT, "step02_simulation_benchmark/results/22_step02C_ablation_summary_by_component.tsv")
MODEL_FILE  <- file.path(ROOT, "step08_algorithm_strengthening/results/48_step08B_realdata_baseline_score_models.tsv")

PDF_OUT <- file.path(OUTDIR, "Figure2_simulation_signed_scoring_v3.pdf")
PNG_OUT <- file.path(OUTDIR, "Figure2_simulation_signed_scoring_v3.png")

PLOTDATA_DIR <- file.path(OUTDIR, "plotdata_v3")
dir.create(PLOTDATA_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------- style ----------
COL_TEXT   <- "#263238"
COL_MUTED  <- "#6B7280"
COL_LIGHT  <- "#D7DEE6"
COL_GRID   <- "#E8EDF2"

COL_GREEN  <- "#4F9A73"
COL_BLUE   <- "#4777A8"
COL_ORANGE <- "#D96C4A"
COL_PURPLE <- "#7A5AA6"
COL_YELLOW <- "#D8B34C"
COL_RED    <- "#B64A4A"
COL_GREY   <- "#8A96A3"

font_family <- "Helvetica"

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.30, colour = COL_TEXT),
      axis.text = element_text(colour = COL_TEXT, size = base_size),
      axis.title = element_text(colour = COL_TEXT, size = base_size + 1),
      strip.background = element_rect(fill = "white", colour = COL_LIGHT, linewidth = 0.35),
      strip.text = element_text(face = "bold", colour = COL_TEXT, size = base_size),
      legend.title = element_blank(),
      legend.text = element_text(colour = COL_TEXT, size = base_size - 0.5),
      legend.key.size = unit(0.34, "cm"),
      plot.tag = element_text(face = "bold", size = 16, colour = "black"),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(6, 7, 6, 7)
    )
}

rescale01 <- function(x, to = c(0.18, 1)) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    return(rep(mean(to), length(x)))
  }
  (x - rng[1]) / diff(rng) * diff(to) + to[1]
}

pretty_method <- function(x) {
  y <- as.character(x)
  y[y == "random"] <- "Random"
  y[y == "mean_signature"] <- "Mean signature"
  y[y == "embedding_nn"] <- "Embedding NN"
  y[y == "simple_embedding_nn"] <- "Embedding NN"
  y[y == "simple_ot"] <- "Simple OT"
  y[y == "graph_ot"] <- "Graph OT"
  y[y == "neurotrace_no_embedding"] <- "No embedding"
  y[y == "neurotrace_no_graph"] <- "No graph"
  y[y == "neurotrace_no_risk"] <- "No risk"
  y[y == "neurotrace_no_invariance"] <- "No invariance"
  y[y == "neurotrace_full"] <- "NeuroTRACE"
  y[y == "neurotrace_lite"] <- "NeuroTRACE"
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

## ============================================================
## A. Benchmark design matrix
## ============================================================

config <- fread(CONFIG_FILE)

cond_levels_raw <- c(
  "weak_signal",
  "low_overlap",
  "composition_confounded",
  "high_dropout",
  "false_prior_stress",
  "combined_hard"
)

cond_label_map <- c(
  weak_signal = "Weak\nsignal",
  low_overlap = "Low\noverlap",
  composition_confounded = "Composition\nconfounded",
  high_dropout = "High\ndropout",
  false_prior_stress = "False-prior\nstress",
  combined_hard = "Combined\nhard"
)

design_long <- rbindlist(list(
  data.table(
    condition = config$condition,
    feature = "Weak signal",
    value = rescale01(max(config$signal_strength, na.rm = TRUE) - config$signal_strength)
  ),
  data.table(
    condition = config$condition,
    feature = "Low state overlap",
    value = rescale01(max(config$state_overlap, na.rm = TRUE) - config$state_overlap)
  ),
  data.table(
    condition = config$condition,
    feature = "Composition shift",
    value = rescale01(config$composition_sd + config$adversarial_composition)
  ),
  data.table(
    condition = config$condition,
    feature = "Dropout",
    value = rescale01(config$dropout_rate)
  ),
  data.table(
    condition = config$condition,
    feature = "False prior",
    value = rescale01(config$false_prior_rate)
  )
), fill = TRUE)

design_long[, condition_label := cond_label_map[condition]]
design_long[, condition_label := factor(condition_label, levels = cond_label_map[cond_levels_raw])]
design_long[, feature := factor(feature, levels = rev(c(
  "Weak signal",
  "Low state overlap",
  "Composition shift",
  "Dropout",
  "False prior"
)))]

fwrite(design_long, file.path(PLOTDATA_DIR, "Figure2A_design_matrix_plotdata.tsv"), sep = "\t")

panelA <- ggplot(design_long, aes(x = condition_label, y = feature, fill = value)) +
  geom_tile(width = 0.91, height = 0.78, colour = "white", linewidth = 0.45) +
  scale_fill_gradient(
    low = "#EEF2F6",
    high = COL_BLUE,
    limits = c(0, 1),
    name = "Stress\nintensity"
  ) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(size = 7.4, angle = 0, hjust = 0.5, vjust = 1),
    axis.text.y = element_text(size = 7.8),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7.4, colour = COL_TEXT),
    legend.text = element_text(size = 7.0, colour = COL_TEXT),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.22, "cm")
  )

## ============================================================
## B. Hard-simulation performance across key methods
## ============================================================

sim <- fread(SIM_FILE)

selected_methods <- c(
  "neurotrace_full",
  "graph_ot",
  "simple_ot",
  "mean_signature",
  "random"
)

sim <- sim[method %in% selected_methods]

metric_long <- melt(
  sim,
  id.vars = c("condition", "method"),
  measure.vars = c(
    "alignment_top1_mean",
    "true_transport_mass_mean",
    "gene_auprc_mean",
    "precision_at_50_mean"
  ),
  variable.name = "metric",
  value.name = "value"
)

metric_long[, method_label := pretty_method(method)]
metric_long[, method_label := factor(method_label, levels = pretty_method(selected_methods))]

metric_long[, metric_label := fifelse(
  metric == "alignment_top1_mean", "Top-1\nalignment",
  fifelse(
    metric == "true_transport_mass_mean", "True transport\nmass",
    fifelse(
      metric == "gene_auprc_mean", "Gene\nAUPRC",
      "Precision\nat 50"
    )
  )
)]

metric_long[, metric_label := factor(metric_label, levels = c(
  "Top-1\nalignment",
  "True transport\nmass",
  "Gene\nAUPRC",
  "Precision\nat 50"
))]

main_metric <- metric_long[
  condition != "combined_hard",
  .(
    value = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = .N
  ),
  by = .(method_label, metric_label)
]
main_metric[, se := sd / sqrt(n)]

fwrite(metric_long, file.path(PLOTDATA_DIR, "Figure2B_hard_simulation_all_conditions_plotdata.tsv"), sep = "\t")
fwrite(main_metric, file.path(PLOTDATA_DIR, "Figure2B_key_method_average_plotdata.tsv"), sep = "\t")

method_cols <- c(
  "NeuroTRACE" = COL_GREEN,
  "Graph OT" = COL_BLUE,
  "Simple OT" = COL_ORANGE,
  "Mean signature" = COL_YELLOW,
  "Random" = COL_GREY
)

panelB <- ggplot(main_metric, aes(x = method_label, y = value, fill = method_label)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.18) +
  geom_errorbar(
    aes(ymin = pmax(0, value - se), ymax = pmin(1, value + se)),
    width = 0.18,
    linewidth = 0.25,
    colour = COL_TEXT
  ) +
  facet_wrap(~ metric_label, nrow = 1) +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(x = NULL, y = "Mean performance", tag = "B") +
  theme_pub(7.7) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.30, "cm"),
    panel.spacing = unit(0.55, "lines")
  )

## ============================================================
## C. Component-ablation loss
## ============================================================

abl <- fread(ABL_FILE)

abl[, component := pretty_method(method)]

component_levels <- c("No embedding", "No graph", "No risk", "No invariance")
abl <- abl[component %in% component_levels]
abl[, component := factor(component, levels = component_levels)]

abl_long <- melt(
  abl,
  id.vars = c("method", "component"),
  measure.vars = c(
    "mean_loss_alignment",
    "mean_loss_transport",
    "mean_decoy_gain_when_removed"
  ),
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

abl_long[, quantity_label := factor(quantity_label, levels = c(
  "Alignment loss",
  "Transport loss",
  "Decoy gain"
))]

fwrite(abl_long, file.path(PLOTDATA_DIR, "Figure2C_ablation_loss_plotdata.tsv"), sep = "\t")

panelC <- ggplot(abl_long, aes(x = value, y = component, colour = quantity_label)) +
  geom_vline(xintercept = 0, colour = COL_LIGHT, linewidth = 0.4) +
  geom_segment(
    aes(x = 0, xend = value, yend = component),
    linewidth = 0.65,
    alpha = 0.85
  ) +
  geom_point(size = 2.35) +
  scale_colour_manual(values = c(
    "Alignment loss" = COL_GREEN,
    "Transport loss" = COL_BLUE,
    "Decoy gain" = COL_RED
  )) +
  labs(x = "Change after component removal", y = NULL, tag = "C") +
  theme_pub(8.2) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8.2),
    axis.title.x = element_text(size = 8.4)
  )

## ============================================================
## D. Real-data NTM2 signed vs unsigned scoring forest plot
## ============================================================

model <- fread(MODEL_FILE)

d <- model[
  program == "NTM2_ASD_down" &
    top_n == max(model[program == "NTM2_ASD_down"]$top_n, na.rm = TRUE)
]

d[, method_label := pretty_score(scoring_method)]

method_levels_d <- c(
  "Abs-weight\nunsigned",
  "Unsigned\nmean",
  "Sign only",
  "NeuroTRACE\nweighted signed"
)
d[, method_label := factor(method_label, levels = method_levels_d)]

d[, cohort := factor(cohort, levels = c("GSE102741", "GSE64018"))]
d[, ci_low := beta_ASD_vs_Control - 1.96 * se]
d[, ci_high := beta_ASD_vs_Control + 1.96 * se]

xlim_d <- range(c(d$ci_low, d$ci_high), na.rm = TRUE)
xpad <- diff(xlim_d) * 0.08
xlim_d <- c(xlim_d[1] - xpad, xlim_d[2] + xpad)

fwrite(d, file.path(PLOTDATA_DIR, "Figure2D_NTM2_realdata_scoring_plotdata.tsv"), sep = "\t")

score_cols <- c(
  "NeuroTRACE\nweighted signed" = COL_GREEN,
  "Sign only" = COL_BLUE,
  "Unsigned\nmean" = COL_ORANGE,
  "Abs-weight\nunsigned" = COL_PURPLE
)

panelD <- ggplot(d, aes(x = beta_ASD_vs_Control, y = method_label, colour = method_label)) +
  geom_vline(xintercept = 0, colour = COL_LIGHT, linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, linewidth = 0.38) +
  geom_point(size = 2.45) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = score_cols, drop = FALSE) +
  coord_cartesian(xlim = xlim_d) +
  labs(x = "ASD effect on NTM2 score", y = NULL, tag = "D") +
  theme_pub(8.2) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7.9),
    panel.spacing = unit(0.8, "lines")
  )

## ============================================================
## Export combined figure
## ============================================================

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 2,
    widths = unit(c(0.95, 1.35), "null"),
    heights = unit(c(0.95, 1.05), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 9.4, height = 7.0, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 2820, height = 2100, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Figure 2 v3 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plot data written to:\n")
cat(PLOTDATA_DIR, "\n")
