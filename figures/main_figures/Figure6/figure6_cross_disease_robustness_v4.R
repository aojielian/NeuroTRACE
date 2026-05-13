suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Figure6")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

CROSS_MODEL_FILE <- file.path(ROOT, "step08_algorithm_strengthening/results/82_step08E_cross_disease_NTM_models.tsv")
DISEASE_SUMMARY_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/02_step11E_v2_disease_level_summary.tsv")
MODULE_SUMMARY_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/03_step11E_v2_module_level_summary.tsv")

PDF_OUT <- file.path(OUTDIR, "Figure6_cross_disease_robustness_v4.pdf")
PNG_OUT <- file.path(OUTDIR, "Figure6_cross_disease_robustness_v4.png")

PLOTDATA_DIR <- file.path(OUTDIR, "plotdata_v4")
dir.create(PLOTDATA_DIR, recursive = TRUE, showWarnings = FALSE)

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
      axis.title = element_text(colour = COL_TEXT, size = base_size + 1),
      strip.background = element_rect(fill = "white", colour = COL_LIGHT, linewidth = 0.35),
      strip.text = element_text(face = "bold", colour = COL_TEXT, size = base_size),
      legend.title = element_text(colour = COL_TEXT, size = base_size - 0.2),
      legend.text = element_text(colour = COL_TEXT, size = base_size - 0.7),
      legend.key.size = unit(0.33, "cm"),
      plot.tag = element_text(face = "bold", size = 16, colour = "black"),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(8, 9, 10, 9)
    )
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

disease_levels <- c("SCZ", "BD", "MDD")

## ============================================================
## A. Cross-disease NTM effects heatmap, top500
## ============================================================

module_summary <- fread(MODULE_SUMMARY_FILE)

a_dt <- module_summary[top_n == 500 & disease %in% disease_levels]
a_dt[, module := pretty_module(module_family)]
a_dt <- a_dt[module %in% module_levels]

a_plot <- a_dt[
  ,
  .(
    beta = mean(beta, na.rm = TRUE),
    min_global_fdr = min(global_fdr, na.rm = TRUE),
    n_datasets = uniqueN(dataset)
  ),
  by = .(disease, module)
]

a_plot[, disease := factor(disease, levels = disease_levels)]
a_plot[, module := factor(module, levels = rev(module_levels))]
a_plot[, sig_label := fifelse(min_global_fdr < 0.05, "*", "")]

fwrite(a_plot, file.path(PLOTDATA_DIR, "Figure6A_cross_disease_effect_heatmap_top500.tsv"), sep = "\t")

lim_a <- max(abs(a_plot$beta), na.rm = TRUE)
if (!is.finite(lim_a) || lim_a == 0) lim_a <- 1

panelA <- ggplot(a_plot, aes(x = disease, y = module, fill = beta)) +
  geom_tile(colour = "white", linewidth = 0.45, width = 0.92, height = 0.78) +
  geom_text(aes(label = sig_label), colour = COL_TEXT, size = 3.2, family = font_family) +
  scale_fill_gradient2(
    low = COL_BLUE,
    mid = "#F4F6F8",
    high = COL_RED,
    midpoint = 0,
    limits = c(-lim_a, lim_a),
    name = "Mean\neffect"
  ) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub(8.2) +
  theme(
    axis.text.x = element_text(size = 8.0),
    axis.text.y = element_text(size = 7.8),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

## ============================================================
## B. Disease-level directional concordance, fixed label spacing
## ============================================================

disease_summary <- fread(DISEASE_SUMMARY_FILE)
b_dt <- disease_summary[disease %in% disease_levels]

b_dt[, disease := factor(disease, levels = disease_levels)]
b_dt[, concordance := positive_rate]
b_dt[, label := paste0(n_positive, "/", n_tests)]
b_dt[, neglog_stouffer_fdr := -log10(pmax(stouffer_fdr_across_diseases, 1e-300))]
b_dt[, boundary := disease == "MDD"]

## Put labels safely inside the plotting region.
b_dt[, label_y := fifelse(concordance <= 0.03, 0.08, pmin(1.08, concordance + 0.045))]

fwrite(b_dt, file.path(PLOTDATA_DIR, "Figure6B_disease_level_concordance_dot.tsv"), sep = "\t")

panelB <- ggplot(b_dt, aes(x = disease, y = concordance)) +
  geom_hline(yintercept = 0.5, colour = COL_LIGHT, linewidth = 0.4, linetype = "dashed") +
  geom_segment(
    aes(xend = disease, y = 0.5, yend = concordance),
    colour = COL_LIGHT,
    linewidth = 0.7
  ) +
  geom_point(
    aes(size = neglog_stouffer_fdr, fill = boundary),
    shape = 21,
    colour = "white",
    stroke = 0.35
  ) +
  geom_text(
    aes(y = label_y, label = label),
    size = 2.9,
    colour = COL_TEXT,
    family = font_family
  ) +
  scale_fill_manual(values = c("TRUE" = COL_GREY, "FALSE" = COL_GREEN), guide = "none") +
  scale_size_continuous(range = c(2.8, 6.5), name = "-log10\nStouffer FDR") +
  scale_y_continuous(
    limits = c(-0.07, 1.15),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Directional concordance", tag = "B") +
  coord_cartesian(clip = "off") +
  theme_pub(8.2) +
  theme(
    axis.text.x = element_text(size = 8.0),
    legend.position = "right",
    plot.margin = margin(8, 14, 18, 10)
  )

## ============================================================
## C. MDD available-dataset boundary forest plot, top500
## ============================================================

cross_model <- fread(CROSS_MODEL_FILE)

c_dt <- cross_model[target_dx == "MDD" & top_n == 500]
c_dt[, module := pretty_module(program)]
c_dt <- c_dt[module %in% module_levels]
c_dt[, module := factor(module, levels = rev(module_levels))]
c_dt[, ci_low := beta_target_vs_Control - 1.96 * se]
c_dt[, ci_high := beta_target_vs_Control + 1.96 * se]

fwrite(c_dt, file.path(PLOTDATA_DIR, "Figure6C_MDD_boundary_top500.tsv"), sep = "\t")

panelC <- ggplot(c_dt, aes(x = beta_target_vs_Control, y = module, colour = module)) +
  geom_vline(xintercept = 0, colour = COL_LIGHT, linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    height = 0.16,
    linewidth = 0.38
  ) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = module_cols, drop = FALSE) +
  labs(x = "MDD effect on NTM score", y = NULL, tag = "C") +
  theme_pub(8.2) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8.0),
    plot.margin = margin(8, 12, 10, 10)
  )

## ============================================================
## Export three-panel figure
## ============================================================

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 2,
    widths = unit(c(1.05, 1.15), "null"),
    heights = unit(c(1.00, 0.90), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:2))

  popViewport()
}

pdf(PDF_OUT, width = 9.6, height = 6.4, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 2880, height = 1920, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Figure 6 v4 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plot data written to:\n")
cat(PLOTDATA_DIR, "\n")
