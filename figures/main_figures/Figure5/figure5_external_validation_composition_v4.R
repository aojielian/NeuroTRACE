suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Figure5")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EXT_MODEL_FILE <- file.path(ROOT, "step06_cross_cohort_validation/results/35_step06B_external_NTM_score_models.tsv")
RAND_SUMMARY_FILE <- file.path(ROOT, "step08_algorithm_strengthening/results/44_step08A_exprvar_matched_random_specificity_summary.tsv")
RAND_NULL_FILE <- file.path(ROOT, "step08_algorithm_strengthening/results/43_step08A_exprvar_matched_random_null_betas.tsv.gz")
COMPOSITION_MODEL_FILE <- file.path(ROOT, "step10_composition_adjustment/results/08_step10E_NTM_models_before_after_composition.tsv")
NNLS_MODEL_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/06_step11H_Gandal_NNLS_model_comparison_summary.tsv")

PDF_OUT <- file.path(OUTDIR, "Figure5_external_validation_composition_v4.pdf")
PNG_OUT <- file.path(OUTDIR, "Figure5_external_validation_composition_v4.png")
PLOTDATA_DIR <- file.path(OUTDIR, "plotdata_v4")
dir.create(PLOTDATA_DIR, recursive = TRUE, showWarnings = FALSE)

COL_TEXT   <- "#263238"
COL_MUTED  <- "#6B7280"
COL_LIGHT  <- "#D7DEE6"
COL_GREEN  <- "#4F9A73"
COL_BLUE   <- "#4777A8"
COL_ORANGE <- "#D96C4A"
COL_PURPLE <- "#7A5AA6"

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
      legend.text = element_text(colour = COL_TEXT, size = base_size - 0.6),
      legend.key.size = unit(0.32, "cm"),
      plot.tag = element_text(face = "bold", size = 16, colour = "black"),
      plot.tag.position = c(0.015, 0.985),
      plot.margin = margin(6, 7, 6, 7)
    )
}

read_dt <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("zcat", shQuote(path)))
  } else {
    fread(path)
  }
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

## ============================================================
## A. External ASD validation, top500 only
## ============================================================

ext <- fread(EXT_MODEL_FILE)

a_dt <- ext[top_n == 500]
a_dt[, module := pretty_module(program)]
a_dt[, module := factor(module, levels = rev(module_levels))]
a_dt[, cohort := factor(cohort, levels = c("GSE102741", "GSE64018"))]
a_dt[, ci_low := beta_ASD_vs_Control - 1.96 * se]
a_dt[, ci_high := beta_ASD_vs_Control + 1.96 * se]

fwrite(a_dt, file.path(PLOTDATA_DIR, "Figure5A_external_ASD_validation_top500.tsv"), sep = "\t")

panelA <- ggplot(a_dt, aes(x = beta_ASD_vs_Control, y = module, colour = module)) +
  geom_vline(xintercept = 0, colour = COL_LIGHT, linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    height = 0.16,
    linewidth = 0.38
  ) +
  geom_point(size = 2.45) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = module_cols, drop = FALSE) +
  labs(x = "ASD effect on NTM score", y = NULL, tag = "A") +
  theme_pub(8.0) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7.7),
    panel.spacing = unit(0.8, "lines")
  )

## ============================================================
## B. Expression/variance-matched random specificity, top500 only
## ============================================================

rand_null <- read_dt(RAND_NULL_FILE)
rand_sum <- fread(RAND_SUMMARY_FILE)

b_null <- rand_null[top_n == 500]
b_null[, module := pretty_module(program)]
b_null <- b_null[module %in% module_levels]
b_null[, module := factor(module, levels = rev(module_levels))]
b_null[, cohort := factor(cohort, levels = c("GSE102741", "GSE64018"))]

b_obs <- rand_sum[top_n == 500]
b_obs[, module := pretty_module(program)]
b_obs <- b_obs[module %in% module_levels]
b_obs[, module := factor(module, levels = rev(module_levels))]
b_obs[, cohort := factor(cohort, levels = c("GSE102741", "GSE64018"))]

fwrite(b_null, file.path(PLOTDATA_DIR, "Figure5B_matched_random_null_top500.tsv"), sep = "\t")
fwrite(b_obs, file.path(PLOTDATA_DIR, "Figure5B_observed_effects_top500.tsv"), sep = "\t")

panelB <- ggplot(b_null, aes(x = null_beta, y = module)) +
  geom_vline(xintercept = 0, colour = COL_LIGHT, linewidth = 0.45) +
  geom_violin(
    fill = "#DDE6EE",
    colour = COL_MUTED,
    linewidth = 0.32,
    scale = "width",
    trim = TRUE
  ) +
  geom_point(
    data = b_obs,
    aes(x = observed_beta, y = module),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.65,
    fill = COL_GREEN,
    colour = "white",
    stroke = 0.25
  ) +
  facet_wrap(~ cohort, nrow = 1) +
  labs(x = "ASD effect relative to matched random modules", y = NULL, tag = "B") +
  theme_pub(8.0) +
  theme(
    axis.text.y = element_text(size = 7.7),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  )

## ============================================================
## C. Marker-adjustment log2 attenuation ratio, top500
## ============================================================

comp <- fread(COMPOSITION_MODEL_FILE)

c_raw <- comp[grepl("_top500$", ntm_score)]
c_raw[, program := gsub("_top500$", "", ntm_score)]
c_raw[, module := pretty_module(program)]
c_raw <- c_raw[module %in% module_levels]

## Match base and adjusted effects within the same cohort and NTM.
base_dt <- c_raw[model_type == "M0_base_covariates", .(
  cohort,
  module,
  base_beta = beta
)]

adj_dt <- c_raw[model_type %in% c("M1_marker_PC1_PC2", "M2_all_marker_scores"), .(
  cohort,
  module,
  model_type,
  beta
)]

adj_dt[, adjust_model := fifelse(
  model_type == "M1_marker_PC1_PC2",
  "Marker PCs",
  "All marker scores"
)]

c_per_cohort <- merge(
  adj_dt,
  base_dt,
  by = c("cohort", "module"),
  all.x = TRUE,
  allow.cartesian = FALSE
)

## Ratio is interpretable only when base and adjusted beta have the same positive sign.
## For this dataset, the effects are direction-preserved; keep a defensive filter.
c_per_cohort <- c_per_cohort[
  !is.na(base_beta) &
    !is.na(beta) &
    base_beta != 0 &
    (beta / base_beta) > 0
]

c_per_cohort[, ratio := beta / base_beta]
c_per_cohort[, log2_ratio := log2(ratio)]

## Collapse to one clear point per NTM × adjustment model.
c_dt <- c_per_cohort[
  ,
  .(
    log2_ratio = median(log2_ratio, na.rm = TRUE),
    log2_min = min(log2_ratio, na.rm = TRUE),
    log2_max = max(log2_ratio, na.rm = TRUE),
    ratio_median = median(ratio, na.rm = TRUE),
    n_cohorts = uniqueN(cohort)
  ),
  by = .(module, adjust_model)
]

c_dt[, module := factor(module, levels = module_levels)]
c_dt[, adjust_model := factor(adjust_model, levels = c("Marker PCs", "All marker scores"))]

fwrite(
  c_per_cohort,
  file.path(PLOTDATA_DIR, "Figure5C_marker_adjustment_log2_ratio_per_cohort_top500.tsv"),
  sep = "\t"
)
fwrite(
  c_dt,
  file.path(PLOTDATA_DIR, "Figure5C_marker_adjustment_log2_ratio_top500.tsv"),
  sep = "\t"
)

adj_cols <- c(
  "Marker PCs" = COL_GREEN,
  "All marker scores" = COL_PURPLE
)

pd <- position_dodge(width = 0.45)

c_dt[, seg_ymin := pmin(0, log2_ratio)]
c_dt[, seg_ymax := pmax(0, log2_ratio)]

panelC <- ggplot(c_dt, aes(x = module, y = log2_ratio, colour = adjust_model, shape = adjust_model)) +
  geom_hline(yintercept = 0, colour = COL_LIGHT, linewidth = 0.45, linetype = "dashed") +
  geom_linerange(
    aes(ymin = seg_ymin, ymax = seg_ymax),
    position = pd,
    linewidth = 0.55,
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(ymin = log2_min, ymax = log2_max),
    position = pd,
    width = 0.12,
    linewidth = 0.32,
    alpha = 0.75
  ) +
  geom_point(
    position = pd,
    size = 2.8,
    stroke = 0.25
  ) +
  scale_colour_manual(values = adj_cols, drop = FALSE) +
  scale_shape_manual(values = c("Marker PCs" = 16, "All marker scores" = 17), drop = FALSE) +
  labs(x = NULL, y = "log2(adjusted / base beta ratio)", tag = "C") +
  theme_pub(8.0) +
  theme(
    axis.text.x = element_text(size = 7.6),
    legend.position = "bottom"
  )


## ============================================================
## D. NNLS deconvolution adjustment, top500
## ============================================================

nnls <- fread(NNLS_MODEL_FILE)

d_dt <- nnls[top_n == 500]
d_dt[, module := pretty_module(program)]
d_dt <- d_dt[module %in% module_levels]

d_long <- melt(
  d_dt,
  id.vars = c("dataset", "program", "module", "top_n"),
  measure.vars = c("base_beta", "NNLS_PC12_beta", "NNLS_props_beta"),
  variable.name = "model_type",
  value.name = "beta"
)

d_long[, model_label := fifelse(
  model_type == "base_beta", "Base model",
  fifelse(
    model_type == "NNLS_PC12_beta", "NNLS-PC adjusted",
    "NNLS-proportion adjusted"
  )
)]

d_long[, model_label := factor(
  model_label,
  levels = c("Base model", "NNLS-PC adjusted", "NNLS-proportion adjusted")
)]
d_long[, module := factor(module, levels = module_levels)]

fwrite(d_long, file.path(PLOTDATA_DIR, "Figure5D_NNLS_adjustment_top500.tsv"), sep = "\t")

panelD <- ggplot(d_long, aes(x = model_label, y = beta, group = module, colour = module)) +
  geom_hline(yintercept = 0, colour = COL_LIGHT, linewidth = 0.45) +
  geom_line(linewidth = 0.55, alpha = 0.9) +
  geom_point(size = 2.35) +
  scale_colour_manual(values = module_cols, drop = FALSE) +
  labs(x = NULL, y = "ASD effect after NNLS adjustment", tag = "D") +
  theme_pub(8.0) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 7.3),
    legend.position = "bottom"
  )

## ============================================================
## Export combined figure
## ============================================================

draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 2,
    widths = unit(c(1.10, 1.12), "null"),
    heights = unit(c(1.00, 1.04), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 9.4, height = 7.1, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 2820, height = 2130, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Figure 5 v4 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plot data written to:\n")
cat(PLOTDATA_DIR, "\n")
