suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
IN_DIR <- file.path(ROOT, "paper_figures", "Figure3", "plotdata_v2")
OUTDIR <- file.path(ROOT, "paper_figures", "Figure3", "FROZEN")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PDF_OUT <- file.path(OUTDIR, "Figure3_OFFICIAL_v1.pdf")
PNG_OUT <- file.path(OUTDIR, "Figure3_OFFICIAL_v1.png")
SCRIPT_COPY <- file.path(OUTDIR, "Figure3_OFFICIAL_v1_script.R")

## ---------------- style ----------------
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

stage_levels <- c(
  "Early prenatal",
  "Mid prenatal",
  "Late prenatal",
  "Infancy",
  "Childhood",
  "Adolescence",
  "Adulthood"
)

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

required_files <- c(
  A = file.path(IN_DIR, "Figure3A_stage_affinity_profile.tsv"),
  B = file.path(IN_DIR, "Figure3_transport_plotdata_clean.tsv"),
  C = file.path(IN_DIR, "Figure3C_random_top_stage_frequency_clean.tsv"),
  C_obs = file.path(IN_DIR, "Figure3C_observed_top_stage_clean.tsv"),
  D = file.path(IN_DIR, "Figure3D_NTM2_random_transport_clean.tsv"),
  D_obs = file.path(IN_DIR, "Figure3D_NTM2_observed_transport_clean.tsv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing plotdata_v2 files:\n",
    paste(missing_files, collapse = "\n"),
    "\nPlease rerun figure3_native_modules_developmental_transport_v2.R first."
  )
}

## ---------------- read plot data ----------------
a_dt <- fread(required_files["A"])
b_dt <- fread(required_files["B"])
c_dt <- fread(required_files["C"])
c_obs <- fread(required_files["C_obs"])
d_dt <- fread(required_files["D"])
d_obs <- fread(required_files["D_obs"])

## enforce factors
a_dt[, module := factor(module, levels = module_levels)]
a_dt[, stage := factor(stage, levels = stage_levels)]

b_dt[, module := factor(module, levels = module_levels)]
b_dt[, stage := factor(stage, levels = stage_levels)]

c_dt[, module := factor(module, levels = c("NTM1\nASD-up", "NTM3\nASD-signed"))]
c_dt[, stage := factor(stage, levels = stage_levels)]

c_obs[, module := factor(module, levels = c("NTM1\nASD-up", "NTM3\nASD-signed"))]
c_obs[, stage := factor(stage, levels = stage_levels)]
c_obs[, legend_label := "Observed top stage"]

d_dt[, stage := factor(stage, levels = stage_levels)]
d_obs[, stage := factor(stage, levels = stage_levels)]
d_obs[, legend_label := "Observed NTM2"]

## ---------------- Panel A ----------------
panelA <- ggplot(a_dt, aes(x = stage, y = rel_value, group = module, colour = module)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.9) +
  facet_wrap(~ module, ncol = 1) +
  scale_colour_manual(values = module_cols, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(x = NULL, y = "Relative stage affinity", tag = "A") +
  theme_pub(7.8) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 7.6),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.8)
  )

## ---------------- Panel B ----------------
panelB <- ggplot(
  b_dt,
  aes(x = stage, y = factor(module, levels = rev(module_levels)), fill = rel_value)
) +
  geom_tile(colour = "white", linewidth = 0.45, width = 0.94, height = 0.78) +
  scale_fill_gradient(
    low = "#EEF2F6",
    high = COL_GREEN,
    limits = c(0, 1),
    name = "Relative\ntransport"
  ) +
  labs(x = NULL, y = NULL, tag = "B") +
  theme_pub(8.2) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.2),
    axis.text.y = element_text(size = 7.7),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7.2, colour = COL_TEXT),
    legend.text = element_text(size = 6.8, colour = COL_TEXT)
  )

## ---------------- Panel C ----------------
panelC <- ggplot(c_dt, aes(x = stage, y = frequency)) +
  geom_col(
    fill = "#C8D4DF",
    width = 0.72,
    colour = "white",
    linewidth = 0.18
  ) +
  geom_point(
    data = c_obs,
    aes(x = stage, y = y, shape = legend_label),
    inherit.aes = FALSE,
    fill = COL_GREEN,
    colour = COL_TEXT,
    size = 2.4,
    stroke = 0.25
  ) +
  facet_wrap(~ module, ncol = 1) +
  scale_shape_manual(values = c("Observed top stage" = 25)) +
  scale_y_continuous(
    limits = c(0, 1.10),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Random top-stage frequency", tag = "C") +
  guides(shape = guide_legend(override.aes = list(fill = COL_GREEN, colour = COL_TEXT))) +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.8),
    strip.text = element_text(size = 7.6),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0)
  )

## ---------------- Panel D ----------------
panelD <- ggplot(d_dt, aes(x = stage, y = value)) +
  geom_boxplot(
    width = 0.62,
    outlier.size = 0.35,
    outlier.alpha = 0.25,
    fill = "#DDE6EE",
    colour = COL_MUTED,
    linewidth = 0.32
  ) +
  geom_point(
    data = d_obs,
    aes(x = stage, y = observed_value, shape = legend_label),
    inherit.aes = FALSE,
    fill = COL_BLUE,
    colour = "white",
    size = 2.7,
    stroke = 0.25
  ) +
  scale_shape_manual(values = c("Observed NTM2" = 23)) +
  guides(shape = guide_legend(override.aes = list(fill = COL_BLUE, colour = "white"))) +
  labs(x = NULL, y = "Transport mass", tag = "D") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.8),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0)
  )

## ---------------- export ----------------
draw_combined <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 2,
    widths = unit(c(1.00, 1.18), "null"),
    heights = unit(c(1.00, 1.05), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 9.3, height = 7.1, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 2790, height = 2130, res = 300, type = "cairo")
draw_combined()
dev.off()

## copy script into frozen folder
this_script <- "step13_main_figures/scripts/figure3_native_modules_developmental_transport_FINAL_v1.R"
if (file.exists(this_script)) {
  file.copy(this_script, SCRIPT_COPY, overwrite = TRUE)
}

cat("Figure 3 OFFICIAL v1 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Script copy:\n")
cat(SCRIPT_COPY, "\n")
