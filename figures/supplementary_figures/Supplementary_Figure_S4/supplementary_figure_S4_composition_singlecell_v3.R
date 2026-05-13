suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Supplementary_FigureS4")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

COMPOSITION_MODEL_FILE <- file.path(ROOT, "step10_composition_adjustment/results/08_step10E_NTM_models_before_after_composition.tsv")
NNLS_MODEL_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/06_step11H_Gandal_NNLS_model_comparison_summary.tsv")

PE_LOCALIZATION_FILE <- file.path(ROOT, "step10_single_cell_validation/results/PsychENCODE_12_step10F_sc_NTM_celltype_localization_summary.tsv")
VEL_LOCALIZATION_FILE <- file.path(ROOT, "step10_single_cell_validation/results/Velmeshev_12_step10F_sc_NTM_celltype_localization_summary.tsv")

PDF_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S4_v3.pdf")
PNG_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S4_v3.png")
DIAG_OUT <- file.path(OUTDIR, "Supplementary_Figure_S4_v3_diagnostics.txt")
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

clean_celltype <- function(x) {
  y <- as.character(x)
  y <- gsub("_", " ", y)

  y[grepl("^EXN$|^EX$", y, ignore.case = TRUE)] <- "Excitatory neurons"
  y[grepl("^INN$|^IN$", y, ignore.case = TRUE)] <- "Inhibitory neurons"
  y[grepl("^AST$|astro", y, ignore.case = TRUE)] <- "Astrocytes"
  y[grepl("^OPC$|OPC", y, ignore.case = TRUE)] <- "OPCs"
  y[grepl("^ODC$|oligo", y, ignore.case = TRUE)] <- "Oligodendrocytes"
  y[grepl("^MG$|micro", y, ignore.case = TRUE)] <- "Microglia"
  y[grepl("^END$|endo", y, ignore.case = TRUE)] <- "Endothelial"

  y
}

celltype_levels <- c(
  "Excitatory neurons",
  "Inhibitory neurons",
  "Astrocytes",
  "OPCs",
  "Oligodendrocytes",
  "Microglia",
  "Endothelial"
)

## ============================================================
## Panel A: marker-adjustment attenuation line plot
## ============================================================

comp <- fread(COMPOSITION_MODEL_FILE)

a_raw <- comp[grepl("_top500$", ntm_score)]
a_raw[, program := gsub("_top500$", "", ntm_score)]
a_raw[, module := pretty_module(program)]
a_raw <- a_raw[module %in% module_levels]

base_dt <- a_raw[model_type == "M0_base_covariates", .(
  cohort,
  module,
  base_beta = beta
)]

adj_dt <- a_raw[model_type %in% c("M1_marker_PC1_PC2", "M2_all_marker_scores"), .(
  cohort,
  module,
  model_type,
  beta
)]

adj_dt[, adjustment := fifelse(
  model_type == "M1_marker_PC1_PC2",
  "Marker PCs",
  "All marker scores"
)]

a_dt <- merge(adj_dt, base_dt, by = c("cohort", "module"), all.x = TRUE)
a_dt <- a_dt[!is.na(base_beta) & base_beta != 0 & (beta / base_beta) > 0]
a_dt[, ratio := beta / base_beta]
a_dt[, log2_ratio := log2(ratio)]

a_plot <- a_dt[
  ,
  .(
    log2_ratio = median(log2_ratio, na.rm = TRUE),
    log2_min = min(log2_ratio, na.rm = TRUE),
    log2_max = max(log2_ratio, na.rm = TRUE),
    ratio = median(ratio, na.rm = TRUE),
    n_cohorts = uniqueN(cohort)
  ),
  by = .(module, adjustment)
]

a_plot[, module := factor(module, levels = module_levels)]
a_plot[, adjustment := factor(adjustment, levels = c("Marker PCs", "All marker scores"))]
a_plot[, seg_ymin := pmin(0, log2_ratio)]
a_plot[, seg_ymax := pmax(0, log2_ratio)]

fwrite(a_plot, file.path(PLOTDIR, "S4A_marker_adjustment_log2_ratio.tsv"), sep = "\t")

adj_cols <- c(
  "Marker PCs" = COL_GREEN,
  "All marker scores" = COL_PURPLE
)

pd_a <- position_dodge(width = 0.48)

panelA <- ggplot(a_plot, aes(x = module, y = log2_ratio, colour = adjustment, shape = adjustment)) +
  geom_hline(yintercept = 0, colour = COL_LIGHT, linewidth = 0.45, linetype = "dashed") +
  geom_linerange(
    aes(ymin = seg_ymin, ymax = seg_ymax),
    position = pd_a,
    linewidth = 0.65,
    alpha = 0.9
  ) +
  geom_errorbar(
    aes(ymin = log2_min, ymax = log2_max),
    position = pd_a,
    width = 0.13,
    linewidth = 0.34,
    alpha = 0.72
  ) +
  geom_point(position = pd_a, size = 2.8, stroke = 0.25) +
  scale_colour_manual(values = adj_cols, name = NULL) +
  scale_shape_manual(values = c("Marker PCs" = 16, "All marker scores" = 17), name = NULL) +
  labs(x = NULL, y = expression(log[2]~"(adjusted/base beta)"), tag = "A") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(size = 7.2),
    legend.position = "bottom"
  )


## ============================================================
## Panel B: NNLS attenuation ratio line plot
## ============================================================

nnls <- fread(NNLS_MODEL_FILE)
b_dt <- nnls[top_n == 500]
b_dt[, module := pretty_module(program)]
b_dt <- b_dt[module %in% module_levels]

b_plot <- melt(
  b_dt,
  id.vars = c("dataset", "program", "module", "top_n"),
  measure.vars = c("PC12_beta_ratio", "props_beta_ratio"),
  variable.name = "adjustment",
  value.name = "beta_ratio"
)

b_plot[, adjustment := fifelse(
  adjustment == "PC12_beta_ratio",
  "NNLS-PC",
  "NNLS-proportion"
)]

b_plot[, module := factor(module, levels = module_levels)]
b_plot[, adjustment := factor(adjustment, levels = c("NNLS-PC", "NNLS-proportion"))]
b_plot[, seg_ymin := pmin(1, beta_ratio)]
b_plot[, seg_ymax := pmax(1, beta_ratio)]

fwrite(b_plot, file.path(PLOTDIR, "S4B_NNLS_beta_ratio.tsv"), sep = "\t")

nnls_cols <- c(
  "NNLS-PC" = COL_GREEN,
  "NNLS-proportion" = COL_BLUE
)

pd_b <- position_dodge(width = 0.48)

panelB <- ggplot(b_plot, aes(x = module, y = beta_ratio, colour = adjustment, shape = adjustment)) +
  geom_hline(yintercept = 1, colour = COL_LIGHT, linewidth = 0.45, linetype = "dashed") +
  geom_linerange(
    aes(ymin = seg_ymin, ymax = seg_ymax),
    position = pd_b,
    linewidth = 0.65,
    alpha = 0.9
  ) +
  geom_point(position = pd_b, size = 2.8, stroke = 0.25) +
  scale_colour_manual(values = nnls_cols, name = NULL) +
  scale_shape_manual(values = c("NNLS-PC" = 16, "NNLS-proportion" = 17), name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Adjusted/base beta ratio", tag = "B") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(size = 7.2),
    legend.position = "bottom"
  )


## ============================================================
## Shared single-nucleus localization data
## ============================================================

standardize_loc <- function(path, dataset_label) {
  dt <- fread(path)

  out <- data.table(
    dataset = dataset_label,
    celltype = clean_celltype(dt$broad_celltype),
    module = pretty_module(dt$ntm_score_name),
    value = as.numeric(dt$mean_score),
    n_donors = as.numeric(dt$n_donors)
  )

  out <- out[module %in% module_levels]
  out[, module := factor(module, levels = module_levels)]
  out[, celltype := factor(celltype, levels = rev(celltype_levels))]
  out
}

pe_sc <- standardize_loc(PE_LOCALIZATION_FILE, "PsychENCODE")
vel_sc <- standardize_loc(VEL_LOCALIZATION_FILE, "Velmeshev")

sc_all <- rbindlist(list(pe_sc, vel_sc), fill = TRUE)
lim_sc <- max(abs(sc_all$value), na.rm = TRUE)
if (!is.finite(lim_sc) || lim_sc == 0) lim_sc <- 1

fwrite(pe_sc, file.path(PLOTDIR, "S4C_PsychENCODE_single_nucleus_localization.tsv"), sep = "\t")
fwrite(vel_sc, file.path(PLOTDIR, "S4D_Velmeshev_single_nucleus_localization.tsv"), sep = "\t")

## ============================================================
## Panel C: PsychENCODE localization
## ============================================================

panelC <- ggplot(pe_sc, aes(x = module, y = celltype, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.40, width = 0.90, height = 0.82) +
  scale_fill_gradient2(
    low = COL_BLUE,
    mid = "#F4F6F8",
    high = COL_RED,
    midpoint = 0,
    limits = c(-lim_sc, lim_sc),
    name = "Mean\nscore"
  ) +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_pub(7.6) +
  theme(
    axis.text.x = element_text(size = 7.0),
    axis.text.y = element_text(size = 7.0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

## ============================================================
## Panel D: Velmeshev localization
## ============================================================

panelD <- ggplot(vel_sc, aes(x = module, y = celltype, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.40, width = 0.90, height = 0.82) +
  scale_fill_gradient2(
    low = COL_BLUE,
    mid = "#F4F6F8",
    high = COL_RED,
    midpoint = 0,
    limits = c(-lim_sc, lim_sc),
    name = "Mean\nscore"
  ) +
  labs(x = NULL, y = NULL, tag = "D") +
  theme_pub(7.6) +
  theme(
    axis.text.x = element_text(size = 7.0),
    axis.text.y = element_text(size = 7.0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

## ============================================================
## Diagnostics
## ============================================================

diag <- c(
  "Supplementary Figure S4 v3 diagnostics",
  paste("Generated:", as.character(Sys.time())),
  "",
  paste("COMPOSITION_MODEL_FILE:", COMPOSITION_MODEL_FILE),
  paste("NNLS_MODEL_FILE:", NNLS_MODEL_FILE),
  paste("PE_LOCALIZATION_FILE:", PE_LOCALIZATION_FILE),
  paste("VEL_LOCALIZATION_FILE:", VEL_LOCALIZATION_FILE),
  "",
  paste("a_plot rows:", nrow(a_plot)),
  paste("b_plot rows:", nrow(b_plot)),
  paste("pe_sc rows:", nrow(pe_sc)),
  paste("vel_sc rows:", nrow(vel_sc)),
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
    widths = unit(c(1.02, 1.02), "null"),
    heights = unit(c(0.88, 1.12), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 10.4, height = 7.8, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 3120, height = 2340, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Supplementary Figure S4 v3 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plotdata:\n")
cat(PLOTDIR, "\n")
cat("Diagnostics:\n")
cat(DIAG_OUT, "\n")
