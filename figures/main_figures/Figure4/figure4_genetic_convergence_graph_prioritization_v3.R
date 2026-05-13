suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Figure4")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

FULL_TRANSPORT_FILE <- file.path(ROOT, "step05_optimal_transport_alignment/results/10_step05B_NTM_stage_graph_transport.tsv")
NO_RISK_TRANSPORT_FILE <- file.path(ROOT, "step10_risk_free_ablation/results/91_step10A_no_risk_NTM_stage_graph_transport.tsv")
SFARI_STABILITY_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/04b_step11A_no_risk_SFARI_stability_summary.corrected.tsv")
PRIORITY_SFARI_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/08_step11B_gene_priority_SFARI_summary.tsv")
GNOMAD_SUMMARY_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12D_gnomad_constraint/04_step12D1_constraint_validation_summary.tsv")

PDF_OUT <- file.path(OUTDIR, "Figure4_genetic_convergence_graph_prioritization_v3.pdf")
PNG_OUT <- file.path(OUTDIR, "Figure4_genetic_convergence_graph_prioritization_v3.png")
PLOTDATA_DIR <- file.path(OUTDIR, "plotdata_v3")
dir.create(PLOTDATA_DIR, recursive = TRUE, showWarnings = FALSE)

COL_TEXT   <- "#263238"
COL_MUTED  <- "#6B7280"
COL_LIGHT  <- "#D7DEE6"
COL_GREEN  <- "#4F9A73"
COL_BLUE   <- "#4777A8"
COL_ORANGE <- "#D96C4A"
COL_RED    <- "#B64A4A"

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
      plot.margin = margin(6, 7, 6, 7)
    )
}

module_levels <- c("NTM1\nASD-up", "NTM2\nASD-down", "NTM3\nASD-signed")

stage_levels_6 <- c(
  "Early prenatal",
  "Mid prenatal",
  "Late prenatal",
  "Childhood",
  "Adolescence",
  "Adulthood"
)

pretty_module <- function(x) {
  y <- as.character(x)
  y <- gsub("module:", "", y)
  y <- gsub("\\|.*$", "", y)
  out <- y
  out[grepl("NTM1", y, ignore.case = TRUE)] <- "NTM1\nASD-up"
  out[grepl("NTM2", y, ignore.case = TRUE)] <- "NTM2\nASD-down"
  out[grepl("NTM3", y, ignore.case = TRUE)] <- "NTM3\nASD-signed"
  out
}

clean_risk_short <- function(x) {
  y <- as.character(x)
  y <- gsub("^risk_set:", "", y)
  y <- gsub("_fixed$", "", y)
  y <- gsub("^is_", "", y)

  out <- y
  out[grepl("SFARI_S_plus_1_2", y)] <- "S+1+2"
  out[grepl("SFARI_S_plus_1$", y)] <- "S+1"
  out[grepl("SFARI_all", y)] <- "All"
  out[grepl("SFARI_score_le1", y)] <- "Score <=1"
  out[grepl("SFARI_score_le2", y)] <- "Score <=2"
  out[grepl("SFARI_score_le3", y)] <- "Score <=3"
  out[grepl("SFARI_syndromic", y)] <- "Syndromic"
  out
}

clean_source_short <- function(x) {
  y <- as.character(x)
  y[y == "native_module_genes"] <- "Native"
  y[y == "native_module"] <- "Native"
  y[y == "full_graph_priority"] <- "Full graph"
  y[y == "graph_priority"] <- "Graph priority"
  y[y == "no_risk_graph_priority"] <- "No-risk graph"
  y
}

stage_label_from_transport <- function(dt) {
  out <- rep(NA_character_, nrow(dt))

  raw <- if ("stage_node" %in% names(dt)) as.character(dt$stage_node) else rep("", nrow(dt))
  raw_low <- tolower(raw)

  out[grepl("early.*prenatal|early_prenatal", raw_low)] <- "Early prenatal"
  out[grepl("mid.*prenatal|mid_prenatal", raw_low)] <- "Mid prenatal"
  out[grepl("late.*prenatal|late_prenatal", raw_low)] <- "Late prenatal"
  out[grepl("child", raw_low)] <- "Childhood"
  out[grepl("adolesc|ado", raw_low)] <- "Adolescence"
  out[grepl("adult", raw_low)] <- "Adulthood"

  if ("stage_order" %in% names(dt)) {
    ord <- suppressWarnings(as.integer(dt$stage_order))
    ord_levels <- sort(unique(ord[is.finite(ord)]))
    map <- data.table(
      stage_order = ord_levels,
      stage = stage_levels_6[seq_len(min(length(ord_levels), length(stage_levels_6)))]
    )
    fallback <- map$stage[match(ord, map$stage_order)]
    out[is.na(out)] <- fallback[is.na(out)]
  }

  out
}

make_transport_dt <- function(path, graph_label) {
  dt <- fread(path)
  stopifnot(all(c("module_node", "top_n", "stage_order", "transport_probability") %in% names(dt)))

  dt[, module := pretty_module(module_node)]
  dt[, stage := stage_label_from_transport(.SD)]
  dt <- dt[module %in% module_levels & stage %in% stage_levels_6]

  top_use <- if (500 %in% dt$top_n) 500 else max(dt$top_n, na.rm = TRUE)
  dt <- dt[top_n == top_use]

  out <- dt[
    ,
    .(transport_probability = mean(transport_probability, na.rm = TRUE)),
    by = .(module, stage)
  ]
  out[, graph := graph_label]
  out[, rel_transport := transport_probability / max(transport_probability, na.rm = TRUE), by = module]
  out[, module := factor(module, levels = module_levels)]
  out[, stage := factor(stage, levels = stage_levels_6)]
  out
}

## ============================================================
## A. Full graph vs no-risk transport
## ============================================================

full_t <- make_transport_dt(FULL_TRANSPORT_FILE, "Full graph")
norisk_t <- make_transport_dt(NO_RISK_TRANSPORT_FILE, "No-risk graph")
a_dt <- rbindlist(list(full_t, norisk_t), fill = TRUE)

fwrite(a_dt, file.path(PLOTDATA_DIR, "Figure4A_full_vs_norisk_transport.tsv"), sep = "\t")

panelA <- ggplot(
  a_dt,
  aes(x = stage, y = factor(module, levels = rev(module_levels)), fill = rel_transport)
) +
  geom_tile(colour = "white", linewidth = 0.42, width = 0.94, height = 0.78) +
  facet_wrap(~ graph, nrow = 1) +
  scale_fill_gradient(low = "#EEF2F6", high = COL_GREEN, limits = c(0, 1), name = "Relative\ntransport") +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.7),
    axis.text.y = element_text(size = 7.2),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    panel.spacing = unit(0.7, "lines")
  )

## ============================================================
## B. No-risk NTM2 SFARI stability
## ============================================================

sfari <- fread(SFARI_STABILITY_FILE)
b_dt <- sfari[module_family == "NTM2_ASD_down" & module_top_n %in% c(200, 500)]

b_dt[, risk_label := clean_risk_short(risk_set)]

risk_order <- c("Syndromic", "S+1", "S+1+2", "Score <=1", "Score <=2", "Score <=3", "All")
risk_order <- risk_order[risk_order %in% unique(b_dt$risk_label)]

b_dt[, risk_label := factor(risk_label, levels = rev(risk_order))]
b_dt[, cutoff := factor(paste0("top", module_top_n), levels = c("top200", "top500"))]

fwrite(b_dt, file.path(PLOTDATA_DIR, "Figure4B_no_risk_SFARI_stability.tsv"), sep = "\t")

panelB <- ggplot(b_dt, aes(x = median_or, y = risk_label, colour = cutoff)) +
  geom_vline(xintercept = 1, colour = COL_LIGHT, linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(group = risk_label), colour = COL_LIGHT, linewidth = 0.35) +
  geom_point(aes(size = fdr_lt_0_05_fraction), alpha = 0.95) +
  scale_colour_manual(values = c("top200" = COL_BLUE, "top500" = COL_GREEN), name = "Gene cutoff") +
  scale_size_continuous(range = c(1.8, 4.2), limits = c(0, 1), name = "FDR<0.05\nfraction") +
  labs(x = "Median SFARI enrichment odds ratio", y = NULL, tag = "B") +
  theme_pub(7.8) +
  theme(
    axis.text.y = element_text(size = 7.2),
    legend.position = "right"
  )

## ============================================================
## C. Held-out gnomAD constraint convergence, top500 only
## ============================================================

gnomad <- fread(GNOMAD_SUMMARY_FILE)

c_dt <- gnomad[top_n == 500]

## For graph-priority rows, use a fixed priority_cutoff=500 for comparability.
## Native-module rows often do not need the same priority_cutoff logic.
if ("priority_cutoff" %in% names(c_dt)) {
  c_dt <- c_dt[
    (source_type %in% c("graph_priority", "no_risk_graph_priority") & priority_cutoff == 500) |
      source_type %in% c("native_module", "native_module_genes")
  ]
}

c_dt[, module := pretty_module(module_family)]
c_dt <- c_dt[module %in% module_levels]

c_dt[, source_label := clean_source_short(source_type)]
source_levels_c <- c("Graph priority", "No-risk graph", "Native")
source_levels_c <- source_levels_c[source_levels_c %in% unique(c_dt$source_label)]

c_dt[, source_label := factor(source_label, levels = source_levels_c)]
c_dt[, module := factor(module, levels = rev(module_levels))]
c_dt[, label_n := paste0(n_fdr_lt_0_05, "/", n_constraint_sets)]

fwrite(c_dt, file.path(PLOTDATA_DIR, "Figure4C_gnomAD_constraint_convergence_top500.tsv"), sep = "\t")

panelC <- ggplot(c_dt, aes(x = source_label, y = module, fill = median_or)) +
  geom_tile(colour = "white", linewidth = 0.42, width = 0.88, height = 0.78) +
  geom_text(aes(label = label_n), size = 2.35, colour = COL_TEXT, family = font_family) +
  scale_fill_gradient2(
    low = COL_BLUE,
    mid = "#F4F6F8",
    high = COL_RED,
    midpoint = 1,
    name = "Median\nOR"
  ) +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_pub(7.8) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 7.2),
    axis.text.y = element_text(size = 7.3),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

## ============================================================
## D. Graph-prioritized vs native SFARI enrichment
## ============================================================

priority <- fread(PRIORITY_SFARI_FILE)
d_dt <- priority[
  module_family == "NTM2_ASD_down" &
    module_top_n %in% c(200, 500) &
    gene_set_source %in% c("full_graph_priority", "no_risk_graph_priority", "native_module_genes")
]

d_dt[, source_label := clean_source_short(gene_set_source)]
d_dt[, source_label := factor(source_label, levels = c("Native", "Full graph", "No-risk graph"))]
d_dt[, cutoff := factor(paste0("top", module_top_n), levels = c("top200", "top500"))]
d_dt[, label_n := paste0(n_fdr_global_lt_0_05, "/", n_risk_sets)]

fwrite(d_dt, file.path(PLOTDATA_DIR, "Figure4D_graph_priority_vs_native_SFARI.tsv"), sep = "\t")

source_cols <- c(
  "Native" = COL_BLUE,
  "Full graph" = COL_GREEN,
  "No-risk graph" = COL_ORANGE
)

panelD <- ggplot(d_dt, aes(x = source_label, y = median_or, fill = source_label)) +
  geom_hline(yintercept = 1, colour = COL_LIGHT, linewidth = 0.4, linetype = "dashed") +
  geom_col(width = 0.65, colour = "white", linewidth = 0.22) +
  geom_text(aes(label = label_n), vjust = -0.55, size = 2.4, colour = COL_TEXT, family = font_family) +
  facet_wrap(~ cutoff, nrow = 1) +
  scale_fill_manual(values = source_cols, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Median SFARI enrichment odds ratio", tag = "D") +
  theme_pub(8.0) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 7.2),
    legend.position = "none",
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
    widths = unit(c(1.22, 1.03), "null"),
    heights = unit(c(1.00, 1.05), "null")
  )))

  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panelB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(panelC, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 9.5, height = 7.1, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 2850, height = 2130, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Figure 4 v3 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plot data written to:\n")
cat(PLOTDATA_DIR, "\n")
