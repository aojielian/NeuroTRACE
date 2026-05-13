suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

ROOT <- "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project"
OUTDIR <- file.path(ROOT, "paper_figures", "Supplementary_FigureS2")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

TRANSPORT_FILE <- file.path(ROOT, "step05_optimal_transport_alignment/results/10_step05B_NTM_stage_graph_transport.tsv")
RANDOM_FREQ_FILE <- file.path(ROOT, "step10_transport_null/results/105_step10B_random_top_stage_frequency.tsv")
NTM2_RANDOM_FILE <- file.path(ROOT, "step11_robustness_sensitivity/results/02_step11D_NTM2_random_transport.tsv.gz")
STEP12C_DECISION_FILE <- file.path(ROOT, "step12_revision_strengthening/results/step12C_stage_alignment_weight_sensitivity_v3/05_step12C_v3_decision_table.tsv")

PDF_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S2_v8.pdf")
PNG_OUT  <- file.path(OUTDIR, "Supplementary_Figure_S2_v8.png")
DIAG_OUT <- file.path(OUTDIR, "Supplementary_Figure_S2_v8_diagnostics.txt")
PLOTDIR  <- file.path(OUTDIR, "plotdata_v8")
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
    warning("Missing file: ", path)
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
  y <- gsub("module:", "", y)
  y <- gsub("\\|.*$", "", y)
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

stage_levels <- c(
  "Early prenatal",
  "Mid prenatal",
  "Late prenatal",
  "Childhood",
  "Adolescence",
  "Adulthood"
)

stage_clean <- function(x) {
  y0 <- as.character(x)
  y <- tolower(y0)
  y <- gsub("^stage_", "", y)
  y <- gsub("^transport_", "", y)
  y <- gsub("^prob_", "", y)
  y <- gsub("\\.", "_", y)
  y <- gsub("-", "_", y)
  y <- gsub("\\s+", "_", y)

  out <- y0
  out[grepl("early.*prenatal|early_prenatal|^early$|^e$", y)] <- "Early prenatal"
  out[grepl("mid.*prenatal|mid_prenatal|^mid$|^m$", y)] <- "Mid prenatal"
  out[grepl("late.*prenatal|late_prenatal|^late$|^l$", y)] <- "Late prenatal"
  out[grepl("child", y)] <- "Childhood"
  out[grepl("adolesc|ado", y)] <- "Adolescence"
  out[grepl("adult|adulthood|^ad$", y)] <- "Adulthood"

  out
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
      stage = stage_levels[seq_len(min(length(ord_levels), length(stage_levels)))]
    )
    fallback <- map$stage[match(ord, map$stage_order)]
    out[is.na(out)] <- fallback[is.na(out)]
  }

  out
}

make_placeholder <- function(tag, label) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.55,
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

extract_stage_long <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  module_col <- pick_col(dt, c("^module$", "^ntm$", "module_family", "module_node", "^program$", "signature"))
  topn_col <- pick_col(dt, c("^top_n$", "^topn$", "^k$", "cutoff", "module_top_n"))
  stage_col <- pick_col(dt, c("^stage$", "stage_node", "top_stage", "transport_stage", "brainspan_stage", "developmental_stage"))
  value_col <- pick_col(
    dt,
    c("transport_probability", "transport.*mass", "transport.*prob", "^probability$", "^mass$", "affinity", "score", "^value$"),
    numeric_only = TRUE,
    exclude = c(topn_col)
  )

  if (!is.na(module_col) && !is.na(stage_col) && !is.na(value_col)) {
    out <- data.table(
      module_raw = dt[[module_col]],
      stage_raw = dt[[stage_col]],
      value = as.numeric(dt[[value_col]]),
      top_n = if (!is.na(topn_col)) suppressWarnings(as.numeric(dt[[topn_col]])) else NA_real_
    )
    out <- out[!is.na(value)]
    out[, module := pretty_module(module_raw)]
    out[, stage := stage_clean(stage_raw)]
    out <- out[module %in% module_levels & stage %in% stage_levels]
    return(out)
  }

  return(NULL)
}

## ============================================================
## Read inputs
## ============================================================

transport <- read_dt(TRANSPORT_FILE)
random_freq <- read_dt(RANDOM_FREQ_FILE)
ntm2_random <- read_dt(NTM2_RANDOM_FILE)
step12c <- read_dt(STEP12C_DECISION_FILE)

diag <- c(
  "Supplementary Figure S2 v8 diagnostics",
  paste("Generated:", as.character(Sys.time())),
  "",
  paste("TRANSPORT_FILE:", TRANSPORT_FILE),
  paste("RANDOM_FREQ_FILE:", RANDOM_FREQ_FILE),
  paste("NTM2_RANDOM_FILE:", NTM2_RANDOM_FILE),
  paste("STEP12C_DECISION_FILE:", STEP12C_DECISION_FILE),
  "",
  paste("transport rows:", if (is.null(transport)) "NA" else nrow(transport)),
  paste("random_freq rows:", if (is.null(random_freq)) "NA" else nrow(random_freq)),
  paste("ntm2_random rows:", if (is.null(ntm2_random)) "NA" else nrow(ntm2_random)),
  paste("step12c rows:", if (is.null(step12c)) "NA" else nrow(step12c)),
  "",
  "transport columns:",
  if (is.null(transport)) "NULL" else paste(names(transport), collapse = ", "),
  "",
  "random_freq columns:",
  if (is.null(random_freq)) "NULL" else paste(names(random_freq), collapse = ", "),
  "",
  "ntm2_random columns:",
  if (is.null(ntm2_random)) "NULL" else paste(names(ntm2_random), collapse = ", "),
  "",
  "step12c columns:",
  if (is.null(step12c)) "NULL" else paste(names(step12c), collapse = ", ")
)

## ============================================================
## Panel A: top200/top500 transport sensitivity
## ============================================================

panelA <- NULL

if (!is.null(transport)) {
  t_long <- extract_stage_long(transport)

  if (!is.null(t_long) && nrow(t_long) > 0) {
    t_long <- t_long[top_n %in% c(200, 500)]
    if (nrow(t_long) == 0) {
      top_keep <- sort(unique(transport$top_n))[sort(unique(transport$top_n)) %in% c(100, 200, 500, 1000)]
      t_long <- extract_stage_long(transport)[top_n %in% tail(top_keep, 2)]
    }

    t_plot <- t_long[
      ,
      .(value = mean(value, na.rm = TRUE)),
      by = .(module, top_n, stage)
    ]

    t_plot[, rel_value := value / max(value, na.rm = TRUE), by = .(module, top_n)]
    t_plot[, module := factor(module, levels = rev(module_levels))]
    t_plot[, stage := factor(stage, levels = stage_levels)]
    t_plot[, cutoff := factor(paste0("top", top_n), levels = paste0("top", sort(unique(t_plot$top_n))))]

    fwrite(t_plot, file.path(PLOTDIR, "S2A_transport_topn_sensitivity.tsv"), sep = "\t")

    panelA <- ggplot(t_plot, aes(x = stage, y = module, fill = rel_value)) +
      geom_tile(colour = "white", linewidth = 0.42, width = 0.94, height = 0.78) +
      facet_wrap(~ cutoff, nrow = 1) +
      scale_fill_gradient(low = "#EEF2F6", high = COL_GREEN, limits = c(0, 1), name = "Relative\ntransport") +
      labs(x = NULL, y = NULL, tag = "A") +
      theme_pub(7.8) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.8),
        axis.text.y = element_text(size = 7.3),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.8, "lines")
      )
  }
}

if (is.null(panelA)) {
  panelA <- make_placeholder("A", "Top-n transport sensitivity\ncould not be generated.\nSee diagnostics.")
}

## Observed top stage for B/C overlays
obs_top <- NULL
if (exists("t_plot") && nrow(t_plot) > 0) {
  obs_top <- t_plot[
    ,
    .SD[which.max(value)],
    by = .(module, top_n)
  ][
    ,
    .(module, top_n, stage, observed_value = value, observed_rel_value = rel_value)
  ]
}

## ============================================================
## Panel B: matched-random top-stage frequency for NTM1/NTM3
## ============================================================

panelB <- NULL

if (!is.null(random_freq) && nrow(random_freq) > 0) {
  module_col <- pick_col(random_freq, c("^module$", "^ntm$", "module_family", "^program$", "signature"))
  topn_col <- pick_col(random_freq, c("^top_n$", "^topn$", "^k$", "cutoff", "module_top_n"))
  stage_col <- pick_col(random_freq, c("top_stage", "^stage$", "transport_stage", "brainspan_stage", "developmental_stage"))
  freq_col <- pick_col(random_freq, c("^frequency$", "freq", "proportion", "fraction", "prob", "rate"), numeric_only = TRUE)
  count_col <- pick_col(random_freq, c("^n$", "count"), numeric_only = TRUE)

  if (!is.na(module_col) && !is.na(stage_col)) {
    b_dt <- copy(random_freq)
    b_dt[, module := pretty_module(get(module_col))]
    b_dt <- b_dt[module %in% c("NTM1\nASD-up", "NTM3\nASD-signed")]
    b_dt[, stage := stage_clean(get(stage_col))]
    b_dt <- b_dt[stage %in% stage_levels]

    if (!is.na(topn_col)) {
      b_dt[, top_n := suppressWarnings(as.numeric(get(topn_col)))]
      b_dt <- b_dt[top_n %in% c(200, 500)]
    } else {
      b_dt[, top_n := 500]
    }

    if (!is.na(freq_col)) {
      b_dt[, frequency := as.numeric(get(freq_col))]
    } else if (!is.na(count_col)) {
      b_dt[, frequency := as.numeric(get(count_col))]
      b_dt[, frequency := frequency / sum(frequency, na.rm = TRUE), by = .(module, top_n)]
    } else {
      b_dt <- b_dt[, .N, by = .(module, top_n, stage)]
      b_dt[, frequency := N / sum(N, na.rm = TRUE), by = .(module, top_n)]
    }

    b_dt <- b_dt[
      ,
      .(frequency = sum(frequency, na.rm = TRUE)),
      by = .(module, top_n, stage)
    ]

    b_dt[, module := factor(module, levels = c("NTM1\nASD-up", "NTM3\nASD-signed"))]
    b_dt[, stage := factor(stage, levels = stage_levels)]
    b_dt[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]

    obs_b <- NULL
    if (!is.null(obs_top)) {
      obs_b <- obs_top[module %in% c("NTM1\nASD-up", "NTM3\nASD-signed") & top_n %in% c(200, 500)]
      obs_b[, module := factor(module, levels = c("NTM1\nASD-up", "NTM3\nASD-signed"))]
      obs_b[, stage := factor(stage, levels = stage_levels)]
      obs_b[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]
      obs_b[, y := 1.04]
    }

    fwrite(b_dt, file.path(PLOTDIR, "S2B_random_top_stage_frequency.tsv"), sep = "\t")
    if (!is.null(obs_b)) fwrite(obs_b, file.path(PLOTDIR, "S2B_observed_top_stage.tsv"), sep = "\t")

    panelB <- ggplot(b_dt, aes(x = stage, y = frequency)) +
      geom_col(fill = "#C8D4DF", width = 0.72, colour = "white", linewidth = 0.18) +
      facet_grid(module ~ cutoff) +
      scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Random top-stage frequency", tag = "B") +
      theme_pub(7.4) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.4),
        strip.text = element_text(size = 7.0),
        legend.position = "none",
        panel.spacing = unit(0.45, "lines")
      )

    if (!is.null(obs_b) && nrow(obs_b) > 0) {
      panelB <- panelB +
        geom_point(
          data = obs_b,
          aes(x = stage, y = y),
          inherit.aes = FALSE,
          shape = 25,
          fill = COL_GREEN,
          colour = COL_TEXT,
          size = 2.0,
          stroke = 0.25
        )
    }
  }
}

if (is.null(panelB)) {
  panelB <- make_placeholder("B", "Matched-random top-stage frequency\ncould not be generated.\nSee diagnostics.")
}

## ============================================================
## Panel C: Stage-alignment stability fraction
## ============================================================

panelD <- NULL

if (!is.null(step12c) && nrow(step12c) > 0) {
  d_dt <- copy(step12c)

  d_dt[, module := pretty_module(module_family)]
  d_dt <- d_dt[module %in% module_levels]
  d_dt <- d_dt[top_n %in% c(200, 500)]

  ## Use the expected stability fraction for each module:
  ## NTM1/NTM3: late-prenatal top fraction;
  ## NTM2: postnatal top fraction.
  d_dt[, stability_fraction := fifelse(
    module == "NTM2\nASD-down",
    postnatal_top_fraction,
    late_prenatal_top_fraction
  )]

  d_dt[, stability_axis := fifelse(
    module == "NTM2\nASD-down",
    "Postnatal/adult-like top-stage fraction",
    "Late-prenatal top-stage fraction"
  )]

  d_dt[, module := factor(module, levels = rev(module_levels))]
  d_dt[, cutoff := factor(paste0("top", top_n), levels = c("top200", "top500"))]

  fwrite(d_dt, file.path(PLOTDIR, "S2D_stage_alignment_stability_fraction.tsv"), sep = "\t")

  panelD <- ggplot(d_dt, aes(x = cutoff, y = module, fill = stability_fraction)) +
    geom_tile(colour = "white", linewidth = 0.45, width = 0.88, height = 0.78) +
    geom_text(
      aes(label = sprintf("%.2f", stability_fraction)),
      colour = COL_TEXT,
      size = 2.7,
      family = font_family
    ) +
    scale_fill_gradient(
      low = "#EEF2F6",
      high = COL_GREEN,
      limits = c(0, 1),
      name = "Stability\nfraction"
    ) +
    labs(x = NULL, y = NULL, tag = "C") +
    theme_pub(7.8) +
    theme(
      axis.text.x = element_text(size = 7.2),
      axis.text.y = element_text(size = 7.3),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

if (is.null(panelD)) {
  panelD <- make_placeholder("D", "Stage-alignment stability fraction\ncould not be generated.\nSee diagnostics.")
}


## ============================================================
## Diagnostics
## ============================================================

diag <- c(
  diag,
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
    widths = unit(c(1.05, 1.00), "null"),
    heights = unit(c(1.00, 1.00), "null")
  )))

  ## A: full-width top panel
  print(panelA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))

  ## B and C: bottom row
  print(panelB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(panelD, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  popViewport()
}

pdf(PDF_OUT, width = 10.4, height = 7.4, useDingbats = FALSE)
draw_combined()
dev.off()

png(PNG_OUT, width = 3120, height = 2220, res = 300, type = "cairo")
draw_combined()
dev.off()

cat("Supplementary Figure S2 v8 written to:\n")
cat(PDF_OUT, "\n")
cat(PNG_OUT, "\n")
cat("Plotdata:\n")
cat(PLOTDIR, "\n")
cat("Diagnostics:\n")
cat(DIAG_OUT, "\n")

