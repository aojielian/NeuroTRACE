suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(cowplot); library(grid)
  library(dplyr); library(tidyr); library(readr); library(scales); library(ggrepel)
})

BASE <- Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd())
OUT <- Sys.getenv("NEUROTRACE_FIGURE_OUTPUT", unset=file.path(BASE,"generated_figures"))
PD <- file.path(BASE,"figure_source_data")
FONT <- "Helvetica"

MODULE_COLORS <- c(NTM1_ASD_up="#7B3294", NTM2_ASD_down="#2166AC", NTM3_ASD_signed="#1B7837")
MODULE_LABELS <- c(NTM1_ASD_up="NTM1 upregulated", NTM2_ASD_down="NTM2 downregulated", NTM3_ASD_signed="NTM3 mixed")
STAGE_ORDER <- c("early_prenatal","mid_prenatal","late_prenatal","early_postnatal","childhood_adolescence","adulthood")
STAGE_LABELS <- c(early_prenatal="EPr",mid_prenatal="MPr",late_prenatal="LPr",early_postnatal="EPost",childhood_adolescence="C/A",adulthood="Adult")
STAGE_COLORS <- c(early_prenatal="#8C6BB1", mid_prenatal="#6C83C5", late_prenatal="#2B5DA8",
                  early_postnatal="#E6A532", childhood_adolescence="#E87532", adulthood="#C73E3A")
METHOD_LABELS <- c(
  DUAL_HEAD_GENE_FIRST_PPR="Dual-head PPR", DIRECT_NATIVE_PROJECTION="Native projection",
  SIGNED_GENE_FIRST_PPR_SINGLE_HEAD="Signed single-head", UNSIGNED_GENE_PPR="Unsigned PPR",
  GENE_ONLY_PPR_STAGE3="Gene-only PPR", HETERO_PPR_STANDARD="Heterogeneous PPR",
  GRAPH_OT="Graph OT", SIMPLE_OT="Simple OT", MEAN_SIGNATURE="Mean signature",
  STAGE_CORRELATION="Stage correlation", RANDOM="Random", LEGACY_NEUROTRACE_SIM_PROXY="Simulation proxy"
)
METHOD_COLORS <- setNames(rep("#424242", length(METHOD_LABELS)), names(METHOD_LABELS))
METHOD_SHAPES <- setNames(c(21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6), names(METHOD_LABELS))

read_panel <- function(name) readr::read_tsv(file.path(PD, name), show_col_types=FALSE, progress=FALSE)

theme_pub <- function(base_size=7.0) {
  theme_classic(base_family=FONT, base_size=base_size) +
    theme(
      text=element_text(family=FONT, colour="#202020"),
      axis.title=element_text(size=7.2), axis.text=element_text(size=7.0, colour="#202020"),
      axis.line=element_line(linewidth=.45, colour="#303030"), axis.ticks=element_line(linewidth=.4),
      strip.background=element_rect(fill="#F2F2F2", colour=NA), strip.text=element_text(size=7.2, face="bold"),
      legend.title=element_text(size=7.0, face="bold"),
      legend.title.position="top", legend.key.width=unit(5.5,"mm"), legend.text=element_text(size=7.0),
      plot.title=element_text(size=7.6, face="bold", hjust=0, lineheight=.96,
                              margin=margin(0,0,1.5,0)),
      plot.margin=margin(2.2, 2.2, 2.2, 2.2, unit="mm"),
      panel.spacing=unit(1.8, "mm"), legend.key.height=unit(3.2, "mm"),
       legend.spacing.x=unit(1.2, "mm"),
      legend.box.spacing=unit(1.2, "mm")
    )
}

theme_heat <- function() theme_pub() + theme(axis.line=element_blank(), axis.ticks=element_blank())

tag_panel <- function(p, tag) {
  cowplot::ggdraw() + cowplot::draw_plot(p, x=.052, y=.024, width=.925, height=.946) +
    cowplot::draw_label(tag, x=.004, y=.997, hjust=0, vjust=1,
                       fontfamily=FONT, fontface="bold", size=11)
}

short_module_top <- function(program, top_n) {
  paste0(MODULE_LABELS[program], "\nTop ", top_n)
}

mm_to_in <- function(x) x / 25.4

save_pub <- function(plot, stem, width_mm, height_mm, subdir) {
  dir <- file.path(OUT, subdir); dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  pdf_path <- file.path(dir, paste0(stem, ".pdf"))
  png_path <- file.path(dir, paste0(stem, ".png"))
  tif_path <- file.path(dir, paste0(stem, ".tiff"))
  grDevices::cairo_pdf(pdf_path, width=mm_to_in(width_mm), height=mm_to_in(height_mm), family=FONT,
                       onefile=FALSE, bg="white")
  print(plot); grDevices::dev.off()
  ragg::agg_png(png_path, width=width_mm, height=height_mm, units="mm", res=600, background="white")
  print(plot); grDevices::dev.off()
  ragg::agg_tiff(tif_path, width=width_mm, height=height_mm, units="mm", res=600,
                 compression="lzw", background="white")
  print(plot); grDevices::dev.off()
  message(stem, " -> PDF/PNG/TIFF")
}

module_factor <- function(x) factor(x, levels=names(MODULE_LABELS), labels=unname(MODULE_LABELS))
stage_factor <- function(x) factor(x, levels=STAGE_ORDER, labels=unname(STAGE_LABELS))
clean_condition <- function(x) {
  z <- gsub("_", " ", x)
  z <- sub("combined stress", "combined stress", z)
  tools::toTitleCase(z)
}

forest_base <- function(d, title=NULL) {
  d <- d %>% mutate(module_top=paste0(MODULE_LABELS[program], " · Top ", top_n),
                    module_top=factor(module_top, levels=rev(unique(module_top))))
  ggplot(d, aes(x=standardized_beta, y=module_top, colour=program, shape=factor(top_n))) +
    geom_vline(xintercept=0, linewidth=.45, colour="#777777") +
    geom_errorbarh(aes(xmin=ci_low, xmax=ci_high), height=.12, linewidth=.55,
                   position=position_dodge(width=.35)) +
    geom_point(size=2.3, stroke=.55, position=position_dodge(width=.35)) +
    facet_wrap(~cohort, nrow=1) + scale_colour_manual(values=MODULE_COLORS, guide="none") +
    scale_shape_manual(values=c(`200`=21, `500`=24), name="Threshold") +
    labs(x="Standardized beta (95% CI)", y=NULL, title=title) + theme_pub()
}
