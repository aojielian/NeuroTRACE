source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
# S2 condition-level simulation heatmaps
method_order_s2 <- names(METHOD_LABELS)
METHOD_LABELS["LEGACY_NEUROTRACE_SIM_PROXY"] <- "Simulation proxy"
metric_heat <- function(file, value_col, title) {
  d <- read_panel(file) %>% mutate(method=factor(METHOD_LABELS[method],levels=unname(METHOD_LABELS[method_order_s2])),condition=clean_condition(condition),
                                  benchmark=recode(benchmark,step02A="High signal",step02B="Challenging"))
  ggplot(d,aes(method,condition,fill=.data[[value_col]])) + geom_tile(colour="white",linewidth=.3) +
    facet_grid(benchmark~.,scales="free_y",space="free_y") + scale_fill_gradient(low="#EFF3FF",high="#08519C",limits=c(0,1),na.value="#ECECEC",breaks=c(0,.5,1),name="Mean (higher = darker)") +
    labs(title=title,x=NULL,y=NULL) + theme_heat() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7),axis.text.y=element_text(size=7),legend.position="bottom",
          panel.spacing=unit(1.2,"mm"),plot.title=element_text(size=7.4))
}
p2A <- metric_heat("Supplementary_FigureS2A_alignment_hit_top1_mean.tsv","alignment_hit_top1_mean","Top-stage accuracy")
p2B <- metric_heat("Supplementary_FigureS2B_true_stage_mass_mean.tsv","true_stage_mass_mean","True-stage mass")
p2C <- metric_heat("Supplementary_FigureS2C_gene_auprc_mean.tsv","gene_auprc_mean","Gene AUPRC")
p2D <- metric_heat("Supplementary_FigureS2D_precision_at_50_mean.tsv","precision_at_50_mean","Precision@50")
figS2 <- (tag_panel(p2A,"A") | tag_panel(p2B,"B")) / (tag_panel(p2C,"C") | tag_panel(p2D,"D"))
save_pub(figS2,"Supplementary_FigureS2",180,175,".")


