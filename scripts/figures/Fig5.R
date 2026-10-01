source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
METHOD_LABELS["LEGACY_NEUROTRACE_SIM_PROXY"] <- "Simulation proxy"
QUANT_LOW <- "#EFF3FF"
QUANT_HIGH <- "#08519C"
ACCENT <- "#D97706"
METHOD_COLORS <- setNames(rep("#666666", length(METHOD_LABELS)), names(METHOD_LABELS))
METHOD_COLORS["DUAL_HEAD_GENE_FIRST_PPR"] <- ACCENT
METHOD_FILLS <- setNames(rep(NA_character_, length(METHOD_LABELS)), names(METHOD_LABELS))
METHOD_FILLS["DUAL_HEAD_GENE_FIRST_PPR"] <- ACCENT
quant_scale <- function(title="Mean (higher = darker)") scale_fill_gradient(
  low=QUANT_LOW, high=QUANT_HIGH, limits=c(0,1), na.value="#ECECEC", name=title,
  breaks=c(0,.5,1), guide=guide_colorbar(title.position="top",barwidth=unit(24,"mm"),barheight=unit(2,"mm"),order=1))

dA <- read_panel("Fig5A_condition_matrix.tsv") %>%
  mutate(condition=clean_condition(condition),
         design=recode(design,`challenging-signal design`="Challenging",`high-signal design`="High signal"),
         perturbation_dimension=recode(perturbation_dimension,
           adversarial_composition="Adversarial", batch_sd="Batch SD",
           cohort_noise_sd="Cohort noise", composition_sd="Composition SD",
           dropout_rate="Dropout", false_prior_rate="False prior",
           signal_strength="Signal", state_overlap="State overlap",
           true_prior_recall="Prior recall"))
pA <- ggplot(dA,aes(x=perturbation_dimension,y=condition,fill=within_design_scaled_value)) +
  geom_tile(colour="white",linewidth=.45) + facet_wrap(~design,ncol=1,scales="free_y") +
  quant_scale("Relative\nintensity") + guides(fill=guide_colorbar(direction="vertical",title.position="top",barwidth=unit(2.5,"mm"),barheight=unit(18,"mm"))) +
  labs(title="Simulation condition matrix",x=NULL,y=NULL) + theme_heat() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7),legend.position="right",
        panel.spacing=unit(1.2,"mm"),plot.title=element_text(size=7.5))

dG <- read_panel("Fig5B_developmental_metrics.tsv") %>%
  mutate(method_label=METHOD_LABELS[method], benchmark=recode(benchmark,step02A="High signal",step02B="Challenging"),
         method_label=factor(method_label,levels=rev(unname(METHOD_LABELS))))
dev <- dG %>% select(benchmark,method,method_label,alignment_hit_top1_mean,true_stage_mass_mean) %>%
  pivot_longer(c(alignment_hit_top1_mean,true_stage_mass_mean),names_to="metric",values_to="value") %>%
  mutate(metric=recode(metric,alignment_hit_top1_mean="Top-stage accuracy",true_stage_mass_mean="True-stage mass"))
pB <- ggplot(dev,aes(x=value,y=method_label,fill=value,shape=benchmark)) +
  geom_segment(aes(x=0,xend=value,yend=method_label),colour="#D6D6D6",linewidth=.35,na.rm=TRUE) +
  geom_point(size=2.0,stroke=.65,colour="#505050",na.rm=TRUE) + facet_wrap(~metric,nrow=1,scales="free_x") +
  quant_scale() + scale_shape_manual(values=c(`High signal`=21,Challenging=24),name="Design",guide=guide_legend(order=2,override.aes=list(fill="white",colour="#505050"))) +
  labs(title="Developmental metrics",x="Mean",y=NULL) + theme_pub() + theme(panel.spacing.x=unit(4,"mm")) + theme(legend.position="bottom",legend.box="vertical",plot.title=element_text(size=7.5))

gene <- dG %>% select(benchmark,method,method_label,gene_auprc_mean,precision_at_50_mean) %>%
  pivot_longer(c(gene_auprc_mean,precision_at_50_mean),names_to="metric",values_to="value") %>%
  mutate(metric=recode(metric,gene_auprc_mean="Gene AUPRC",precision_at_50_mean="Precision@50"))
pC <- ggplot(gene,aes(x=value,y=method_label,fill=value,shape=benchmark)) +
  geom_segment(aes(x=0,xend=value,yend=method_label),colour="#D6D6D6",linewidth=.35) +
  geom_point(size=2.0,stroke=.65,colour="#505050") + facet_wrap(~metric,nrow=1,scales="free_x") +
  quant_scale() + scale_shape_manual(values=c(`High signal`=21,Challenging=24),name="Design",guide=guide_legend(order=2,override.aes=list(fill="white",colour="#505050"))) +
  labs(title="Gene-prioritization metrics",x="Mean",y=NULL) + theme_pub() + theme(panel.spacing.x=unit(4,"mm")) + theme(legend.position="bottom",legend.box="vertical",plot.title=element_text(size=7.5))

dD <- read_panel("Fig5D_tradeoff.tsv") %>% mutate(
  label=recode(method,DUAL_HEAD_GENE_FIRST_PPR="Dual-head",DIRECT_NATIVE_PROJECTION="Native",
    SIGNED_GENE_FIRST_PPR_SINGLE_HEAD="Signed 1-head",UNSIGNED_GENE_PPR="Unsigned",
    GENE_ONLY_PPR_STAGE3="Gene-only",HETERO_PPR_STANDARD="Heterog.",GRAPH_OT="Graph OT",
    SIMPLE_OT="Simple OT",MEAN_SIGNATURE="Mean sig.",STAGE_CORRELATION="Stage corr.",
    RANDOM="Random",LEGACY_NEUROTRACE_SIM_PROXY="Sim. proxy"),
  label_x=recode(method,DUAL_HEAD_GENE_FIRST_PPR=.58,DIRECT_NATIVE_PROJECTION=.275,
    SIGNED_GENE_FIRST_PPR_SINGLE_HEAD=.31,UNSIGNED_GENE_PPR=.58,GENE_ONLY_PPR_STAGE3=.56,
    HETERO_PPR_STANDARD=.48,GRAPH_OT=.49,SIMPLE_OT=.275,MEAN_SIGNATURE=.275,
    STAGE_CORRELATION=.30,RANDOM=.20,LEGACY_NEUROTRACE_SIM_PROXY=.63),
  label_y=recode(method,DUAL_HEAD_GENE_FIRST_PPR=.252,DIRECT_NATIVE_PROJECTION=.300,
    SIGNED_GENE_FIRST_PPR_SINGLE_HEAD=.250,UNSIGNED_GENE_PPR=.190,GENE_ONLY_PPR_STAGE3=.365,
    HETERO_PPR_STANDARD=.285,GRAPH_OT=.405,SIMPLE_OT=.348,MEAN_SIGNATURE=.378,
    STAGE_CORRELATION=.178,RANDOM=.078,LEGACY_NEUROTRACE_SIM_PROXY=.292))
# Annotation anchors only revised at plot specification; all data-point coordinates remain frozen.
ref <- dD %>% filter(!method %in% c("RANDOM","DUAL_HEAD_GENE_FIRST_PPR")) %>%
  summarise(x=median(alignment_hit_top1_mean,na.rm=TRUE),y=median(gene_auprc_mean,na.rm=TRUE))
pD <- ggplot(dD,aes(x=alignment_hit_top1_mean,y=gene_auprc_mean,colour=method,shape=method,size=precision_at_50_mean)) +
  geom_vline(xintercept=ref$x,linetype=2,linewidth=.45,colour="#777777") +
  geom_hline(yintercept=ref$y,linetype=2,linewidth=.45,colour="#777777") +
  geom_point(aes(fill=method),stroke=.7) +
  geom_segment(data=dD,
    aes(x=alignment_hit_top1_mean,y=gene_auprc_mean,xend=label_x,yend=label_y),
    inherit.aes=FALSE,linewidth=.25,colour="#999999",show.legend=FALSE) +
  geom_label(aes(x=label_x,y=label_y,label=label),family=FONT,size=2.5,fill="white",linewidth=0,label.padding=unit(.06,"lines"),label.r=unit(0,"lines"),show.legend=FALSE) +
  scale_colour_manual(values=METHOD_COLORS,guide="none") + scale_fill_manual(values=METHOD_FILLS,guide="none",na.value=NA) + scale_shape_manual(values=METHOD_SHAPES,guide="none") +
  scale_size_continuous(range=c(2.0,4.1),name="Precision@50") +
  labs(title="Challenging-design trade-off",x="Top-stage accuracy",y="Gene AUPRC") + theme_pub() + theme(panel.spacing.x=unit(4,"mm")) +
  coord_cartesian(clip="off") + theme(legend.position="bottom",plot.title=element_text(size=7.5),plot.margin=margin(3,5,3,3,"mm"))

fig <- (tag_panel(pA,"A") | tag_panel(pB,"B")) / (tag_panel(pC,"C") | tag_panel(pD,"D")) +
  plot_layout(widths=c(1,1),heights=c(.55,.45))
save_pub(fig,"Fig5",180,180,".")
