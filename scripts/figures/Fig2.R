source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/null_display.R"))
dA <- read_panel("Fig2A_stage_profiles.tsv") %>%
  mutate(module=module_factor(program),
         threshold=paste0("Top ",top_n), stage=stage_factor(state_id))
ymax <- max(dA$magnitude_cosine,na.rm=TRUE)*1.08
pA <- ggplot(dA,aes(x=stage,y=magnitude_cosine,fill=state_id)) +
  geom_hline(yintercept=0,linewidth=.35,colour="#777777") + geom_col(width=.74) +
  facet_grid(module~threshold,labeller=labeller(module=function(x) sub(" ","\n",x,fixed=TRUE))) + scale_fill_manual(values=STAGE_COLORS,name="Developmental stage",
                                                   breaks=STAGE_ORDER, labels=STAGE_LABELS) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  coord_cartesian(ylim=c(min(0,min(dA$magnitude_cosine)),ymax)) +
  labs(title="Six-stage magnitude-head profiles",x=NULL,y="Cosine localization") + theme_pub() +
  theme(axis.text.x=element_text(angle=28,hjust=1,size=7),legend.position="bottom",
        legend.direction="horizontal",panel.spacing=unit(1.4,"mm"),
        plot.title=element_text(size=7.6),strip.text.y=element_text(angle=0,size=6.5))

dB <- read_panel("Fig2B_contrast_null_interval.tsv") %>%
  mutate(module_top=short_module_top(program, top_n),
         module_top=factor(module_top,levels=rev(unique(module_top))))
pB <- ggplot(dB,aes(y=module_top)) + geom_vline(xintercept=0,colour="#777777",linewidth=.4) +
  geom_errorbarh(aes(xmin=null_localization_contrast_q025,xmax=null_localization_contrast_q975),height=0,linewidth=.8,colour="#7A7A7A") +
  geom_point(aes(x=null_localization_contrast_median),shape=0,size=2.5,stroke=.7,colour="#4D4D4D") +
  geom_point(aes(x=observed_localization_contrast,fill=program),shape=21,size=2.7,stroke=.55,colour="black") +
  scale_fill_manual(values=MODULE_COLORS,guide="none") +
  labs(title="Observed vs matched-null\ncontrasts",x="Prenatal − postnatal contrast",y=NULL) + theme_pub() +
  theme(plot.title=element_text(size=7.4))


legendB <- ggplot() +
 annotate("point",x=.05,y=3,shape=21,size=2.3,fill=MODULE_COLORS[1],colour="black") +
 annotate("point",x=.05,y=2,shape=0,size=2.3,colour="#4D4D4D") +
 annotate("segment",x=.02,xend=.08,y=1,yend=1,linewidth=.7,colour="#7A7A7A") +
 annotate("text",x=.12,y=c(3,2,1),label=c("Observed contrast (module color)","Null median","Null 2.5–97.5 percentiles"),hjust=0,family=FONT,size=2.1) +
 coord_cartesian(xlim=c(0,1),ylim=c(.6,3.4),clip="off") + theme_void() + theme(plot.margin=margin(0,2,0,2,"mm"))
pB <- pB / legendB + plot_layout(heights=c(.85,.15))

dC <- read_panel("Fig2C_empirical_null_replicates.tsv")
pC <- null_display(dC, bins=24, density=TRUE, supplementary=FALSE) +
 labs(title="Empirical matched-null distributions",x="Null contrast",y="Density") +
 theme(plot.title=element_text(size=7.2),axis.text=element_text(size=6.4),
       legend.text=element_text(size=5.8),strip.text.y=element_text(size=6.2,angle=0),
       panel.spacing=unit(1.2,"mm"))

dD <- read_panel("Fig2D_observed_stage_frequency.tsv") %>%
  mutate(module_top=paste0(MODULE_LABELS[program]," · Top ",top_n),
         module_top=factor(module_top,levels=rev(unique(module_top))),
         stage_lab=STAGE_LABELS[observed_top_stage])
pD <- ggplot(dD,aes(x=null_top_stage_frequency_observed,y=module_top,colour=program)) +
  geom_segment(aes(x=0,xend=null_top_stage_frequency_observed,yend=module_top),colour="#D0D0D0",linewidth=.5) +
  geom_point(size=2.4) + geom_text(aes(label=stage_lab),hjust=ifelse(dD$null_top_stage_frequency_observed>.72,1.16,-.16),
                                  family=FONT,size=2.5,colour="#333333") +
  scale_colour_manual(values=MODULE_COLORS,guide="none") +
  scale_x_continuous(limits=c(0,1),labels=percent_format(accuracy=1),expand=expansion(mult=c(.02,.04))) +
  labs(title="Null frequency of\nobserved top stage",x="Frequency",y=NULL) + theme_pub() +
  theme(axis.text.y=element_text(size=6.4),plot.title=element_text(size=7.2))

bottom <- tag_panel(pB,"B") | (tag_panel(pC,"C") / tag_panel(pD,"D") + plot_layout(heights=c(.58,.42)))
fig <- tag_panel(pA,"A") / bottom + plot_layout(heights=c(.43,.57),widths=1)
fig <- fig & plot_layout()
save_pub(fig,"Fig2",180,180,".")
