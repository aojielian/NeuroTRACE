source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
# S7 descriptive cross-disorder context; frozen single-model MDD intervals.
d7A <- read_panel("Supplementary_FigureS7A_cross_disorder_heatmap.tsv") %>%
  mutate(module=sub(" ASD.*","",module),module=factor(module,levels=c("NTM1","NTM2","NTM3")),disease=factor(disease,levels=c("MDD","BD","SCZ")))
lim7 <- max(abs(d7A$beta))
p7A <- ggplot(d7A,aes(module,disease,fill=beta)) + geom_tile(colour="white",linewidth=.65) +
  geom_text(aes(label=sprintf("%.3f",beta)),family=FONT,size=2.8,colour="#202020") +
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-lim7,lim7),name="Mean standardized beta") +
  labs(title="Cross-disorder native-NTM effects",subtitle="Top 500 coefficient means",x=NULL,y=NULL) + theme_heat() +
  theme(axis.text.x=element_text(angle=0,hjust=.5),legend.position="bottom",plot.subtitle=element_text(size=7))
d7B <- read_panel("Supplementary_FigureS7B_disease_summary.tsv") %>% mutate(disease=factor(disease,levels=c("SCZ","BD","MDD")))
p7B <- ggplot(d7B,aes(disease,positive_rate,fill=disease)) + geom_col(width=.58) +
  geom_text(aes(label=label),vjust=-.45,size=3,family=FONT) +
  scale_fill_manual(values=c(SCZ="#555555",BD="#555555",MDD="#555555"),guide="none") +
  scale_y_continuous(limits=c(0,1.15),breaks=c(0,.5,1),labels=percent_format()) +
  labs(title="Direction concordance",x=NULL,y="Concordant proportion") + theme_pub()
d7C <- read_panel("Supplementary_FigureS7C_MDD_forest.tsv") %>%
  mutate(module=sub(" ASD.*","",module),module=factor(module,levels=rev(c("NTM1","NTM2","NTM3"))),
         label=sprintf("FDR %.3f",within_disease_fdr))
p7C <- ggplot(d7C,aes(beta_target_vs_Control,module,colour=program)) +
  geom_vline(xintercept=0,colour="#777777",linewidth=.45) +
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.14,linewidth=.6) + geom_point(size=2.8) +
  geom_text(aes(x=.30,label=label),hjust=0,size=2.5,family=FONT,colour="#333333") +
  scale_colour_manual(values=MODULE_COLORS,guide="none") + scale_x_continuous(limits=c(-.33,.52),breaks=c(-.3,0,.3)) +
  labs(title="MDD: GSE53987, Top 500",subtitle="Within-disease FDR; single-model intervals",x="Standardized beta (95% CI)",y=NULL) + theme_pub() +
  theme(plot.subtitle=element_text(size=7))
figS7 <- tag_panel(p7A,"A") | (tag_panel(p7B,"B") / tag_panel(p7C,"C") + plot_layout(heights=c(.48,.52))) +
  plot_layout(widths=c(.47,.53))
save_pub(figS7,"Supplementary_FigureS7",180,155,".")
