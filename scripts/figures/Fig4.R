source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
dA <- read_panel("Fig4A_forest.tsv") %>% mutate(method=recode(method,NATIVE_SIGNED_MODULE="Native",SIGNED_GENE_PPR="Signed PPR"),
  row=short_module_top(program, top_n),row=factor(row,levels=rev(unique(row))))
pA <- ggplot(dA,aes(standardized_beta,row,colour=program,shape=method)) + geom_vline(xintercept=0,colour="#777777",linewidth=.45) +
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.12,linewidth=.55,position=position_dodge(width=.42)) +
  geom_point(size=2.4,stroke=.65,position=position_dodge(width=.42)) + facet_wrap(~cohort,nrow=1) +
  scale_colour_manual(values=MODULE_COLORS,guide="none") + scale_shape_manual(values=c(Native=21,`Signed PPR`=24),name="Score") +
  labs(title="Standardized effects",x="Standardized beta (95% CI)",y=NULL) + theme_pub() +
  guides(shape=guide_legend(nrow=1,byrow=TRUE)) + theme(legend.position="bottom",plot.title=element_text(size=7.4))

dB <- read_panel("Fig4B_paired_beta.tsv") %>% mutate(module=module_factor(program),threshold=factor(top_n))
pB <- ggplot(dB,aes(NATIVE_SIGNED_MODULE_standardized_beta,SIGNED_GENE_PPR_standardized_beta,colour=program,shape=threshold)) +
  geom_abline(slope=1,intercept=0,linetype=2,linewidth=.45,colour="#777777") + geom_point(size=2.7,stroke=.65) +
  facet_wrap(~cohort,nrow=1) + scale_colour_manual(values=MODULE_COLORS,name="Module",breaks=names(MODULE_LABELS),labels=unname(MODULE_LABELS)) +
  scale_shape_manual(values=c(`200`=21,`500`=24),labels=c(`200`="Top 200",`500`="Top 500"),name="Threshold") + coord_equal() +
  labs(title="Native vs signed-PPR effects",x="Native standardized beta",y="Signed-PPR standardized beta") + theme_pub() +
  guides(colour=guide_legend(ncol=1,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE)) +
  theme(legend.position="bottom",legend.box="vertical",plot.title=element_text(size=7.4))

dC <- read_panel("Fig4C_effect_FDR_matrix.tsv") %>%
  mutate(row=paste0(MODULE_LABELS[program]," · Top ",top_n),row=factor(row,levels=rev(unique(row))),
         method=factor(method,levels=c("native signed module","signed gene PPR"),labels=c("Native signed module","Signed gene PPR")),sig=fdr<.05)
lim <- max(abs(dC$standardized_beta),na.rm=TRUE)
pC <- ggplot(dC,aes(method,row,fill=standardized_beta)) + geom_tile(colour="white",linewidth=.5) +
  facet_grid(cohort~.,scales="free_y",space="free_y") +
  geom_point(aes(shape=sig),size=2.1,colour="black",fill="black",stroke=.55) +
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-lim,lim),name="Std. beta") +
  scale_shape_manual(values=c(`TRUE`=21,`FALSE`=1),name="FDR < 0.05",labels=c(`FALSE`="No",`TRUE`="Yes")) +
  labs(title="Effect/FDR matrix",x=NULL,y=NULL) + theme_heat() +
  theme(axis.text.x=element_text(angle=20,hjust=1),axis.text.y=element_text(size=6.4),legend.position="bottom",legend.box="vertical") +
  guides(fill=guide_colorbar(barwidth=unit(28,"mm")),shape=guide_legend(nrow=1))

dD <- read_panel("Fig4D_cohort_summary.tsv")
bars <- dD %>% select(cohort,direction_concordant_n,native_fdr_significant_n,signed_ppr_fdr_significant_n) %>%
  pivot_longer(-cohort,names_to="metric",values_to="count") %>%
  mutate(metric=recode(metric,direction_concordant_n="Concordant",native_fdr_significant_n="Native FDR < 0.05",signed_ppr_fdr_significant_n="Signed PPR FDR < 0.05"))
pD1 <- ggplot(bars,aes(cohort,count,fill=metric)) + geom_col(position=position_dodge(.72),width=.65) +
  geom_text(aes(label=paste0(count,"/6")),position=position_dodge(.72),vjust=-.25,size=2.5,family=FONT) +
  scale_fill_manual(values=c("#6B6B6B","#B2182B","#D98E8E"),name=NULL) + scale_y_continuous(limits=c(0,6.8),breaks=0:6) +
  labs(title="Concordance and FDR",x=NULL,y="Comparisons") + theme_pub() +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) + theme(legend.position="bottom",plot.title=element_text(size=7.2))
pD2 <- ggplot(dD,aes(cohort,median_signed_minus_native_beta,group=1)) + geom_hline(yintercept=0,colour="#777777",linewidth=.45) +
  geom_line(colour="#8B1A1A",linewidth=.6) + geom_point(shape=21,fill="#8B1A1A",size=2.8) +
  labs(title="Signed-PPR effect change\nrelative to native",x=NULL,y="Median difference") + theme_pub() +
  theme(plot.title=element_text(size=7.1))
pD <- pD1 / pD2 + plot_layout(heights=c(.62,.38))

fig <- (tag_panel(pA,"A") | tag_panel(pB,"B")) / (tag_panel(pC,"C") | tag_panel(pD,"D")) +
  plot_layout(heights=c(.49,.51))
save_pub(fig,"Fig4",180,158,".")
