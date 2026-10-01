source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
# S5 complete external models
s5forest <- function(file,title) {
  d <- read_panel(file) %>% mutate(row=short_module_top(program,top_n),row=factor(row,levels=rev(unique(row))))
  ggplot(d,aes(standardized_beta,row,colour=program,shape=factor(top_n))) + geom_vline(xintercept=0,colour="#777777",linewidth=.45) +
    geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.12,linewidth=.55) + geom_point(size=2.35,stroke=.65) + facet_wrap(~cohort,nrow=1) +
    scale_colour_manual(values=MODULE_COLORS,guide="none") + scale_shape_manual(values=c(`200`=21,`500`=24),name="Threshold",labels=c(`200`="Top 200",`500`="Top 500")) +
    labs(title=title,x="Standardized beta (95% CI)",y=NULL) + theme_pub() + theme(legend.position="bottom")
}
p5A <- s5forest("Supplementary_FigureS5A_NATIVE_SIGNED_MODULE.tsv","Native signed module")
p5B <- s5forest("Supplementary_FigureS5B_SIGNED_GENE_PPR.tsv","Signed gene PPR")
p5C <- s5forest("Supplementary_FigureS5C_UNSIGNED_GENE_PPR.tsv","Unsigned gene PPR")
d5D <- read_panel("Supplementary_FigureS5D_mapping_and_score_sd.tsv") %>%
  mutate(row=short_module_top(program, top_n),row=factor(row,levels=rev(unique(row))),
         method=recode(method,NATIVE_SIGNED_MODULE="Native",SIGNED_GENE_PPR="Signed PPR",UNSIGNED_GENE_PPR="Unsigned PPR"))
p5D1 <- ggplot(d5D,aes(n_common_genes,row,colour=method,shape=cohort)) + geom_point(size=2.2,stroke=.6,position=position_dodge(.45)) +
  scale_colour_manual(values=c(Native="#333333",`Signed PPR`="#777777",`Unsigned PPR`="#AAAAAA")) +
  labs(title="Mapped genes",x="Common genes",y=NULL,colour="Method",shape="Cohort") + theme_pub() +
  scale_x_continuous(breaks=c(500,5000)) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=1,byrow=TRUE)) +
  theme(legend.position="bottom",legend.box="vertical")
p5D2 <- ggplot(d5D,aes(score_SD_before_standardization,row,colour=method,shape=cohort)) + geom_point(size=2.2,stroke=.6,position=position_dodge(.45)) +
  scale_colour_manual(values=c(Native="#333333",`Signed PPR`="#777777",`Unsigned PPR`="#AAAAAA")) +
  labs(title="Score SD",x="Pre-standardization",y=NULL) + theme_pub() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position="none")
bottomS5 <- (tag_panel(p5C,"C") | tag_panel(p5D1|p5D2,"D")) + plot_layout(widths=c(.48,.52))
figS5 <- (tag_panel(p5A,"A") | tag_panel(p5B,"B")) / bottomS5
save_pub(figS5,"Supplementary_FigureS5",180,160,".")

