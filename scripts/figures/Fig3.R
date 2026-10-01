source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
dA <- read_panel("Fig3A_NTM3_channel_profiles.tsv") %>%
  mutate(stage=stage_factor(state_id),threshold=paste0("Top ",top_n),channel=recode(channel,positive="Increased channel\n(q_pos)",negative="Decreased channel\n(q_neg)"))
pA <- ggplot(dA,aes(stage,channel_cosine,fill=channel)) + geom_hline(yintercept=0,colour="#777777",linewidth=.35) +
  geom_col(width=.72,na.rm=TRUE) + facet_grid(threshold~channel) + scale_fill_manual(values=c(`Increased channel\n(q_pos)`="#7B3294",`Decreased channel\n(q_neg)`="#2166AC"),guide="none") +
  labs(title="NTM3 mixed: channel profiles",x=NULL,y="Cosine localization") + theme_pub() +
  theme(axis.text.x=element_text(angle=34,hjust=1,size=7),panel.spacing=unit(1.4,"mm"))

dB <- read_panel("Fig3B_all_channel_heatmap.tsv") %>%
  mutate(row=paste0(MODULE_LABELS[program]," · Top ",top_n," · ",recode(channel,positive="Increased",negative="Decreased")),
         row=factor(row,levels=rev(unique(row))),stage=stage_factor(state_id))
pB <- ggplot(dB,aes(stage,row,fill=channel_cosine)) + geom_tile(colour="white",linewidth=.35) +
  scale_fill_viridis_c(option="C",na.value="#F1F1F1",name="Cosine") +
  labs(title="All module/channel profiles",x=NULL,y=NULL) + theme_heat() +
  theme(axis.text.x=element_text(angle=34,hjust=1,size=7),axis.text.y=element_text(size=6.4),legend.position="bottom",
        plot.title=element_text(size=7.4))

dC <- read_panel("Fig3C_readout_contrasts.tsv") %>%
  mutate(row=short_module_top(program, top_n),row=factor(row,levels=rev(unique(row))),
         readout=recode(readout,q_pos="Increased channel (q_pos)",q_neg="Decreased channel (q_neg)",q_magnitude="Combined magnitude",q_signed="Signed orientation"),
         readout=factor(readout,levels=c("Increased channel (q_pos)","Decreased channel (q_neg)","Combined magnitude","Signed orientation")))
pC <- ggplot(dC,aes(localization_contrast,row,colour=readout,shape=readout)) +
  geom_vline(xintercept=0,colour="#777777",linewidth=.45) + geom_point(size=2.5,stroke=.65,position=position_dodge(width=.42),na.rm=TRUE) +
  scale_colour_manual(values=c("#7B3294","#2166AC","#7F7F7F","#202020"),name="Readout") + scale_shape_manual(values=c(21,22,24,4),name="Readout") +
  labs(title="Localization contrast\nby readout",x="Prenatal − postnatal contrast",y=NULL) + theme_pub() +
  guides(shape=guide_legend(ncol=1),colour=guide_legend(ncol=1)) +
  theme(legend.position="bottom",legend.box="vertical",plot.title=element_text(size=7.4))

dD <- read_panel("Fig3D_channel_composition.tsv") %>%
  mutate(row=short_module_top(program, top_n),row=factor(row,levels=rev(unique(row))),
         channel=recode(channel,positive="Increased",negative="Decreased"))
pD1 <- ggplot(dD,aes(channel_mass_share_of_total_abs_weight,row,fill=channel)) +
  geom_col(position="stack",width=.68) + scale_fill_manual(values=c(Increased="#7B3294",Decreased="#2166AC"),name=NULL) +
  scale_x_continuous(labels=percent_format(),breaks=c(0,.5,1)) + labs(title="Weight share",x="Share",y=NULL) + theme_pub() +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + theme(legend.position="bottom",plot.title=element_text(size=7.1))
pD2 <- ggplot(dD,aes(n_restart_genes,row,colour=channel,shape=channel)) +
  geom_line(aes(group=row),colour="#B8B8B8",linewidth=.45) + geom_point(size=2.3,stroke=.65) +
  scale_colour_manual(values=c(Increased="#7B3294",Decreased="#2166AC"),guide="none") + scale_shape_manual(values=c(21,22),guide="none") +
  scale_x_continuous(breaks=c(0,100,200)) + labs(title="Restart genes",x="Genes",y=NULL) + theme_pub() +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=7.1))
pD <- (pD1 | pD2) + plot_layout(widths=c(.57,.43),guides="collect") & theme(legend.position="bottom")

fig <- (tag_panel(pA,"A") | tag_panel(pB,"B")) / (tag_panel(pC,"C") | tag_panel(pD,"D")) +
  plot_layout(heights=c(.49,.51))
save_pub(fig,"Fig3",180,155,".")
