# Display-only preparation. Histogram maxima set annotation spacing; no inference is calculated.
null_display <- function(d,bins,density=FALSE,supplementary=FALSE) {
 labels <- c(NTM1_ASD_up="NTM1 upregulated",NTM2_ASD_down="NTM2 downregulated",NTM3_ASD_signed="NTM3 mixed")
 d <- d %>% mutate(module=factor(program,levels=names(labels),labels=unname(labels)),threshold=paste0("Top ",top_n))
 ln <- d %>% distinct(program,module,top_n,threshold,observed_localization_contrast,
   null_localization_contrast_median,null_localization_contrast_q025,null_localization_contrast_q975,
   empirical_P_two_sided,BH_FDR_two_sided) %>% arrange(module,top_n) %>%
   mutate(panel_letter=LETTERS[1:6],ann=if(supplementary) sprintf("obs %.3f\nP %.3g; FDR %.3g",observed_localization_contrast,empirical_P_two_sided,BH_FDR_two_sided) else sprintf("P %.3g; FDR %.3g",empirical_P_two_sided,BH_FDR_two_sided))
 h <- ggplot(d,aes(localization_contrast)) +
   geom_histogram(aes(y=if(.env$density) after_stat(density) else after_stat(count)),bins=bins,fill=if(density) "#CED7E0" else "#CBD5DF",colour="white",linewidth=if(density) .18 else .2) +
   geom_vline(data=ln,aes(xintercept=observed_localization_contrast),alpha=0) +
   facet_grid(module~threshold,scales="free_y",labeller=labeller(module=function(x) if(supplementary) x else sub(" ","\n",x,fixed=TRUE)))
 built <- ggplot_build(h)
 maxima <- built$data[[1]] %>% group_by(PANEL) %>% summarise(hist_top=max(y),.groups="drop")
 maxima <- left_join(built$layout$layout %>% select(PANEL,module,threshold),maxima,by="PANEL") %>%
   group_by(module) %>% mutate(hist_top=max(hist_top)) %>% ungroup()
 ln <- left_join(ln,maxima %>% select(module,threshold,hist_top),by=c("module","threshold")) %>%
   mutate(line_top=hist_top*1.06,ann_y=hist_top*if(supplementary) 1.27 else 1.22,canvas_top=hist_top*if(supplementary) 1.46 else 1.40)
 p <- ggplot(d,aes(localization_contrast))
 if(supplementary) p <- p + geom_rect(data=ln,aes(xmin=null_localization_contrast_q025,xmax=null_localization_contrast_q975,ymin=0,ymax=line_top),inherit.aes=FALSE,fill="#DCE5EE",alpha=.55)
 p <- p + geom_histogram(aes(y=if(.env$density) after_stat(density) else after_stat(count)),bins=bins,fill=if(density) "#CED7E0" else "#CBD5DF",colour="white",linewidth=if(density) .18 else .2) +
   geom_segment(data=ln,aes(x=null_localization_contrast_median,xend=null_localization_contrast_median,y=0,yend=line_top,linetype="Null median"),inherit.aes=FALSE,colour="#666666",linewidth=.5) +
   geom_segment(data=ln,aes(x=null_localization_contrast_q025,xend=null_localization_contrast_q025,y=0,yend=line_top,linetype="Null 2.5–97.5 percentiles"),inherit.aes=FALSE,colour="#777777",linewidth=.45) +
   geom_segment(data=ln,aes(x=null_localization_contrast_q975,xend=null_localization_contrast_q975,y=0,yend=line_top,linetype="Null 2.5–97.5 percentiles"),inherit.aes=FALSE,colour="#777777",linewidth=.45) +
   geom_segment(data=ln,aes(x=observed_localization_contrast,xend=observed_localization_contrast,y=0,yend=line_top,colour=program,linetype="Observed contrast (module color)"),inherit.aes=FALSE,linewidth=.8) +
   geom_text(data=ln,aes(x=Inf,y=ann_y,label=ann),inherit.aes=FALSE,hjust=1.04,size=if(supplementary) 2.5 else 2.1,family=FONT,lineheight=1) +
   geom_blank(data=ln,aes(x=0,y=canvas_top),inherit.aes=FALSE) +
   facet_grid(module~threshold,scales="free_y",labeller=labeller(module=function(x) if(supplementary) x else sub(" ","\n",x,fixed=TRUE))) +
   scale_colour_manual(values=MODULE_COLORS,guide="none") +
   scale_linetype_manual(name=NULL,breaks=c("Observed contrast (module color)","Null median","Null 2.5–97.5 percentiles"),values=c("Observed contrast (module color)"=1,"Null median"=2,"Null 2.5–97.5 percentiles"=3)) +
   guides(linetype=guide_legend(ncol=1,byrow=TRUE,override.aes=list(colour=c(MODULE_COLORS[1],"#666666","#777777"),linewidth=c(.8,.5,.45)))) +
   scale_y_continuous(n.breaks=if(supplementary) 5 else 3,expand=expansion(mult=c(0,.01))) + theme_pub() +
   theme(legend.position="bottom",legend.key.height=unit(2.1,"mm"),legend.key.width=unit(7,"mm"),
         legend.spacing.y=unit(0,"mm"),legend.margin=margin(0,0,0,0),legend.box.spacing=unit(.6,"mm"),legend.text=element_text(size=6.5))
 if(supplementary) p <- p+geom_text(data=ln,aes(x=-Inf,y=ann_y,label=panel_letter),inherit.aes=FALSE,hjust=-.25,fontface="bold",size=3.6,family=FONT)
 p
}
