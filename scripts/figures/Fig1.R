source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
# Schematic geometry only. Statistical sources and propagation settings remain frozen.
arrow_std <- arrow(length=unit(1.7,"mm"),type="closed")
nodes <- expand.grid(ix=1:5,iy=1:4) %>% mutate(x=3+(ix-1)*.92+ifelse(iy%%2==0,.13,0),y=3.7+(iy-1)*.55,id=row_number())
edges <- bind_rows(lapply(1:(nrow(nodes)-1),function(i)data.frame(x=nodes$x[i],y=nodes$y[i],xend=nodes$x[i+1],yend=nodes$y[i+1])),lapply(1:(nrow(nodes)-5),function(i)data.frame(x=nodes$x[i],y=nodes$y[i],xend=nodes$x[i+5],yend=nodes$y[i+5])))
pA <- ggplot() +
 annotate("rect",xmin=.8,xmax=9.2,ymin=8.6,ymax=10.8,fill="#FAFAFA",colour="#777777",linewidth=.5) +
 annotate("text",x=5,y=10.4,label="Signed gene weights within a module",family=FONT,fontface="bold",size=3) +
 annotate("segment",x=c(2.3,4.1,5.9,7.7),xend=c(2.3,4.1,5.9,7.7),y=9.25,yend=c(9.88,9.65,8.94,8.75),colour=c("#7B3294","#7B3294","#2166AC","#2166AC"),linewidth=4) +
 annotate("segment",x=1.65,xend=8.35,y=9.25,yend=9.25,colour="#777777",linewidth=.35) +
 annotate("text",x=c(2.3,4.1,5.9,7.7),y=c(10,9.8,9.68,9.68),label=c("+","+","−","−"),family=FONT,size=3.1,fontface="bold",colour=c("#7B3294","#7B3294","#2166AC","#2166AC")) +
 annotate("rect",xmin=.6,xmax=4.3,ymin=6.6,ymax=8,fill="#F3EAF6",colour="#7B3294",linewidth=.6) +
 annotate("rect",xmin=5.7,xmax=9.4,ymin=6.6,ymax=8,fill="#E9F0F8",colour="#2166AC",linewidth=.6) +
 annotate("text",x=2.45,y=7.55,label="Increased expression",family=FONT,fontface="bold",colour="#7B3294",size=2.85) +
 annotate("text",x=7.55,y=7.55,label="Decreased expression",family=FONT,fontface="bold",colour="#2166AC",size=2.85) +
 annotate("text",x=c(2.45,7.55),y=6.98,label=c("Restart: |w|, w > 0","Restart: |w|, w < 0"),family=FONT,size=2.65) +
 annotate("segment",x=c(3,7),xend=c(2.45,7.55),y=8.6,yend=8.08,arrow=arrow_std,linewidth=.5,colour="#555555") +
 annotate("segment",x=c(2.45,7.55),xend=c(2.9,7.1),y=6.6,yend=5.58,arrow=arrow_std,linewidth=.5,colour="#555555") +
 geom_segment(data=edges,aes(x=x,y=y,xend=xend,yend=yend),colour="#C6C6C6",linewidth=.4) +
 geom_point(data=nodes,aes(x=x,y=y),shape=21,size=2.4,fill="#F3F3F3",colour="#777777",stroke=.45) +
 annotate("text",x=5,y=5.98,label="Gene-only developmental\nsimilarity network",family=FONT,fontface="bold",size=2.7,lineheight=.93) +
 annotate("segment",x=c(3.3,6.7),xend=c(2.45,7.55),y=3.55,yend=2.98,arrow=arrow_std,linewidth=.5,colour="#555555") +
 annotate("rect",xmin=.6,xmax=4.3,ymin=1.4,ymax=2.9,fill="#F3EAF6",colour="#7B3294",linewidth=.6) +
 annotate("rect",xmin=5.7,xmax=9.4,ymin=1.4,ymax=2.9,fill="#E9F0F8",colour="#2166AC",linewidth=.6) +
 annotate("text",x=2.45,y=2.37,label="Increased channel\n(q_pos)",family=FONT,fontface="bold",colour="#7B3294",size=3,lineheight=.95) +
 annotate("text",x=7.55,y=2.37,label="Decreased channel\n(q_neg)",family=FONT,fontface="bold",colour="#2166AC",size=3,lineheight=.95) +
 annotate("text",x=c(2.45,7.55),y=1.65,label="Gene-level PPR",family=FONT,size=2.5) +
 annotate("text",x=5,y=.8,label="alpha = 0.35  ·  tol = 10⁻¹⁰  ·  max_iter = 120",family=FONT,size=2.5,colour="#555555") +
 annotate("text",x=5,y=.3,label="Schematic",family=FONT,size=2.25,colour="#777777") +
 coord_cartesian(xlim=c(0,10),ylim=c(0,11),clip="off")+theme_void(base_family=FONT)
# Six stage bars occur only in the combined-magnitude developmental readout.
bars <- data.frame(x=1:6,y=c(.18,.25,.55,.31,.16,.11),stage=STAGE_ORDER)
pB <- ggplot() +
 annotate("text",x=3.5,y=7.9,label="Dual-head outputs (schematic)",family=FONT,fontface="bold",size=2.85) +
 annotate("text",x=3.5,y=7.25,label="Combined magnitude",family=FONT,fontface="bold",size=2.75) +
 annotate("text",x=3.5,y=6.72,label="q_magnitude = q_pos + q_neg",family=FONT,size=2.45) +
 annotate("text",x=3.5,y=6.2,label="Developmental localization",family=FONT,fontface="bold",size=2.6) +
 geom_rect(data=bars,aes(xmin=x-.33,xmax=x+.33,ymin=4.8,ymax=4.8+y*2.1,fill=stage),colour=NA) +
 scale_fill_manual(values=STAGE_COLORS,guide="none") +
 annotate("text",x=1:6,y=4.4,label=c("EPr","MPr","LPr","EPost","C/A","Adult"),family=FONT,size=2.1,lineheight=.9,colour="#555555") +
 annotate("segment",x=.4,xend=6.6,y=3.95,yend=3.95,colour="#D0D0D0",linewidth=.45) +
 annotate("text",x=3.5,y=3.55,label="Signed orientation",family=FONT,fontface="bold",size=2.75) +
 annotate("text",x=3.5,y=3.05,label="q_signed = q_pos − q_neg",family=FONT,size=2.45) +
 annotate("text",x=c(1.65,5.0),y=2.38,label=c("Gene priority","External sample\nscoring"),family=FONT,fontface="bold",size=2.5,lineheight=.92) +
 annotate("segment",x=c(.8,1.3,1.8,2.3),xend=c(.8,1.3,1.8,2.3),y=1.2,yend=c(1.85,1.6,.8,.65),colour=c("#7B3294","#7B3294","#2166AC","#2166AC"),linewidth=2.3) +
 annotate("segment",x=.45,xend=2.65,y=1.2,yend=1.2,colour="#999999",linewidth=.35) +
 annotate("segment",x=3.0,xend=3.7,y=1.2,yend=1.2,arrow=arrow_std,colour="#555555",linewidth=.45) +
 annotate("segment",x=4,xend=6.15,y=1.2,yend=1.2,colour="#999999",linewidth=.4) +
 annotate("point",x=c(4.4,5.65),y=1.2,shape=21,fill=c("#2166AC","#7B3294"),colour=c("#2166AC","#7B3294"),size=2.3) +
 annotate("text",x=c(4.35,5.7),y=.7,label=c("−","+"),colour=c("#2166AC","#7B3294"),family=FONT,size=3) +
 annotate("text",x=1.6,y=.25,label="|q_signed|",family=FONT,size=2.35) +
 coord_cartesian(xlim=c(0,7),ylim=c(0,8.2),clip="off")+theme_void(base_family=FONT)
pC <- ggplot() +
 annotate("text",x=0,y=5.9,label="Definitions and data flow",hjust=0,family=FONT,fontface="bold",size=2.85) +
 annotate("text",x=0,y=5.28,label="Localization: cosine(q_magnitude, p_t)",hjust=0,family=FONT,size=2.5) +
 annotate("text",x=0,y=4.75,label="Gene priority: |q_signed(g)|",hjust=0,family=FONT,size=2.5) +
 annotate("text",x=0,y=4.06,label="Modules",hjust=0,family=FONT,fontface="bold",size=2.75) +
 annotate("text",x=.1,y=c(3.5,3,2.5),label=unname(MODULE_LABELS),hjust=0,family=FONT,size=2.7,colour=unname(MODULE_COLORS)) +
 annotate("text",x=0,y=1.75,label="External cohorts",hjust=0,family=FONT,fontface="bold",size=2.75) +
 annotate("text",x=.1,y=1.13,label="GSE102741: 13 ASD / 39 control",hjust=0,family=FONT,size=2.5) +
 annotate("text",x=.1,y=.55,label="GSE64018: 12 ASD / 12 control",hjust=0,family=FONT,size=2.5) +
 coord_cartesian(xlim=c(0,7),ylim=c(0,6.2),clip="off")+theme_void(base_family=FONT)
fig <- tag_panel(pA,"A") | (tag_panel(pB,"B") / tag_panel(pC,"C") + plot_layout(heights=c(.57,.43)))
fig <- fig + plot_layout(widths=c(.62,.38))
save_pub(fig,"Fig1",180,118,".")
