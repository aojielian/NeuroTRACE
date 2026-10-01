source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/theme.R"))
source(file.path(Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd()),"scripts/figures/null_display.R"))
s3 <- read_panel("Supplementary_FigureS3_empirical_null_distributions.tsv")
pS3 <- null_display(s3,bins=32,density=FALSE,supplementary=TRUE) +
 labs(x="Prenatal − postnatal contrast",y="Null replicate count") +
 theme(strip.text=element_text(size=8),panel.spacing=unit(2.4,"mm"),
       plot.margin=margin(3,3,3,3,"mm"))
save_pub(pS3,"Supplementary_FigureS3",180,155,".")
