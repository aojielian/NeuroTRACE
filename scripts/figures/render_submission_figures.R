ROOT <- Sys.getenv("NEUROTRACE_PROJECT_ROOT", unset=getwd())
for (n in c("Fig1","Fig2","Fig3","Fig4","Fig5","S2","S3","S5","S7")) source(file.path(ROOT,"scripts/figures",paste0(n,".R")))
