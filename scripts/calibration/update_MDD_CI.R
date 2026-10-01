d <- read.delim("processed_results/MDD_individual_models.tsv", check.names=FALSE)
stopifnot(nrow(d)==6, all(d$df_resid==30))
stopifnot(max(abs(2*pt(-abs(d$t),df=d$df_resid)-d$p_value))<1e-13)
lo <- d$beta_target_vs_Control-qt(.975,d$df_resid)*d$se
hi <- d$beta_target_vs_Control+qt(.975,d$df_resid)*d$se
stopifnot(max(abs(lo-d$ci_low))<1e-14,max(abs(hi-d$ci_high))<1e-14)
