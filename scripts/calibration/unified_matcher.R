#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Matrix); library(data.table)})
options(stringsAsFactors=FALSE)
ROOT <- Sys.getenv("NEUROTRACE_PROJECT_ROOT"); stopifnot(nzchar(ROOT), dir.exists(ROOT))
RUN <- Sys.getenv("NEUROTRACE_OUTPUT_RUN"); stopifnot(nzchar(RUN))
OUT <- file.path(RUN, "output/03_matching")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
N_PERM <- 1000L; BASE_SEED <- 202609301L
ALPHA <- 0.35; TOL <- 1e-10; MAX_ITER <- 120L
STAGE_ORDER <- c("early_prenatal","mid_prenatal","late_prenatal","early_postnatal","childhood_adolescence","adulthood")
PROGRAMS <- c("NTM1_ASD_up","NTM2_ASD_down","NTM3_ASD_signed"); TOPS <- c(200L,500L)
clean <- function(x) toupper(trimws(as.character(x)))
# Rank-bin rule: rank average, then ceil(rank/max*5), bounded 1..5.
make_bins <- function(x, n_bins=5L) {
  r <- rank(x, ties.method="average", na.last="keep"); out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
  if (any(ok)) out[ok] <- pmax(1L, pmin(n_bins, ceiling(r[ok]/max(r[ok])*n_bins)))
  out
}
bh <- function(p) { n <- length(p); o <- order(p); out <- numeric(n); prev <- 1; for (k in n:1) {i<-o[k]; prev<-min(prev,p[i]*n/k);out[i]<-min(1,prev)}; out }
qsafe <- function(x,p) as.numeric(quantile(x,p,type=7,names=FALSE,na.rm=TRUE))
# Inputs: processed NeuroTRACE reference and fixed graph.
node_file <- file.path(ROOT,"processed_inputs/graph_nodes.tsv")
edge_file <- file.path(ROOT,"processed_inputs/graph_edges.tsv.gz")
prof_file <- file.path(ROOT,"processed_inputs/BrainSpan_stage_profiles.tsv.gz")
wt_file <- file.path(ROOT,"processed_inputs/native_module_weights.tsv")
obs_file <- file.path(ROOT,"processed_inputs/observed_developmental_contrasts.tsv")
nodes <- fread(node_file, showProgress=FALSE)
gn <- nodes[node_type=="gene"]
universe <- clean(gn$feature_id); universe <- universe[!duplicated(universe)]
stopifnot(length(universe)==5000L)
idx <- setNames(seq_along(universe), universe)
# Stage profiles aligned to frozen graph universe.
prof <- fread(prof_file, showProgress=FALSE); prof[, sym:=clean(gene_key)]
Z <- matrix(0, nrow=6, ncol=length(universe), dimnames=list(STAGE_ORDER,universe))
for (si in seq_along(STAGE_ORDER)) {
  z <- prof[state_id==STAGE_ORDER[si]]
  jj <- idx[clean(z$gene_key)]; ok <- !is.na(jj)
  Z[si, jj[ok]] <- as.numeric(z$stage_profile_z[ok])
}
stopifnot(all(rowSums(abs(Z))>0))
expr_mean <- colMeans(Z); expr_var <- apply(Z,2,var)
expr_bin <- make_bins(expr_mean); var_bin <- make_bins(expr_var)
gene_stats <- data.table(gene=universe, expr_bin=expr_bin, var_bin=var_bin)
# Symmetric gene-only graph, row-normalized exactly as frozen Stage4.
edges <- fread(edge_file, showProgress=FALSE)
edges <- edges[edge_type=="gene_gene_embedding_knn"]
si <- idx[clean(sub("^gene:","",edges[["source"]]))]; sj <- idx[clean(sub("^gene:","",edges[["target"]]))]
ok <- !is.na(si)&!is.na(sj)&is.finite(edges[["weight"]])&edges[["weight"]]>0
si<-si[ok];sj<-sj[ok];ww<-as.numeric(edges$weight[ok])
A <- sparseMatrix(i=c(si,sj),j=c(sj,si),x=c(ww,ww),dims=c(5000L,5000L),giveCsparse=TRUE)
rs <- as.numeric(rowSums(A)); inv <- rep(0,length(rs)); inv[rs>0] <- 1/rs[rs>0]
P <- Diagonal(x=inv) %*% A
# Solve the identical PPR fixed point in batch; contraction 0.65 makes the residual
# equivalent to <=120 power iterations at tol=1e-10, and is checked below.
M <- Diagonal(n=5000L) - (1-ALPHA)*t(P)
LU <- lu(M)
# Native module blocks, mapped to frozen universe and sorted by rank.
wt <- fread(wt_file, showProgress=FALSE); wt[, gene:=clean(gene_symbol_fixed)]
obs <- fread(obs_file, showProgress=FALSE)
get_module <- function(program, topn) {
  d <- wt[program==program & top_n==topn]
  d <- d[order(rank_within_program)]
  d <- d[gene %in% universe]
  # exact mapped set used by frozen stage4; de-duplicate conservatively by rank
  d <- d[!duplicated(gene)]
  list(genes=d$gene, weights=as.numeric(d$weight))
}
# Avoid NSE collision in get_module.
get_module <- function(pr, tn) { d<-wt[wt$program==pr & wt$top_n==tn]; d<-d[order(d$rank_within_program)]; d<-d[d$gene %in% universe]; d<-d[!duplicated(d$gene)]; list(genes=d$gene,weights=as.numeric(d$weight)) }
module_stage <- function(q) {
  nq <- sqrt(sum(q*q)); nz <- sqrt(rowSums(Z*Z)); as.numeric((Z %*% q)/(nz*nq))
}
batch_stage <- function(Q) {
  nq <- sqrt(colSums(Q*Q)); nz <- sqrt(rowSums(Z*Z)); C <- t(Q) %*% t(Z); C <- sweep(C,2,nz,"/"); C <- sweep(C,1,nq,"/"); C
}
solve_channels <- function(Spos,Sneg) {
  Qp <- if (sum(abs(Spos))>0) as.matrix(solve(LU, ALPHA*Spos)) else matrix(0,nrow(Spos),ncol(Spos))
  Qn <- if (sum(abs(Sneg))>0) as.matrix(solve(LU, ALPHA*Sneg)) else matrix(0,nrow(Sneg),ncol(Sneg))
  list(q=Qp+Qn, qp=Qp, qn=Qn)
}
# Deterministic independent combo streams: base seed plus fixed stride.
all_rows <- list(); draw_rows <- list(); replicate_rows <- list(); combo_i <- 0L
for (pr in PROGRAMS) for (tn in TOPS) {
  combo_i <- combo_i+1L; combo_seed <- BASE_SEED + combo_i*1000003L; set.seed(combo_seed)
  mod <- get_module(pr,tn); mg <- mod$genes; mw <- mod$weights; m <- length(mg)
  if (m==0) stop("empty mapped module")
  # observed native channels
  ow <- numeric(5000); ow[match(mg,universe)] <- mw
  Sop <- matrix(0,5000,1); Sne <- matrix(0,5000,1); ii<-which(ow>0); jj<-which(ow<0)
  if(length(ii)) Sop[ii,1]<-abs(ow[ii])/sum(abs(ow[ii])); if(length(jj)) Sne[jj,1]<-abs(ow[jj])/sum(abs(ow[jj]))
  qobs <- solve_channels(Sop,Sne)$q[,1]; cobs <- module_stage(qobs); obs_contrast <- mean(cobs[1:3])-mean(cobs[4:6])
  Spos <- matrix(0,5000,N_PERM); Sneg <- matrix(0,5000,N_PERM); tier_count <- c(expr_var_bin=0L,expr_bin_only=0L,var_bin_only=0L,any_nonmodule=0L,fallback=0L)
  module_mask <- universe %in% mg; chosen_all <- vector("list",N_PERM)
  for (b in seq_len(N_PERM)) {
    chosen <- character(m); tiers <- character(m); used <- logical(5000)
    for (k in seq_len(m)) {
      g0 <- mg[k]; eb <- expr_bin[idx[[g0]]]; vb <- var_bin[idx[[g0]]]
      pool1 <- which(!module_mask & !used & expr_bin==eb & var_bin==vb)
      pool2 <- which(!module_mask & !used & expr_bin==eb)
      pool3 <- which(!module_mask & !used & var_bin==vb)
      pool4 <- which(!module_mask & !used)
      if(length(pool1)) {pool<-pool1;tier<-"expr_var_bin"} else if(length(pool2)){pool<-pool2;tier<-"expr_bin_only"} else if(length(pool3)){pool<-pool3;tier<-"var_bin_only"} else if(length(pool4)){pool<-pool4;tier<-"any_nonmodule"} else {pool<-seq_len(5000);tier<-"fallback"}
      pick <- sample(pool,1L); chosen[k]<-universe[pick]; used[pick]<-TRUE; tiers[k]<-tier; tier_count[tier]<-tier_count[tier]+1L
    }
    wp <- sample(mw, size=m, replace=FALSE); jj<-idx[chosen]; pos<-which(wp>0); neg<-which(wp<0)
    if(length(pos)) Spos[jj[pos],b] <- abs(wp[pos])/sum(abs(wp[pos])); if(length(neg)) Sneg[jj[neg],b] <- abs(wp[neg])/sum(abs(wp[neg]))
    chosen_all[[b]] <- paste(chosen,collapse=",")
    if(b<=3L) draw_rows[[length(draw_rows)+1L]] <- data.table(program=pr,top_n=tn,replicate=b,combo_seed=combo_seed,genes=chosen_all[[b]],tier_counts=paste(names(table(tiers)),as.integer(table(tiers)),sep=":",collapse=";"))
  }
  sol <- solve_channels(Spos,Sneg); C <- batch_stage(sol$q); nullc <- rowMeans(C[,1:3])-rowMeans(C[,4:6])
  # residual check against PPR linear equations; equivalent fixed point convergence evidence.
  rp <- if(sum(abs(sol$qp))>0) max(colSums(abs(M %*% sol$qp - ALPHA*Spos))) else 0
  rn <- if(sum(abs(sol$qn))>0) max(colSums(abs(M %*% sol$qn - ALPHA*Sneg))) else 0
  rpv <- if(sum(abs(sol$qp))>0) colSums(abs(M %*% sol$qp - ALPHA*Spos)) else rep(0,N_PERM)
  rnv <- if(sum(abs(sol$qn))>0) colSums(abs(M %*% sol$qn - ALPHA*Sneg)) else rep(0,N_PERM)
  rr <- data.table(program=pr,top_n=tn,replicate=seq_len(N_PERM),localization_contrast=nullc,null_top_stage=STAGE_ORDER[max.col(C,ties.method="first")],positive_residual=rpv,negative_residual=rnv,converged=pmax(rpv,rnv)<=TOL,combo_seed=combo_seed)
  for(si in seq_along(STAGE_ORDER)) rr[[paste0("cosine_",STAGE_ORDER[si])]] <- C[,si]
  replicate_rows[[length(replicate_rows)+1L]] <- rr
  ph <- (1+sum(nullc>=obs_contrast))/(N_PERM+1); pl <- (1+sum(nullc<=obs_contrast))/(N_PERM+1); pt <- min(1,2*min(ph,pl))
  all_rows[[length(all_rows)+1L]] <- data.table(program=pr,top_n=tn,n_mapped_genes=m,combo_seed=combo_seed,observed_contrast_frozen=as.numeric(obs[program==pr & top_n==tn]$localization_contrast[1]),observed_contrast_unified=obs_contrast,observed_delta_unified_minus_frozen=obs_contrast-as.numeric(obs[program==pr & top_n==tn]$localization_contrast[1]),null_median=median(nullc),null_q025=qsafe(nullc,.025),null_q975=qsafe(nullc,.975),null_min=min(nullc),null_max=max(nullc),p_unified=pt,p_high=ph,p_low=pl,n_ge_obs=sum(nullc>=obs_contrast),n_le_obs=sum(nullc<=obs_contrast),fallback_fraction=unname(tier_count["fallback"])/(m*N_PERM),expr_var_fraction=unname(tier_count["expr_var_bin"])/(m*N_PERM),expr_fraction=unname(tier_count["expr_bin_only"])/(m*N_PERM),var_fraction=unname(tier_count["var_bin_only"])/(m*N_PERM),any_fraction=unname(tier_count["any_nonmodule"])/(m*N_PERM),max_positive_residual=rp,max_negative_residual=rn,solver="sparse_lu_fixed_point_equivalent_to_120_iter",convergence_status=ifelse(max(rp,rn)<=TOL,"PASS_RESIDUAL<=1e-10","CHECK_RESIDUAL"),B=N_PERM,seed=BASE_SEED)
  cat(sprintf("[unified] %s top%d: mapped=%d obs=%.9f null median=%.9f P=%.6f residual=(%.3e,%.3e)\n",pr,tn,m,obs_contrast,median(nullc),pt,rp,rn),file=stderr())
}
res <- rbindlist(all_rows); res[,fdr_unified:=bh(p_unified)]; setcolorder(res,c("program","top_n","n_mapped_genes","combo_seed","observed_contrast_frozen","observed_contrast_unified","observed_delta_unified_minus_frozen","null_median","null_q025","null_q975","null_min","null_max","p_unified","fdr_unified","p_high","p_low","n_ge_obs","n_le_obs","expr_var_fraction","expr_fraction","var_fraction","any_fraction","fallback_fraction","max_positive_residual","max_negative_residual","solver","convergence_status","B","seed"))
fwrite(res,file.path(OUT,"unified_matching_sensitivity_results.tsv"),sep="\t",quote=FALSE,na="NA")
if(length(draw_rows)) fwrite(rbindlist(draw_rows),file.path(OUT,"unified_matching_sensitivity_draw_manifest_first3.tsv"),sep="\t",quote=FALSE)
writeLines(c("Unified matching sensitivity completed.",sprintf("B=%d; base seed=%d; combo streams seed + combo_index*1000003",N_PERM,BASE_SEED),sprintf("max PPR fixed-point residual: %.3e",max(res$max_positive_residual,res$max_negative_residual))),file.path(OUT,"unified_matching_sensitivity_run.log"))

fwrite(rbindlist(replicate_rows),file.path(OUT,"unified_matching_null_replicates.tsv.gz"),sep="\t",quote=FALSE,compress="gzip")
