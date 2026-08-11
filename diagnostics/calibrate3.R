suppressMessages({library(mixOmics); library(caret); library(dplyr); library(tibble); library(parallel)})
source("R/generate_data.R"); source("cluster/config.R")
src <- readLines("R/cv_approaches.R")
eval(parse(text=paste(src[!grepl("^library[(]caspoc[)]", src)], collapse="\n")), envir=globalenv())
SP <- dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
N_SUB<-100; N_REPS<-15; KEEPX<-c(5,10,20,50); KEEPY<-c(5,10,20)
STRENGTHS<-c(0,2,4,8,16,32,64)
data("breast.TCGA"); blocks<-breast.TCGA$data.train[c("mirna","mrna","protein")]
prs<-list(c("mrna","protein"),c("mirna","protein"),c("mirna","mrna"))
stat_of<-function(X,Y,seed){ r<-tryCatch(run_naive_cv(X,Y,ncomp=1,num_folds=10,
  keepX_options=KEEPX,keepY_options=KEEPY,seed=seed),error=function(e) NULL)
  if(is.null(r)) NA_real_ else r$observed_stat }
jobs<-expand.grid(rep_id=seq_len(N_REPS), s=STRENGTHS, pair=seq_along(prs))
sim<-mclapply(seq_len(nrow(jobs)),function(k){
  pr<-prs[[jobs$pair[k]]]; rid<-jobs$rep_id[k]; s<-jobs$s[k]
  p<-ncol(blocks[[pr[1]]]); q<-ncol(blocks[[pr[2]]])
  d<-generate_signal_data(N_SUB,p,q,1,20,10,signal_strength=s,structure=BLOCK_STRUCTURE,
                          signal_alignment=SIGNAL_ALIGNMENT,seed=rid*7919+round(s*1000)+p)
  data.frame(label=paste(pr,collapse=" ~ "),signal_strength=s,stat=stat_of(d$X,d$Y,seed=rid))
},mc.cores=max(1,detectCores()-1))
sim<-do.call(rbind,sim[vapply(sim,is.data.frame,logical(1))])
real<-readRDS(file.path(SP,"calib_s.rds"))$real %>% group_by(label) %>%
  summarise(real_stat=mean(stat,na.rm=TRUE),.groups="drop")
ss<-sim %>% group_by(label,signal_strength) %>% summarise(v=mean(stat,na.rm=TRUE),.groups="drop")
cat("=== statistic vs s (final structure, spread) ===\n")
for (L in unique(ss$label)) { a<-ss[ss$label==L,]; a<-a[order(a$signal_strength),]
  cat(sprintf("%-18s",L), sprintf("%6.3f",a$v), "\n") }
cat(sprintf("%-18s","s ="), sprintf("%6g",sort(unique(ss$signal_strength))), "\n")
cat("\n=== matched s ===\n")
print(do.call(rbind, lapply(unique(ss$label), function(L){
  a<-ss[ss$label==L,]; a<-a[order(a$signal_strength),]; t<-real$real_stat[real$label==L]
  s_hat<-if(t>=max(a$v)) Inf else if(t<=min(a$v)) NA_real_ else approx(a$v,a$signal_strength,xout=t,ties="ordered")$y
  data.frame(pair=L, real_stat=round(t,3), matched_s=round(s_hat,1))})), row.names=FALSE)
