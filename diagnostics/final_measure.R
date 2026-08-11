#!/usr/bin/env Rscript
# Final measurements under the fitted decaying structure:
#   (a) how loosely signal_strength maps to the oracle correlation
#   (b) where the power transition sits -> places SIGNAL_STRENGTHS
suppressMessages({ library(mixOmics); library(caret); library(parallel) })
source("R/generate_data.R"); source("cluster/config.R")
options(width = 150)
SP <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))

cat("=== (a) oracle cross-block correlation vs the s/(s+1) identity ===\n")
for (al in c("aligned", "spread")) {
  for (s in c(1, 2, 4)) {
    v <- vapply(1:200, function(r) {
      d <- generate_signal_data(3000, 200, 50, 1, 20, 10, s, structure = BLOCK_STRUCTURE,
                                signal_alignment = al, seed = r)
      cor(d$X %*% d$true_loadings_x[, 1], d$Y %*% d$true_loadings_y[, 1])[1, 1]
    }, numeric(1))
    cat(sprintf("  %-8s s=%-2g  s/(s+1)=%.3f  observed mean=%.3f sd=%.3f  IQR=[%.2f, %.2f]\n",
                al, s, s / (s + 1), mean(v), sd(v), quantile(v, .25), quantile(v, .75)))
  }
}

cat("\n=== (b) power transition under the fitted structure ===\n")
N <- 100; NF <- 10; REPS <- 200
SS <- c(0, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8)
sl <- function(te, tr) scale(te, center = attr(tr, "scaled:center"), scale = attr(tr, "scaled:scale"))
cell <- function(rep_id, s, al) {
  d <- generate_signal_data(N, 200, 50, 1, 20, 10, s, structure = BLOCK_STRUCTURE,
                            signal_alignment = al, seed = rep_id * 7919 + round(s * 100))
  set.seed(rep_id); folds <- createFolds(seq_len(N), k = NF, list = TRUE)
  fc <- vapply(seq_len(NF), function(i) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- scale(d$X[tr, , drop = FALSE]); Ytr <- scale(d$Y[tr, , drop = FALSE])
    f <- tryCatch(spls(Xtr, Ytr, ncomp = 1, keepX = 20, keepY = 10,
                       mode = "regression", scale = TRUE), error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    cor(sl(d$X[te, , drop = FALSE], Xtr) %*% f$loadings$X[, 1],
        sl(d$Y[te, , drop = FALSE], Ytr) %*% f$loadings$Y[, 1], method = "spearman")[1, 1]
  }, numeric(1))
  data.frame(al = al, s = s, stat = mean(fc, na.rm = TRUE))
}
jobs <- expand.grid(rep_id = seq_len(REPS), s = SS, al = c("aligned", "spread"),
                    stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(k) cell(jobs$rep_id[k], jobs$s[k], jobs$al[k]),
                mc.cores = max(1, detectCores() - 1))
out <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])
for (al in c("aligned", "spread")) {
  a <- out[out$al == al, ]; thr <- quantile(a$stat[a$s == 0], 0.95, na.rm = TRUE)
  cat(sprintf("  %-8s null sd=%.3f thr=%.3f\n", al, sd(a$stat[a$s == 0], na.rm = TRUE), thr))
  cat("    s     :", sprintf("%6g", SS[-1]), "\n")
  cat("    power :", sprintf("%6.2f", vapply(SS[-1], function(s)
        mean(a$stat[a$s == s] > thr, na.rm = TRUE), numeric(1))), "\n")
}
saveRDS(out, file.path(SP, "final_measure.rds"))
