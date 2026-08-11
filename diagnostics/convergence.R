#!/usr/bin/env Rscript
# How often does spls fail to converge, and does the structure cause it?
suppressMessages({ library(mixOmics); library(caret); library(parallel) })
source("R/generate_data.R"); source("cluster/config.R")
options(width = 140)

N <- 100; REPS <- 60
CASES <- list(independent = NULL, structured = BLOCK_STRUCTURE)
SS <- c(0, 2, 8)
KX <- c(5, 10, 20, 50); KY <- c(5, 10, 20)

probe <- function(rep_id, s, cn) {
  d <- generate_signal_data(N, 200, 50, 1, 20, 10, s, structure = CASES[[cn]],
                            signal_alignment = SIGNAL_ALIGNMENT, seed = rep_id * 91 + round(s * 10))
  set.seed(rep_id); folds <- createFolds(seq_len(N), k = 10, list = TRUE)
  nwarn <- 0L; nfit <- 0L; nerr <- 0L
  for (kx in KX) for (ky in KY) for (i in 1:3) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- scale(d$X[tr, , drop = FALSE]); Ytr <- scale(d$Y[tr, , drop = FALSE])
    nfit <- nfit + 1L
    r <- withCallingHandlers(
      tryCatch(spls(Xtr, Ytr, ncomp = 1, keepX = kx, keepY = ky,
                    mode = "regression", scale = TRUE), error = function(e) NULL),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w), ignore.case = TRUE)) nwarn <<- nwarn + 1L
        invokeRestart("muffleWarning")
      })
    if (is.null(r)) nerr <- nerr + 1L
  }
  data.frame(cfg = cn, s = s, fits = nfit, warns = nwarn, errs = nerr)
}

jobs <- expand.grid(rep_id = seq_len(REPS), s = SS, cn = names(CASES), stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(k) probe(jobs$rep_id[k], jobs$s[k], jobs$cn[k]),
                mc.cores = max(1, detectCores() - 1))
out <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])
agg <- aggregate(cbind(fits, warns, errs) ~ cfg + s, out, sum)
agg$pct_nonconverged <- round(100 * agg$warns / agg$fits, 2)
agg$pct_error <- round(100 * agg$errs / agg$fits, 2)
cat("=== spls convergence, ", sum(out$fits), " fits ===\n", sep = "")
print(agg[order(agg$cfg, agg$s), ], row.names = FALSE)
