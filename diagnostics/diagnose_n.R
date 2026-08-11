#!/usr/bin/env Rscript
# =============================================================================
# diagnose_n.R
# Is sample size, rather than signal strength, the binding constraint?
# Holds signal_strength at the design point and varies n, measuring support
# recovery and the CV statistic against a matched null (s = 0) at each n.
# =============================================================================

suppressMessages({
  library(MASS); library(mixOmics); library(caret); library(parallel)
})
source("R/generate_data.R")

P <- 200; Q <- 50; N_REL_X <- 20; N_REL_Y <- 10
NUM_FOLDS <- 10
N_REPS    <- 100
N_VALUES  <- c(100, 150, 200, 300, 500)
S_VALUES  <- c(0, 1)
OUT_FILE  <- file.path(dirname(sub("^--file=", "", grep("^--file=",
               commandArgs(FALSE), value = TRUE)[1])), "diag_n.rds")

scale_like <- function(test, tr) {
  scale(test, center = attr(tr, "scaled:center"), scale = attr(tr, "scaled:scale"))
}

run_cell <- function(rep_id, n, s) {
  d <- generate_signal_data(
    n = n, p = P, q = Q, n_comp_true = 1,
    n_relevant_x = N_REL_X, n_relevant_y = N_REL_Y,
    signal_strength = s, seed = rep_id * 7919 + n * 13 + round(s * 1000)
  )
  X <- d$X; Y <- d$Y

  set.seed(rep_id)
  folds <- createFolds(seq_len(n), k = NUM_FOLDS, list = TRUE)
  fold_cors <- numeric(NUM_FOLDS)
  for (i in seq_len(NUM_FOLDS)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(n), te)
    Xtr <- scale(X[tr, , drop = FALSE]); Ytr <- scale(Y[tr, , drop = FALSE])
    fit <- tryCatch(spls(Xtr, Ytr, ncomp = 1, keepX = N_REL_X, keepY = N_REL_Y,
                         mode = "regression", scale = TRUE), error = function(e) NULL)
    if (is.null(fit)) { fold_cors[i] <- NA; next }
    tx <- scale_like(X[te, , drop = FALSE], Xtr) %*% fit$loadings$X[, 1]
    ty <- scale_like(Y[te, , drop = FALSE], Ytr) %*% fit$loadings$Y[, 1]
    fold_cors[i] <- cor(tx, ty, method = "spearman")[1, 1]
  }

  Xs <- scale(X); Ys <- scale(Y)
  full <- tryCatch(spls(Xs, Ys, ncomp = 1, keepX = N_REL_X, keepY = N_REL_Y,
                        mode = "regression", scale = TRUE), error = function(e) NULL)
  recov_x <- if (is.null(full)) NA else
    length(intersect(which(full$loadings$X[, 1] != 0), seq_len(N_REL_X))) / N_REL_X

  data.frame(n = n, signal_strength = s, rep = rep_id,
             cv_meanfold = mean(fold_cors, na.rm = TRUE),
             fold_size = round(n / NUM_FOLDS), recov_x = recov_x)
}

grid <- expand.grid(rep_id = seq_len(N_REPS), n = N_VALUES, s = S_VALUES)
message(sprintf("Running %d cells...", nrow(grid)))
res <- mclapply(seq_len(nrow(grid)),
                function(k) run_cell(grid$rep_id[k], grid$n[k], grid$s[k]),
                mc.cores = max(1, detectCores() - 1))
bad <- !vapply(res, is.data.frame, logical(1))
if (any(bad)) message(sprintf("%d failed. First: %s", sum(bad), as.character(res[bad][[1]])))
out <- do.call(rbind, res[!bad])
saveRDS(out, OUT_FILE)
message(sprintf("Done: %d rows -> %s", nrow(out), OUT_FILE))
