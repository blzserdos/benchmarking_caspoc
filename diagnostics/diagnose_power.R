#!/usr/bin/env Rscript
# =============================================================================
# diagnose_power.R
# Decomposes where power is lost in the sim_signal experiment.
#
# Loss chain, measured stage by stage:
#   (1) pop_ceiling  s/(s+1)      -- generator's population ceiling, analytic
#   (2) oracle_big   true loadings, 1000 independent samples  -> checks (1)
#   (3) fit_big      fitted loadings (n=100), 1000 indep samples -> ESTIMATION loss
#   (4) cv_pooled    fitted per fold, all 100 held-out pooled  -> + CV train-size loss
#   (5) cv_meanfold  fitted per fold, mean of 10 corrs on n=10 -> + SMALL-FOLD noise
#
# (5) is the statistic naive_cv actually reports (at the true HPs).
# Running the same battery at s = 0 gives the null reference, so each statistic
# gets an effect size (mean_signal - mean_null)/sd_null and a predicted power.
# =============================================================================

suppressMessages({
  library(MASS)
  library(mixOmics)
  library(caret)
  library(parallel)
})

source("R/generate_data.R")

N          <- 100     # study sample size (matches config)
P          <- 200
Q          <- 50
N_REL_X    <- 20
N_REL_Y    <- 10
N_BIG      <- 1000    # independent evaluation set
NUM_FOLDS  <- 10
.args      <- commandArgs(trailingOnly = TRUE)
N_REPS     <- if (length(.args) >= 1) as.integer(.args[1]) else 40
STRENGTHS  <- if (length(.args) >= 2) {
  as.numeric(strsplit(.args[2], ",")[[1]])
} else {
  c(0, 0.25, 0.5, 1, 2, 4, 8, 16)
}
.tag       <- if (length(.args) >= 3) .args[3] else "diag_power"
OUT_FILE   <- file.path(dirname(sub("^--file=", "", grep("^--file=",
                commandArgs(FALSE), value = TRUE)[1])),
                paste0(.tag, ".rds"))

# Fit one sPLS component and return the loading pair.
fit_pair <- function(Xtr, Ytr, kx, ky) {
  fit <- tryCatch(
    spls(Xtr, Ytr, ncomp = 1, keepX = kx, keepY = ky,
         mode = "regression", scale = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  list(lx = fit$loadings$X[, 1], ly = fit$loadings$Y[, 1])
}

# Scale test block with train centre/scale, as the wrappers do.
scale_like <- function(test, train_scaled) {
  scale(test,
        center = attr(train_scaled, "scaled:center"),
        scale  = attr(train_scaled, "scaled:scale"))
}

run_rep <- function(rep_id, s) {
  # One draw of loadings shared by the study set and the big evaluation set.
  d <- generate_signal_data(
    n = N + N_BIG, p = P, q = Q,
    n_comp_true = 1, n_relevant_x = N_REL_X, n_relevant_y = N_REL_Y,
    signal_strength = s, seed = rep_id * 7919 + round(s * 1000)
  )
  idx_study <- seq_len(N)
  X  <- d$X[idx_study, , drop = FALSE]; Y  <- d$Y[idx_study, , drop = FALSE]
  Xb <- d$X[-idx_study, , drop = FALSE]; Yb <- d$Y[-idx_study, , drop = FALSE]
  lx_true <- d$true_loadings_x[, 1]; ly_true <- d$true_loadings_y[, 1]

  # (2) Oracle on the big set: true loadings, no estimation, no small folds.
  oracle_big <- cor(Xb %*% lx_true, Yb %*% ly_true, method = "spearman")[1, 1]

  # Oracle scored on 10-sample folds: isolates small-fold noise alone.
  set.seed(rep_id)
  folds <- createFolds(seq_len(N), k = NUM_FOLDS, list = TRUE)
  oracle_fold <- vapply(folds, function(te) {
    cor(X[te, ] %*% lx_true, Y[te, ] %*% ly_true, method = "spearman")[1, 1]
  }, numeric(1))

  # (3) Fitted on all 100, scored on the big set: estimation loss only.
  Xs <- scale(X); Ys <- scale(Y)
  full <- fit_pair(Xs, Ys, N_REL_X, N_REL_Y)
  if (is.null(full)) return(NULL)
  fit_big <- cor(scale_like(Xb, Xs) %*% full$lx,
                 scale_like(Yb, Ys) %*% full$ly, method = "spearman")[1, 1]

  sel_x <- which(full$lx != 0); sel_y <- which(full$ly != 0)
  recov_x <- length(intersect(sel_x, seq_len(N_REL_X))) / N_REL_X
  recov_y <- length(intersect(sel_y, seq_len(N_REL_Y))) / N_REL_Y

  # (4)+(5) 10-fold CV at the true HPs.
  fold_cors <- numeric(NUM_FOLDS)
  pooled_tx <- numeric(0); pooled_ty <- numeric(0)
  ref_lx <- NULL
  for (i in seq_len(NUM_FOLDS)) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- scale(X[tr, , drop = FALSE]); Ytr <- scale(Y[tr, , drop = FALSE])
    f <- fit_pair(Xtr, Ytr, N_REL_X, N_REL_Y)
    if (is.null(f)) { fold_cors[i] <- NA; next }
    tx <- scale_like(X[te, , drop = FALSE], Xtr) %*% f$lx
    ty <- scale_like(Y[te, , drop = FALSE], Ytr) %*% f$ly
    fold_cors[i] <- cor(tx, ty, method = "spearman")[1, 1]

    # Align the (lx, ly) pair to fold 1 before pooling; sPLS fixes the pair
    # only up to a joint sign flip.
    if (is.null(ref_lx)) ref_lx <- f$lx
    flip <- if (sum(f$lx * ref_lx) < 0) -1 else 1
    pooled_tx <- c(pooled_tx, flip * tx); pooled_ty <- c(pooled_ty, flip * ty)
  }

  data.frame(
    signal_strength = s, rep = rep_id,
    pop_ceiling  = s / (s + 1),
    oracle_big   = oracle_big,
    oracle_fold  = mean(oracle_fold, na.rm = TRUE),
    fit_big      = fit_big,
    cv_pooled    = cor(pooled_tx, pooled_ty, method = "spearman"),
    cv_meanfold  = mean(fold_cors, na.rm = TRUE),
    recov_x      = recov_x,
    recov_y      = recov_y
  )
}

grid <- expand.grid(rep_id = seq_len(N_REPS), s = STRENGTHS)
message(sprintf("Running %d cells on %d cores...", nrow(grid), detectCores() - 1))

res <- mclapply(seq_len(nrow(grid)),
                function(k) run_rep(grid$rep_id[k], grid$s[k]),
                mc.cores = max(1, detectCores() - 1))

bad <- !vapply(res, is.data.frame, logical(1))
if (any(bad)) {
  message(sprintf("%d cells failed. First error:\n%s",
                  sum(bad), as.character(res[bad][[1]])))
}
out <- do.call(rbind, res[!bad])
saveRDS(out, OUT_FILE)
message(sprintf("Done: %d rows -> %s", nrow(out), OUT_FILE))
