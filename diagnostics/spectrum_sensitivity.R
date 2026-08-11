#!/usr/bin/env Rscript
# =============================================================================
# spectrum_sensitivity.R
# Does the "3 spikes then a cliff" artifact of equal-rho block compound symmetry
# change the benchmark's behaviour, versus a structure with a realistically
# decaying eigenvalue tail?
#
# Decaying variant: K modules, sizes ~ k^-alpha, per-module correlation
# rho_k = rho1 * k^-beta. Same factor construction, so marginal variances stay 1
# and signal_strength keeps its meaning.
# =============================================================================

suppressMessages({ library(mixOmics); library(caret); library(parallel) })
SP <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))

msize <- function(p, K, alpha, frac) {
  w <- (seq_len(K))^(-alpha); m <- pmax(2, round(w / sum(w) * frac * p))
  while (sum(m) > p) m[which.max(m)] <- m[which.max(m)] - 1
  m
}
# rho_k = rho1 * k^-beta  (beta = 0 and K = 3 reproduces the current structure)
noise_decay <- function(n, p, K, alpha, rho1, beta, frac) {
  m <- msize(p, K, alpha, frac); E <- matrix(rnorm(n * p), n, p); pos <- 0L
  for (k in seq_along(m)) {
    idx <- (pos + 1L):(pos + m[k]); pos <- pos + m[k]
    rk <- min(0.95, rho1 * k^(-beta))
    E[, idx] <- sqrt(rk) * rnorm(n) + sqrt(1 - rk) * E[, idx]
  }
  E
}
spec <- function(M, k = 20) { e <- prcomp(scale(M), center = FALSE, scale. = FALSE)$sdev^2
  (e / sum(e))[seq_len(k)] }

# --- Target: real mRNA spectrum ----------------------------------------------
data("breast.TCGA")
ev_real <- spec(as.matrix(breast.TCGA$data.train$mrna))

# --- Fit the decaying variant to the FULL first-20 spectrum -------------------
grid <- expand.grid(K = c(8, 12, 16, 20, 30), alpha = c(0.5, 1, 1.5, 2),
                    rho1 = c(0.5, 0.6, 0.7, 0.8), beta = c(0.3, 0.5, 0.8, 1.2))
fit <- mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  set.seed(i * 31 + 7)
  ev <- rowMeans(replicate(10,
    spec(noise_decay(150, 200, g$K, g$alpha, g$rho1, g$beta, 1.0))))
  cbind(g, err = sqrt(mean((ev - ev_real)^2)) * 100)
}, mc.cores = max(1, detectCores() - 1))
fit <- do.call(rbind, fit[vapply(fit, is.data.frame, logical(1))])
fit <- fit[order(fit$err), ]
cat("=== best decaying-spectrum fits (RMSE over PC1-20, in % points) ===\n")
print(round(head(fit, 5), 3), row.names = FALSE)
best <- fit[1, ]

cur <- c(K = 3, alpha = 0.25, rho1 = 0.45, beta = 0)
ev_cur <- rowMeans(replicate(10, spec(noise_decay(150, 200, 3, 0.25, 0.45, 0, 1.0))))
ev_new <- rowMeans(replicate(10, spec(noise_decay(150, 200, best$K, best$alpha,
                                                  best$rho1, best$beta, 1.0))))
cat("\n=== spectra (%) ===\nPC       :", sprintf("%5d", 1:10), "\n")
cat("real     :", sprintf("%5.1f", 100 * ev_real[1:10]), "\n")
cat("current  :", sprintf("%5.1f", 100 * ev_cur[1:10]), "\n")
cat("decaying :", sprintf("%5.1f", 100 * ev_new[1:10]), "\n")
cat(sprintf("RMSE(PC1-20) vs real:  current %.3f   decaying %.3f\n",
    sqrt(mean((ev_cur - ev_real)^2)) * 100, sqrt(mean((ev_new - ev_real)^2)) * 100))

# --- Does it change power? ----------------------------------------------------
N <- 100; P <- 200; Q <- 50; RX <- 20; RY <- 10; NF <- 10; REPS <- 150
CFG <- list(current  = list(K = 3, alpha = 0.25, rho1 = 0.45, beta = 0),
            decaying = list(K = best$K, alpha = best$alpha, rho1 = best$rho1, beta = best$beta))
SS <- c(0, 2, 4, 8)
sl <- function(te, tr) scale(te, center = attr(tr, "scaled:center"), scale = attr(tr, "scaled:scale"))

cell <- function(rep_id, s, cn) {
  g <- CFG[[cn]]; set.seed(rep_id * 7919 + round(s * 100) + nchar(cn))
  Z <- rnorm(N)
  lx <- numeric(P); ly <- numeric(Q)
  lx[1:RX] <- rnorm(RX); ly[1:RY] <- rnorm(RY)
  lx <- lx / sqrt(sum(lx^2)); ly <- ly / sqrt(sum(ly^2))
  X <- sqrt(s) * outer(Z, lx) + noise_decay(N, P, g$K, g$alpha, g$rho1, g$beta, 1.0)
  Y <- sqrt(s) * outer(Z, ly) + noise_decay(N, Q, g$K, g$alpha, g$rho1, g$beta, 1.0)
  set.seed(rep_id); folds <- createFolds(seq_len(N), k = NF, list = TRUE)
  fc <- vapply(seq_len(NF), function(i) {
    te <- folds[[i]]; tr <- setdiff(seq_len(N), te)
    Xtr <- scale(X[tr, , drop = FALSE]); Ytr <- scale(Y[tr, , drop = FALSE])
    f <- tryCatch(spls(Xtr, Ytr, ncomp = 1, keepX = RX, keepY = RY,
                       mode = "regression", scale = TRUE), error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    cor(sl(X[te, , drop = FALSE], Xtr) %*% f$loadings$X[, 1],
        sl(Y[te, , drop = FALSE], Ytr) %*% f$loadings$Y[, 1], method = "spearman")[1, 1]
  }, numeric(1))
  data.frame(cfg = cn, s = s, rep = rep_id, stat = mean(fc, na.rm = TRUE))
}
jobs <- expand.grid(rep_id = seq_len(REPS), s = SS, cn = names(CFG), stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(k) cell(jobs$rep_id[k], jobs$s[k], jobs$cn[k]),
                mc.cores = max(1, detectCores() - 1))
out <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])

cat("\n=== power (vs each structure's own null 95% threshold) ===\n")
for (cn in names(CFG)) {
  a <- out[out$cfg == cn, ]; thr <- quantile(a$stat[a$s == 0], 0.95, na.rm = TRUE)
  pw <- vapply(SS[-1], function(s) mean(a$stat[a$s == s] > thr, na.rm = TRUE), numeric(1))
  cat(sprintf("%-9s null_sd=%.3f thr=%.3f | ", cn, sd(a$stat[a$s == 0], na.rm = TRUE), thr),
      paste(sprintf("s=%g: %.3f", SS[-1], pw), collapse = "  "), "\n")
}
saveRDS(list(fit = fit, out = out), file.path(SP, "spectrum_sensitivity.rds"))
