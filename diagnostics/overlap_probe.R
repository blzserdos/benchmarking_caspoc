#!/usr/bin/env Rscript
# =============================================================================
# overlap_probe.R
# Can an OVERLAPPING factor model match the eigenvalue spectrum AND mean|cor|
# simultaneously, where disjoint blocks cannot?
#
# Model:  X_j = sum_k a_jk F_k + sqrt(psi_j) eps_j,   F_k, eps_j ~ N(0,1) iid
#         a_jk = sqrt(w_k) * s_jk,  s_jk in {0,+1,-1}
#         feature j joins factor k with probability pi;  sign is random
#         psi_j = 1 - sum_k a_jk^2   (so every feature keeps unit variance)
#
# Because a feature can join SEVERAL factors, cor(i,j) = sum_k a_ik a_jk takes a
# continuum of values instead of the binary {rho_k, 0} of disjoint blocks.
# Random signs also produce genuine negative correlations, which the
# block-diagonal structure cannot.
# =============================================================================
suppressMessages({ library(mixOmics); library(parallel) })
options(width = 150)
SP <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))

sample_overlap <- function(n, p, K, w1, gamma, pi_join, cap = 0.95) {
  w <- w1 * (seq_len(K))^(-gamma)
  A <- matrix(0, p, K)
  for (k in seq_len(K)) {
    idx <- which(runif(p) < pi_join)
    if (length(idx)) A[idx, k] <- sqrt(w[k]) * sample(c(-1, 1), length(idx), TRUE)
  }
  ss <- rowSums(A^2); over <- ss > cap
  if (any(over)) A[over, ] <- A[over, ] / sqrt(ss[over] / cap)
  psi <- pmax(1e-6, 1 - rowSums(A^2))
  matrix(rnorm(n * K), n, K) %*% t(A) +
    sweep(matrix(rnorm(n * p), n, p), 2, sqrt(psi), "*")
}

stats_of <- function(M) {
  M <- scale(M)
  ev <- prcomp(M, center = FALSE, scale. = FALSE)$sdev^2; ev <- ev / sum(ev)
  o <- cor(M)[upper.tri(diag(ncol(M)))]
  list(ev = ev[1:20], mabs = mean(abs(o)), rms = sqrt(mean(o^2)),
       neg = mean(o < 0), med = median(abs(o)))
}

data("breast.TCGA")
R <- stats_of(as.matrix(breast.TCGA$data.train$mrna))
cat(sprintf("TARGET (real mRNA): mean|c|=%.3f  RMS=%.3f  median|c|=%.3f  %%neg=%.1f%%\n\n",
            R$mabs, R$rms, R$med, 100 * R$neg))

grid <- expand.grid(K = c(5, 10, 20, 40), gamma = c(0.5, 1, 1.5),
                    w1 = c(0.3, 0.5, 0.7), pi_join = c(0.1, 0.2, 0.3, 0.5))
res <- mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; set.seed(i * 17 + 3)
  s <- replicate(8, stats_of(sample_overlap(150, 200, g$K, g$w1, g$gamma, g$pi_join)),
                 simplify = FALSE)
  ev <- rowMeans(vapply(s, function(x) x$ev, numeric(20)))
  mabs <- mean(vapply(s, function(x) x$mabs, numeric(1)))
  spec_rmse <- sqrt(mean((ev - R$ev)^2)) * 100
  cbind(g, spec_rmse = spec_rmse, mabs = mabs,
        mabs_err = abs(mabs - R$mabs),
        med = mean(vapply(s, function(x) x$med, numeric(1))),
        neg = mean(vapply(s, function(x) x$neg, numeric(1))),
        # joint objective: spectrum RMSE (% pts) + mean|c| error scaled to match
        obj = spec_rmse + 100 * abs(mabs - R$mabs) / 5)
}, mc.cores = max(1, detectCores() - 1))
out <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])

cat("=== best by JOINT objective (spectrum + mean|cor|) ===\n")
print(round(head(out[order(out$obj), ], 6), 3), row.names = FALSE)
cat("\n=== best by spectrum alone ===\n")
print(round(head(out[order(out$spec_rmse), ], 3), 3), row.names = FALSE)

b <- out[order(out$obj), ][1, ]
set.seed(99)
s <- replicate(15, stats_of(sample_overlap(150, 200, b$K, b$w1, b$gamma, b$pi_join)),
               simplify = FALSE)
ev <- rowMeans(vapply(s, function(x) x$ev, numeric(20)))
cat("\n=== best overlapping config vs real, and vs the current disjoint structure ===\n")
cat("PC          :", sprintf("%5d", 1:10), "\n")
cat("real   (%)  :", sprintf("%5.1f", 100 * R$ev[1:10]), "\n")
cat("overlap(%)  :", sprintf("%5.1f", 100 * ev[1:10]), "\n")
cat("disjoint(%) :", sprintf("%5.1f", c(19.2, 10.0, 6.6, 5.2, 4.1, 3.6, 3.1, 2.6, 1.2, 1.1)), "\n\n")
cat(sprintf("             mean|c|  RMS   median|c|  %%neg   spectrumRMSE\n"))
cat(sprintf("real         %.3f    %.3f  %.3f      %.1f%%    --\n", R$mabs, R$rms, R$med, 100*R$neg))
cat(sprintf("overlapping  %.3f    %.3f  %.3f      %.1f%%    %.2f\n",
    mean(vapply(s, function(x) x$mabs, numeric(1))),
    mean(vapply(s, function(x) x$rms, numeric(1))),
    mean(vapply(s, function(x) x$med, numeric(1))),
    100 * mean(vapply(s, function(x) x$neg, numeric(1))),
    sqrt(mean((ev - R$ev)^2)) * 100))
cat(sprintf("disjoint     0.134    0.234  0.066      ~50%%*   0.83   (*noise only, no structural neg)\n"))
saveRDS(out, file.path(SP, "overlap_probe.rds"))
