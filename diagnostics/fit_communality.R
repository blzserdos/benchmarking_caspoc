#!/usr/bin/env Rscript
# =============================================================================
# fit_communality.R
# Estimates the COMMUNALITY DISTRIBUTION of a real block, as a ballpark check
# on h2_mean in the overlapping factor structure (see R/generate_data.R).
#
# Communality h2_j = share of feature j's variance explained by the common
# factors. Under a K-factor model extracted by PCA, h2_j = sum over the top K
# components of (eigenvector_jk * sqrt(eigenvalue_k))^2, applied to the
# correlation matrix.
#
# THE BIAS THAT MATTERS: with n = 150 and p = 200 the sample correlation matrix
# is rank-deficient and noisy, so the top K PCs soak up a great deal of pure
# sampling noise and h2 is inflated. This script therefore runs the identical
# procedure on INDEPENDENT data at matched n, p, K to get the noise floor, and
# reports a parallel-analysis-style adjusted estimate:
#
#     h2_adj = (h2_real - h2_floor) / (1 - h2_floor)
#
# i.e. the share of variance explained beyond what K components would explain
# on data with no structure at all. Rescaled so an unstructured feature maps to
# 0 and a fully-explained one to 1.
#
# Usage:  Rscript diagnostics/fit_communality.R [K]
# =============================================================================
suppressMessages({ library(mixOmics) })
options(width = 130)

args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.integer(args[1]) else 40

communality <- function(M, K) {
  M <- scale(M)
  C <- cor(M)
  e <- eigen(C, symmetric = TRUE)
  lam <- pmax(0, e$values[seq_len(K)])
  L <- e$vectors[, seq_len(K), drop = FALSE] * rep(sqrt(lam), each = nrow(C))
  pmin(1, rowSums(L^2))
}

data("breast.TCGA")
blocks <- breast.TCGA$data.train[c("mrna", "mirna", "protein")]

cat(sprintf("K = %d factors\n\n", K))
fits <- list()

for (nm in names(blocks)) {
  M <- as.matrix(blocks[[nm]])
  n <- nrow(M); p <- ncol(M)
  if (K >= min(n - 1, p)) {
    cat(sprintf("%-8s SKIPPED (K=%d not < min(n-1,p)=%d)\n", nm, K, min(n - 1, p)))
    next
  }

  h2_real <- communality(M, K)

  # noise floor: same n, p, K, but no structure whatsoever
  set.seed(1)
  floor_reps <- replicate(10, mean(communality(matrix(rnorm(n * p), n, p), K)))
  h2_floor <- mean(floor_reps)

  h2_adj <- pmax(0, pmin(0.99, (h2_real - h2_floor) / (1 - h2_floor)))

  q <- quantile(h2_adj, c(0.05, 0.25, 0.5, 0.75, 0.95))
  cat(sprintf("%-8s n=%3d p=%3d | raw h2 mean=%.2f | noise floor=%.2f | ADJUSTED mean=%.2f sd=%.2f\n",
              nm, n, p, mean(h2_real), h2_floor, mean(h2_adj), sd(h2_adj)))
  cat(sprintf("         adjusted h2 quantiles (5/25/50/75/95): %s\n",
              paste(sprintf("%.2f", q), collapse = "  ")))

  fits[[nm]] <- h2_adj
}

# --- Summarise on the same scale the generator uses --------------------------
H2_MIN <- 0.02; H2_MAX <- 0.95
cat(sprintf("\n=== suggested structure parameters (Beta on [%.2f, %.2f]) ===\n",
            H2_MIN, H2_MAX))
for (nm in names(fits)) {
  h <- pmin(H2_MAX - 1e-6, pmax(H2_MIN + 1e-6, fits[[nm]]))
  u <- (h - H2_MIN) / (H2_MAX - H2_MIN)
  m <- mean(u); v <- var(u)
  conc <- max(0.5, m * (1 - m) / v - 1)   # method of moments for Beta
  cat(sprintf("  %-8s  h2_mean = %.2f   (spread implies concentration ~%.0f)\n",
              nm, H2_MIN + m * (H2_MAX - H2_MIN), conc))
}

cat("\nThis is a BALLPARK CHECK, not a fitting step. It answers 'is h2_mean set to\n")
cat("roughly the right order of magnitude for omics data?' -- nothing finer. The\n")
cat("estimate is PCA extraction with a crude noise correction, not a maximum\n")
cat("likelihood factor fit, and the simulation is not meant to reproduce any\n")
cat("particular dataset. Only h2_mean is a knob in the generator; the spread is\n")
cat("fixed internally, since values from 5 to 50 move mean|c| by 0.002.\n")
