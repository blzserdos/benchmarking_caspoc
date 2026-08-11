#!/usr/bin/env Rscript
# =============================================================================
# fit_structure.R
# Fits a modular within-block correlation structure to breast.TCGA.
#
# Structure: K disjoint modules of correlated features (co-expression modules /
# co-occurring taxa), sizes decaying as k^-alpha, within-module correlation rho,
# covering frac*p of the features; the remainder is independent.
#
# Generated as a factor model rather than via mvrnorm -- O(n*p) instead of an
# eigendecomposition of a p x p matrix:
#     feature j in module k:  sqrt(rho) * F_k + sqrt(1 - rho) * eps_j
#
# Fit by matching SAMPLE observables (PC1 / PC1-5 / PC1-10 variance fractions
# and mean |correlation|) computed identically on real and simulated blocks at
# the same n and p, so estimation noise cancels out of the comparison.
# =============================================================================

suppressMessages({ library(mixOmics); library(parallel) })

OUT_FILE <- file.path(dirname(sub("^--file=", "", grep("^--file=",
              commandArgs(FALSE), value = TRUE)[1])), "fit_structure.rds")

# --- Observables --------------------------------------------------------------
observables <- function(M) {
  M <- scale(M)
  ev <- prcomp(M, center = FALSE, scale. = FALSE)$sdev^2
  ev <- ev / sum(ev)
  cm <- cor(M)
  c(pc1 = ev[1], pc5 = sum(ev[1:5]), pc10 = sum(ev[1:10]),
    mabs = mean(abs(cm[upper.tri(cm)])))
}

# --- Modular sampler ----------------------------------------------------------
module_sizes <- function(p, K, alpha, frac) {
  w <- (seq_len(K))^(-alpha)
  m <- pmax(2, round(w / sum(w) * frac * p))
  while (sum(m) > p) m[which.max(m)] <- m[which.max(m)] - 1
  m
}

sample_modular <- function(n, p, K, alpha, rho, frac) {
  m <- module_sizes(p, K, alpha, frac)
  X <- matrix(rnorm(n * p), n, p)
  pos <- 0
  for (k in seq_along(m)) {
    idx <- (pos + 1):(pos + m[k]); pos <- pos + m[k]
    Fk <- rnorm(n)
    X[, idx] <- sqrt(rho) * Fk + sqrt(1 - rho) * X[, idx]
  }
  X
}

# --- Real targets -------------------------------------------------------------
data("breast.TCGA")
blocks <- breast.TCGA$data.train[c("mirna", "mrna", "protein")]
real <- t(vapply(blocks, function(b) observables(as.matrix(b)), numeric(4)))
cat("=== Real blocks (n = 150) ===\n"); print(round(real, 3))

# --- Grid search --------------------------------------------------------------
grid <- expand.grid(K = c(3, 5, 10, 20), alpha = c(0, 0.5, 1, 1.5),
                    rho = c(0.2, 0.3, 0.5, 0.7, 0.9), frac = c(0.4, 0.7, 1.0))
N_REP <- 6
TARGET_BLOCK <- "mrna"          # p = 200, matches the simulation design
tgt <- real[TARGET_BLOCK, ]
n_t <- nrow(blocks[[TARGET_BLOCK]]); p_t <- ncol(blocks[[TARGET_BLOCK]])

res <- mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  o <- vapply(seq_len(N_REP), function(r) {
    set.seed(r * 101 + i)
    observables(sample_modular(n_t, p_t, g$K, g$alpha, g$rho, g$frac))
  }, numeric(4))
  o <- rowMeans(o)
  # Relative error on the three spectrum observables (mean |cor| is the most
  # noise-contaminated, so it is reported but down-weighted in the fit).
  err <- mean(abs(o[1:3] - tgt[1:3]) / tgt[1:3]) + 0.25 * abs(o[4] - tgt[4]) / tgt[4]
  cbind(g, pc1 = o[1], pc5 = o[2], pc10 = o[3], mabs = o[4], err = err)
}, mc.cores = max(1, detectCores() - 1))

out <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])
out <- out[order(out$err), ]
cat("\n=== Best 10 fits to", TARGET_BLOCK, "===\n")
print(round(head(out, 10), 3), row.names = FALSE)
cat("\ntarget:", paste(sprintf("%s=%.3f", names(tgt), tgt), collapse = "  "), "\n")

saveRDS(list(real = real, grid = out, target_block = TARGET_BLOCK), OUT_FILE)
message("Saved -> ", OUT_FILE)
