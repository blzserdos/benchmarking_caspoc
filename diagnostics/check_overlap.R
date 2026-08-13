#!/usr/bin/env Rscript
# =============================================================================
# check_overlap.R
# Sanity + regression checks for the overlapping factor structure.
#
#   1. REGRESSION: structure = NULL and the modular structure are byte-identical
#      to the pre-change behaviour (the overlap branch must be purely additive)
#   2. Unit marginal variance is exact; psi is always a valid variance
#   3. Sigma is redrawn per replicate but reproducible from the ambient RNG
#   4. pi_join is the sparsity dial: low pi_join gives unstructured features
#   5. Descriptive statistics against real breast.TCGA mRNA
#   6. The oracle identity s/(s + l'Sigma l) still centres where it should
#   7. SIGNAL_ALIGNMENT is inert under overlap
#
# Usage:  Rscript diagnostics/check_overlap.R
# =============================================================================
suppressMessages({ library(mixOmics) })
source("R/generate_data.R")
options(width = 130)

ok <- function(lab, pass, detail = "") {
  cat(sprintf("  [%s] %-52s %s\n", if (pass) "PASS" else "FAIL", lab, detail))
  invisible(pass)
}
N <- 150; P <- 200

cat("=== 1. regression: existing behaviour unchanged ===\n")

set.seed(11); a <- make_noise(N, P, structure = NULL)
set.seed(11); b <- matrix(rnorm(N * P), N, P)
ok("structure = NULL is plain independent noise", identical(a, b))

# the pre-change modular algorithm, reproduced inline
old_modular <- function(n, p, st) {
  m <- module_sizes(p, st$n_modules, st$alpha, st$frac)
  beta <- if (is.null(st$beta)) 0 else st$beta
  rk <- module_cors(length(m), st$rho, beta)
  E <- matrix(rnorm(n * p), nrow = n, ncol = p)
  pos <- 0L
  for (k in seq_along(m)) {
    idx <- (pos + 1L):(pos + m[k]); pos <- pos + m[k]
    E[, idx] <- sqrt(rk[k]) * rnorm(n) + sqrt(1 - rk[k]) * E[, idx]
  }
  E
}
set.seed(22); a <- make_noise(N, P, structure = DEFAULT_BLOCK_STRUCTURE)
set.seed(22); b <- old_modular(N, P, DEFAULT_BLOCK_STRUCTURE)
ok("modular structure is byte-identical", identical(a, b))

st_typed <- c(DEFAULT_BLOCK_STRUCTURE, list(type = "modules"))
set.seed(22); a2 <- make_noise(N, P, structure = st_typed)
ok("explicit type = 'modules' matches implicit", identical(a2, b))

set.seed(33); a <- make_noise(N, P, cor_within = 0.3)
ok("cor_within path untouched", is.matrix(a) && all(dim(a) == c(N, P)))

cat("\n=== 2. the invariants the signal design depends on ===\n")
OV <- DEFAULT_OVERLAP_STRUCTURE
L <- overlap_loadings(P, OV)
ok("communality hits its target exactly",
   max(abs(rowSums(L$A^2) - L$h2)) < 1e-12,
   sprintf("max err %.1e", max(abs(rowSums(L$A^2) - L$h2))))
ok("psi is always a valid variance (0 <= psi <= 1)",
   all(L$psi >= 0) && all(L$psi <= 1),
   sprintf("psi in [%.2f, %.2f]", min(L$psi), max(L$psi)))
ok("no clipping: h2 never pinned to the ceiling",
   mean(abs(L$h2 - H2_MAX) < 1e-9) < 0.02,
   sprintf("%.1f%% at ceiling", 100 * mean(abs(L$h2 - H2_MAX) < 1e-9)))
S <- L$A %*% t(L$A); diag(S) <- diag(S) + L$psi
ok("population marginal variance is exactly 1",
   max(abs(diag(S) - 1)) < 1e-12, sprintf("max err %.1e", max(abs(diag(S) - 1))))

set.seed(44)
V <- apply(make_noise(20000, P, structure = OV), 2, var)
ok("empirical marginal variance ~ 1 (n = 20000)",
   abs(mean(V) - 1) < 0.02, sprintf("mean %.3f, range [%.2f, %.2f]",
                                    mean(V), min(V), max(V)))

cat("\n=== 3. Sigma is an ensemble, driven by the ambient RNG ===\n")
set.seed(1); R1 <- overlap_loadings(P, OV)
set.seed(2); R2 <- overlap_loadings(P, OV)
ok("Sigma is redrawn per replicate", !identical(R1$A, R2$A))

set.seed(1); R1b <- overlap_loadings(P, OV)
ok("...but reproducible: same seed, same Sigma", identical(R1$A, R1b$A))

set.seed(9); LX <- overlap_loadings(P, OV); LY <- overlap_loadings(50, OV)
ok("X and Y blocks get different structures",
   !isTRUE(all.equal(LX$A[1:50, ], LY$A)))

cat("\n=== 4. pi_join is the sparsity dial ===\n")
cat("  pi_join  fac/feat  unstructured%  mean_psi  mean|c|   RMS   %neg  %|c|<.02\n")
for (pj in c(0.02, 0.05, 0.1, 0.2, 0.35, 0.5)) {
  s <- modifyList(OV, list(pi_join = pj))
  set.seed(100 + round(1000 * pj))
  st <- replicate(6, {
    Lp <- overlap_loadings(P, s)
    Sp <- Lp$A %*% t(Lp$A); diag(Sp) <- diag(Sp) + Lp$psi
    cc <- Sp[upper.tri(Sp)]
    c(nj = mean(rowSums(Lp$A != 0)), un = mean(Lp$h2 == 0), psi = mean(Lp$psi),
      mabs = mean(abs(cc)), rms = sqrt(mean(cc^2)), neg = mean(cc < 0),
      z = mean(abs(cc) < 0.02))
  })
  r <- rowMeans(st)
  cat(sprintf("   %.2f    %6.1f      %6.1f       %6.2f   %6.3f  %.3f  %5.1f    %5.1f\n",
              pj, r["nj"], 100 * r["un"], r["psi"], r["mabs"], r["rms"],
              100 * r["neg"], 100 * r["z"]))
}
cat("  (low pi_join restores a genuinely unstructured stratum; pi_join ~ 1/K\n")
cat("   is roughly one factor per feature, i.e. the disjoint model)\n")

cat("\n=== 5. descriptive statistics vs real breast.TCGA mRNA (n=150, p=200) ===\n")
stats_of <- function(M) {
  M <- scale(M)
  ev <- prcomp(M, center = FALSE, scale. = FALSE)$sdev^2; ev <- ev / sum(ev)
  o <- cor(M)[upper.tri(diag(ncol(M)))]
  list(ev = ev[1:20], mabs = mean(abs(o)), rms = sqrt(mean(o^2)),
       med = median(abs(o)), neg = mean(o < 0))
}
data("breast.TCGA")
R <- stats_of(as.matrix(breast.TCGA$data.train$mrna))
row_for <- function(lab, st) {
  set.seed(555)
  s <- replicate(10, stats_of(make_noise(N, P, structure = st)), simplify = FALSE)
  ev <- rowMeans(vapply(s, function(x) x$ev, numeric(20)))
  g <- function(f) mean(vapply(s, f, numeric(1)))
  cat(sprintf("  %-24s %.3f    %.3f  %.3f     %5.1f%%   %5.2f\n", lab,
              g(function(x) x$mabs), g(function(x) x$rms), g(function(x) x$med),
              100 * g(function(x) x$neg), sqrt(mean((ev - R$ev)^2)) * 100))
}
cat("                           mean|c|   RMS   median|c|   %neg   specRMSE\n")
cat(sprintf("  %-24s %.3f    %.3f  %.3f     %5.1f%%      --\n",
            "real mRNA", R$mabs, R$rms, R$med, 100 * R$neg))
row_for("disjoint (current)", DEFAULT_BLOCK_STRUCTURE)
row_for("overlap", OV)
cat("\n  These need only be in the right BALLPARK -- the goal is omics-like data\n")
cat("  for benchmarking CV methods, not a breast.TCGA replica.\n")

cat("\n=== 6. oracle identity: s / (s + l'Sigma l) at nominal s = 1 ===\n")
oracle_row <- function(lab, Sig, idx) {
  v <- vapply(1:400, function(i) {
    l <- rnorm(length(idx)); l <- l / sqrt(sum(l^2))
    as.numeric(t(l) %*% Sig[idx, idx] %*% l)
  }, numeric(1))
  oc <- 1 / (1 + v)
  cat(sprintf("  %-22s mean=%.3f  sd=%.3f  IQR=[%.2f, %.2f]   (target 0.500)\n",
              lab, mean(oc), sd(oc), quantile(oc, .25), quantile(oc, .75)))
}
m <- module_sizes(P, 8, 0.5, 1.0); rk <- module_cors(8, 0.8, 0.5)
Sd <- matrix(0, P, P); pos <- 0L
for (k in seq_along(m)) { ii <- (pos + 1L):(pos + m[k]); pos <- pos + m[k]
  Sd[ii, ii] <- rk[k] }
diag(Sd) <- 1
set.seed(66)
oracle_row("disjoint / spread", Sd,
           relevant_indices(P, 20, 1, DEFAULT_BLOCK_STRUCTURE, "spread"))
oracle_row("overlap", S, 1:20)

cat("\n=== 7. SIGNAL_ALIGNMENT is inert under overlap ===\n")
i1 <- relevant_indices(P, 20, 1, OV, "spread")
i2 <- relevant_indices(P, 20, 1, OV, "aligned")
ok("spread and aligned coincide", identical(i1, i2), sprintf("both = %d:%d", min(i1), max(i1)))

d1 <- generate_signal_data(50, P, 50, 1, 20, 10, 1, structure = OV,
                           signal_alignment = "spread", seed = 3)
d2 <- generate_signal_data(50, P, 50, 1, 20, 10, 1, structure = OV,
                           signal_alignment = "aligned", seed = 3)
ok("generators agree given the same seed", identical(d1$X, d2$X))

cat("\n=== 8. generators run end to end ===\n")
dn <- generate_null_data(100, 200, 50, structure = OV, seed = 1)
ds <- generate_signal_data(100, 200, 50, 1, 20, 10, 4, structure = OV, seed = 1)
ok("generate_null_data", all(dim(dn$X) == c(100, 200)) && all(dim(dn$Y) == c(100, 50)))
ok("generate_signal_data", all(dim(ds$X) == c(100, 200)) && !any(is.na(ds$X)))
ok("no cross-block leakage in the null",
   abs(mean(cor(dn$X, dn$Y))) < 0.02, sprintf("mean cross-cor %.4f", mean(cor(dn$X, dn$Y))))
