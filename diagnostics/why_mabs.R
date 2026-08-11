#!/usr/bin/env Rscript
# Why does mean|cor| undershoot (0.135 vs 0.200) when the spectrum matches?
#
# Key identity: ||C||_F^2 = sum_k lambda_k^2 = p + sum_{i!=j} c_ij^2.
# So matching the eigenvalue spectrum pins the SUM OF SQUARES of the
# off-diagonal correlations -- i.e. their RMS -- but says nothing about how that
# mass is distributed across pairs. mean|c| depends on the distribution shape.
suppressMessages(library(mixOmics))
source("R/generate_data.R"); source("cluster/config.R")
options(width = 150)

offstats <- function(C) {
  o <- C[upper.tri(C)]
  c(mean_abs = mean(abs(o)), rms = sqrt(mean(o^2)),
    ratio = mean(abs(o)) / sqrt(mean(o^2)),
    q50 = quantile(abs(o), .50), q90 = quantile(abs(o), .90),
    q99 = quantile(abs(o), .99), frac_gt_.3 = mean(abs(o) > .3))
}

data("breast.TCGA")
Mr <- scale(as.matrix(breast.TCGA$data.train$mrna))
Cr <- cor(Mr); er <- prcomp(Mr, center = FALSE, scale. = FALSE)$sdev^2

sim1 <- replicate(10, { X <- generate_signal_data(150, 200, 50, 1, 20, 10, 0,
                          structure = BLOCK_STRUCTURE, seed = NULL)$X
                        offstats(cor(scale(X))) })
old  <- list(n_modules = 3, alpha = 0.25, rho = 0.45, beta = 0, frac = 1.0)
sim0 <- replicate(10, { X <- generate_signal_data(150, 200, 50, 1, 20, 10, 0,
                          structure = old, seed = NULL)$X
                        offstats(cor(scale(X))) })

cat("=== off-diagonal correlation distribution (sample, n=150, p=200) ===\n")
tab <- rbind(real = offstats(Cr), `sim: 8 decaying modules` = rowMeans(sim1),
             `sim: 3 equal modules (old)` = rowMeans(sim0))
print(round(tab, 3))

cat("\n=== Frobenius check: does the spectrum really pin the sum of squares? ===\n")
fro <- function(e) sum(e^2)
e1 <- rowMeans(replicate(10, prcomp(scale(generate_signal_data(150,200,50,1,20,10,0,
        structure = BLOCK_STRUCTURE, seed=NULL)$X), center=FALSE, scale.=FALSE)$sdev^2))
e0 <- rowMeans(replicate(10, prcomp(scale(generate_signal_data(150,200,50,1,20,10,0,
        structure = old, seed=NULL)$X), center=FALSE, scale.=FALSE)$sdev^2))
cat(sprintf("  sum(lambda^2):  real %.0f   8-decaying %.0f   3-equal %.0f\n",
            fro(er), fro(e1), fro(e0)))
cat(sprintf("  implied RMS(c): real %.3f   8-decaying %.3f   3-equal %.3f\n",
    sqrt((fro(er)-200)/(200*199)), sqrt((fro(e1)-200)/(200*199)),
    sqrt((fro(e0)-200)/(200*199))))

cat("\n=== where the correlation mass sits ===\n")
m <- module_sizes(200, BLOCK_STRUCTURE$n_modules, BLOCK_STRUCTURE$alpha, BLOCK_STRUCTURE$frac)
rk <- module_cors(length(m), BLOCK_STRUCTURE$rho, BLOCK_STRUCTURE$beta)
wp <- sum(choose(m, 2)); tp <- choose(200, 2)
cat(sprintf("  8 decaying: %d of %d pairs within a module (%.1f%%), rho = %s\n",
            wp, tp, 100 * wp / tp, paste(sprintf("%.2f", rk), collapse = " ")))
cat(sprintf("     -> POPULATION mean|c| = %.3f ; plus n=150 noise on the %.1f%% true zeros",
            sum(choose(m, 2) * rk) / tp, 100 * (1 - wp / tp)))
cat(sprintf(" (~%.3f each) = %.3f predicted sample mean|c|\n",
            sqrt(2 / pi) / sqrt(149),
            sum(choose(m, 2) * rk) / tp + (1 - wp / tp) * sqrt(2 / pi) / sqrt(149)))
m0 <- module_sizes(200, 3, 0.25, 1.0); wp0 <- sum(choose(m0, 2))
cat(sprintf("  3 equal   : %d of %d pairs (%.1f%%) at rho = 0.45\n", wp0, tp, 100 * wp0 / tp))
cat(sprintf("     -> POPULATION mean|c| = %.3f ; predicted sample mean|c| = %.3f\n",
            0.45 * wp0 / tp, 0.45 * wp0 / tp + (1 - wp0 / tp) * sqrt(2 / pi) / sqrt(149)))
