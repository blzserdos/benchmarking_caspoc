#!/usr/bin/env Rscript
# =============================================================================
# calibrate_s.R
# What signal_strength makes sim_signal look like real omics data?
#
# Applies the repo's own run_naive_cv() to (a) real block pairs from
# breast.TCGA and (b) simulated data at matched n/p/q across a grid of s,
# then reads off the s whose CV statistic matches each real pair.
#
# naive_cv is optimistically biased, but the bias applies identically to real
# and simulated data, so the matched s is still meaningful. Real blocks are
# subsampled to n = 100 to match the simulation design.
# =============================================================================

suppressMessages({
  library(MASS); library(mixOmics); library(caret); library(dplyr)
  library(tibble); library(parallel)
})
source("R/generate_data.R")

# Load the repo's real CV code, minus the caspoc dependency (not installed
# locally and not needed for naive_cv).
src <- readLines("R/cv_approaches.R")
eval(parse(text = paste(src[!grepl("^library\\(caspoc\\)", src)], collapse = "\n")),
     envir = globalenv())

N_SUB     <- 100          # match the simulation's sample size
N_REPS    <- 15
STRENGTHS <- c(0, 0.5, 1, 2, 4, 8, 16, 32)
KEEPX     <- c(5, 10, 20, 50)
KEEPY     <- c(5, 10, 20)
OUT_FILE  <- file.path(dirname(sub("^--file=", "", grep("^--file=",
              commandArgs(FALSE), value = TRUE)[1])), "calib_s.rds")

data("breast.TCGA")
blocks <- breast.TCGA$data.train[c("mirna", "mrna", "protein")]
pairs <- list(
  c("mrna", "protein"),
  c("mirna", "protein"),
  c("mirna", "mrna")
)

stat_of <- function(X, Y, seed) {
  r <- tryCatch(run_naive_cv(X, Y, ncomp = 1, num_folds = 10,
                             keepX_options = KEEPX, keepY_options = KEEPY,
                             seed = seed),
                error = function(e) NULL)
  if (is.null(r)) NA_real_ else r$observed_stat
}

# --- Within-block structure, for the realism caveat ---------------------------
block_structure <- do.call(rbind, lapply(names(blocks), function(b) {
  M <- scale(as.matrix(blocks[[b]]))
  ev <- prcomp(M, center = FALSE, scale. = FALSE)$sdev^2
  cm <- cor(M); offdiag <- cm[upper.tri(cm)]
  data.frame(block = b, p = ncol(M),
             pc1_var_frac = ev[1] / sum(ev),
             pc1_5_var_frac = sum(ev[1:5]) / sum(ev),
             mean_abs_cor = mean(abs(offdiag)))
}))

# --- Real pairs ---------------------------------------------------------------
real_jobs <- expand.grid(rep_id = seq_len(N_REPS), pair = seq_along(pairs))
real_res <- mclapply(seq_len(nrow(real_jobs)), function(k) {
  pr <- pairs[[real_jobs$pair[k]]]; rid <- real_jobs$rep_id[k]
  X <- as.matrix(blocks[[pr[1]]]); Y <- as.matrix(blocks[[pr[2]]])
  set.seed(1000 + rid)
  idx <- sample(nrow(X), N_SUB)
  data.frame(source = "real", label = paste(pr, collapse = " ~ "),
             p = ncol(X), q = ncol(Y), signal_strength = NA_real_,
             rep = rid, stat = stat_of(X[idx, ], Y[idx, ], seed = rid))
}, mc.cores = max(1, detectCores() - 1))

# --- Simulated, at each real pair's dimensions --------------------------------
sim_jobs <- expand.grid(rep_id = seq_len(N_REPS), s = STRENGTHS,
                        pair = seq_along(pairs))
sim_res <- mclapply(seq_len(nrow(sim_jobs)), function(k) {
  pr <- pairs[[sim_jobs$pair[k]]]; rid <- sim_jobs$rep_id[k]; s <- sim_jobs$s[k]
  p <- ncol(blocks[[pr[1]]]); q <- ncol(blocks[[pr[2]]])
  d <- generate_signal_data(n = N_SUB, p = p, q = q, n_comp_true = 1,
                            n_relevant_x = 20, n_relevant_y = 10,
                            signal_strength = s,
                            seed = rid * 7919 + round(s * 1000) + p)
  data.frame(source = "sim", label = paste(pr, collapse = " ~ "),
             p = p, q = q, signal_strength = s,
             rep = rid, stat = stat_of(d$X, d$Y, seed = rid))
}, mc.cores = max(1, detectCores() - 1))

keep <- function(l) do.call(rbind, l[vapply(l, is.data.frame, logical(1))])
saveRDS(list(real = keep(real_res), sim = keep(sim_res),
             structure = block_structure), OUT_FILE)
message("Saved -> ", OUT_FILE)
