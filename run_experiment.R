#!/usr/bin/env Rscript
# =============================================================================
# run_experiment.R
# Pooled-null benchmark of the four CV approaches. Runs on a laptop.
#
#   Rscript run_experiment.R                 # run (resumable)
#   Rscript run_experiment.R --analyse       # re-print results, no computation
#
# THE DESIGN
#
# One factorial sweep, one CV run per cell, no permutations:
#
#   4 approaches x N_DATASETS datasets x (sim_null + SIGNAL_STRENGTHS)
#
# Two readouts come out of the same runs:
#
#   BIAS    mean statistic on sim_null. Cov(X, Y) = 0 there by construction, so
#           whatever a method reports IS its bias. Nothing else needed.
#
#   POWER   each method's own 95th-percentile null statistic is its critical
#           value c_M; power is the share of signal runs that clear it.
#
# WHY NO PERMUTATIONS. A permutation test manufactures null data by shuffling Y.
# We already simulate null data -- that is what sim_null is -- so the null
# distribution is available directly. That removes the 101x multiplier and takes
# the job count from 283,200 (cluster/, ~7,500 CPU-hours) to 3,200. Measured at
# ~18 s/job: ~16 CPU-hours, about 2 h wall clock on 9 cores.
# See cluster/ for the per-dataset permutation design, which answers
# the narrower practitioner question ("what power does ONE analyst get from
# method M plus a permutation test?") and is what a CASPOC user would actually
# run.
#
# WHY c_M IS PER METHOD. Naive CV's null statistics are inflated, so its
# threshold comes out high; CASPOC's comes out low. Each method is judged
# against its own null, exactly as a permutation test would do per dataset. A
# method therefore cannot win by inflating: inflation raises its threshold by
# the same amount. This is what makes power comparable across methods whose
# statistics live on different scales.
#
# WHAT THIS GIVES UP. The pooled threshold is marginal over Sigma draws (Sigma
# is redrawn every iteration), where a permutation test would condition on the
# dataset at hand. The null distribution is therefore slightly wider and power
# comes out mildly conservative -- equally so for every method, so the
# comparison is unaffected.
# =============================================================================

suppressMessages({
  library(dplyr)
  library(parallel)
})

source("R/generate_data.R")
source("R/cv_approaches.R")
source("cluster/config.R")   # APPROACHES, CV_CONFIG, HP_GRID, N_DATASETS,
                             # SIGNAL_STRENGTHS, BLOCK_STRUCTURE, datasets

# Optional overrides, for a quick end-to-end check before committing a machine
# to the full run. Results land in a separate file so a smoke test cannot
# contaminate the real one:
#
#   N_DATASETS=6 STRENGTHS=4,12 Rscript run_experiment.R
#
SMOKE <- nzchar(Sys.getenv("N_DATASETS")) || nzchar(Sys.getenv("STRENGTHS"))
if (nzchar(Sys.getenv("N_DATASETS"))) {
  N_DATASETS <- as.integer(Sys.getenv("N_DATASETS"))
}
if (nzchar(Sys.getenv("STRENGTHS"))) {
  SIGNAL_STRENGTHS <- as.numeric(strsplit(Sys.getenv("STRENGTHS"), ",")[[1]])
}
if (SMOKE) {
  message(sprintf("SMOKE TEST: %d datasets, signal strengths %s",
                  N_DATASETS, paste(SIGNAL_STRENGTHS, collapse = ", ")))
}

RESULTS_DIR <- "results"
RAW_FILE    <- file.path(RESULTS_DIR,
                         if (SMOKE) "smoke_experiment_raw.rds" else "experiment_raw.rds")
ALPHA       <- 0.05
# Under SLURM, detectCores() reports every core on the NODE, not the cores this
# job was allocated -- using it would oversubscribe a shared node badly. Honour
# the allocation when there is one, fall back to the local machine otherwise.
N_CORES     <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (is.na(N_CORES) || N_CORES < 1L) N_CORES <- max(1L, detectCores() - 1L)
BATCH       <- 200L          # checkpoint interval; also the progress granularity

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)


# =============================================================================
# Job grid
# =============================================================================
#
# One row per CV run. sim_signal carries a signal_strength, sim_null does not.
#
# Note the seeding: the generator seed is the iteration, so at a given iteration
# every signal strength shares the same Z, the same true loadings and the same
# noise draw -- only the signal scaling differs. The power curve is therefore
# PAIRED across strengths, which removes a large amount of between-draw noise
# from its shape.

build_grid <- function() {
  g <- rbind(
    data.frame(dataset = "sim_null",
               expand.grid(iteration       = seq_len(N_DATASETS),
                           approach        = APPROACHES,
                           signal_strength = NA_real_,
                           stringsAsFactors = FALSE),
               stringsAsFactors = FALSE),
    data.frame(dataset = "sim_signal",
               expand.grid(iteration       = seq_len(N_DATASETS),
                           approach        = APPROACHES,
                           signal_strength = SIGNAL_STRENGTHS,
                           stringsAsFactors = FALSE),
               stringsAsFactors = FALSE)
  )
  # Shuffle so each parallel batch gets a mix of fast (naive_cv) and slow
  # (nested_cv) jobs rather than a run of one kind.
  set.seed(42)
  g <- g[sample(nrow(g)), ]
  rownames(g) <- NULL
  g
}

job_key <- function(df) {
  paste(df$dataset, df$iteration, df$approach, df$signal_strength, sep = "|")
}


# =============================================================================
# One CV run
# =============================================================================

run_one <- function(job) {
  data <- datasets[[job$dataset]]$generator(
    seed            = job$iteration,
    signal_strength = job$signal_strength
  )

  # CASPOC() reports progress with cat()/print(); capture.output() keeps 3,200
  # runs of it out of the log without hiding real errors.
  res <- tryCatch({
    invisible(capture.output(
      out <- suppressMessages(run_cv_approach(
        approach        = job$approach,
        X               = data$X,
        Y               = data$Y,
        ncomp           = CV_CONFIG$ncomp,
        num_folds       = CV_CONFIG$num_folds,
        num_repeats     = CV_CONFIG$num_repeats,
        num_folds_inner = CV_CONFIG$num_folds_inner,
        keepX_options   = HP_GRID$keepX_options,
        keepY_options   = HP_GRID$keepY_options,
        seed            = job$iteration
      ))
    ))
    out
  }, error = function(e) NULL)

  data.frame(
    dataset         = job$dataset,
    iteration       = job$iteration,
    approach        = job$approach,
    signal_strength = job$signal_strength,
    observed_stat   = if (!is.null(res)) res$observed_stat  else NA_real_,
    selected_keepX  = if (!is.null(res)) res$selected_keepX else NA_real_,
    selected_keepY  = if (!is.null(res)) res$selected_keepY else NA_real_,
    runtime_sec     = if (!is.null(res)) res$runtime        else NA_real_,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Analysis: bias, threshold, power
# =============================================================================

analyse <- function(raw) {

  null_runs   <- raw %>% dplyr::filter(dataset == "sim_null")
  signal_runs <- raw %>% dplyr::filter(dataset == "sim_signal")

  # --- BIAS: truth is 0 on sim_null, so the reported statistic IS the bias ---
  bias <- null_runs %>%
    group_by(approach) %>%
    summarise(
      n_null       = sum(!is.na(observed_stat)),
      bias         = mean(observed_stat, na.rm = TRUE),
      sd_null      = sd(observed_stat, na.rm = TRUE),
      # type = 1 keeps the threshold an actually observed value, so the nominal
      # level is not inflated by interpolation between order statistics.
      crit_value   = quantile(observed_stat, 1 - ALPHA, na.rm = TRUE, type = 1),
      mean_runtime = mean(runtime_sec, na.rm = TRUE),
      .groups = "drop"
    )

  # --- POWER: share of signal runs clearing that method's own threshold ---
  power <- signal_runs %>%
    left_join(bias %>% dplyr::select(approach, crit_value), by = "approach") %>%
    group_by(approach, signal_strength) %>%
    summarise(
      n_signal  = sum(!is.na(observed_stat)),
      mean_stat = mean(observed_stat, na.rm = TRUE),
      power     = mean(observed_stat > crit_value, na.rm = TRUE),
      .groups   = "drop"
    )

  list(bias = bias, power = power)
}

report <- function(a) {
  cat("\n=====================================================================\n")
  cat("BIAS AND CALIBRATION  (sim_null: Cov(X,Y) = 0, so truth = 0.000)\n")
  cat("=====================================================================\n\n")
  cat(sprintf("%-14s %6s %10s %10s %12s %10s\n",
              "Approach", "n", "Bias", "SD", "Crit value", "Time(s)"))
  cat(strrep("-", 68), "\n")
  for (i in seq_len(nrow(a$bias))) with(a$bias[i, ],
    cat(sprintf("%-14s %6d %10.3f %10.3f %12.3f %10.1f\n",
                approach, n_null, bias, sd_null, crit_value, mean_runtime)))
  cat("\n  Bias       = mean statistic where the true association is zero.\n")
  cat("  Crit value = that method's 95th-percentile null statistic, used as its\n")
  cat("               own alpha = 0.05 threshold below. Higher bias -> higher\n")
  cat("               threshold, so inflation buys a method nothing.\n")

  cat("\n=====================================================================\n")
  cat("POWER  (share of signal datasets clearing that method's threshold)\n")
  cat("=====================================================================\n\n")
  strengths <- sort(unique(a$power$signal_strength))
  cat(sprintf("%-14s", "Approach"), sprintf("%7g", strengths), "\n")
  cat(sprintf("%-14s", "  s ="),    strrep(" ", 0), strrep("-", 8 * length(strengths)), "\n")
  for (ap in a$bias$approach) {
    row <- a$power %>% dplyr::filter(approach == ap) %>% arrange(signal_strength)
    cat(sprintf("%-14s", ap), sprintf("%7.2f", row$power), "\n")
  }
  cat("\n  Every method is held to a 0.05 false-positive rate by construction,\n")
  cat("  so these are directly comparable. 0.05 = no better than chance.\n\n")
}


# =============================================================================
# Main
# =============================================================================

if ("--analyse" %in% commandArgs(TRUE)) {
  if (!file.exists(RAW_FILE)) stop("No results at ", RAW_FILE, " -- run without --analyse first.")
  a <- analyse(readRDS(RAW_FILE))
  report(a)
  quit(save = "no")
}

grid <- build_grid()

# Resume: skip anything already on disk. An overnight laptop run should survive
# being interrupted.
done <- if (file.exists(RAW_FILE)) readRDS(RAW_FILE) else NULL
if (!is.null(done)) {
  todo <- grid[!job_key(grid) %in% job_key(done), ]
  message(sprintf("Resuming: %d of %d jobs already complete.", nrow(done), nrow(grid)))
} else {
  todo <- grid
}

message(sprintf("Jobs to run: %d | cores: %d | checkpoint every %d",
                nrow(todo), N_CORES, BATCH))
if (nrow(todo) > 0) {
  message(sprintf("Rough estimate: %.1f CPU-hours -> %.1f h wall clock.",
                  nrow(todo) * 18.3 / 3600, nrow(todo) * 18.3 / 3600 / N_CORES))
}

t0 <- proc.time()["elapsed"]
batches <- split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / BATCH))

for (b in seq_along(batches)) {
  idx <- batches[[b]]
  out <- mclapply(idx, function(i) run_one(todo[i, ]), mc.cores = N_CORES)
  out <- out[vapply(out, is.data.frame, logical(1))]

  done <- rbind(done, do.call(rbind, out))
  saveRDS(done, RAW_FILE)

  el   <- proc.time()["elapsed"] - t0
  frac <- max(idx) / nrow(todo)
  message(sprintf("  batch %d/%d | %d/%d jobs | %.1f min elapsed | ~%.1f min left",
                  b, length(batches), max(idx), nrow(todo),
                  el / 60, (el / frac - el) / 60))
}

message(sprintf("\nDone. %d rows saved to %s", nrow(done), RAW_FILE))

a <- analyse(done)
saveRDS(a, file.path(RESULTS_DIR, "experiment_summary.rds"))
report(a)
