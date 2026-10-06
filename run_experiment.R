#!/usr/bin/env Rscript
# =============================================================================
# run_experiment.R
# Pooled-null benchmark of the four CV approaches. Runs on a laptop.
#
#   Rscript run_experiment.R                 # run (resumable)
#   Rscript run_experiment.R --analyse       # re-print results, no computation
#
# THE DESIGN  (DESIGN = "demo", the default)
#
# A demonstration of CASPOC, not a complete benchmark. Three claims, nothing
# more:
#
#   1. naive_cv and repeated_cv report inflated association strength;
#      nested_cv and caspoc do not.
#   2. caspoc does not achieve that by being insensitive -- its estimate rises
#      with the true signal.
#   3. caspoc costs less than nested_cv (fit counts + recorded runtimes).
#
# One factorial sweep, one CV run per cell:
#
#   4 approaches x N_DATASETS datasets x (sim_null + STRENGTHS_DEMO)
#
# ONE READOUT: the mean reported statistic per method per condition. Claim 1 is
# carried entirely by the sim_null row -- Cov(X, Y) = 0 there by construction,
# so whatever a method prints IS its bias, and training-set-size differences
# between the methods do not matter because the truth is 0 for all of them.
# Claim 2 is carried by reading caspoc across conditions: an estimate that
# tracks s is responsive, which is all that needs showing.
#
# WHY NO PERMUTATIONS OR p-VALUES. Both exist only to turn a statistic into a
# yes/no decision, and a decision needs a threshold. Measuring the MEAN of the
# statistic declares nothing, so it needs no threshold -- and therefore no
# permutation distribution, no critical value, no rejection rate, no FPR. This
# is the single change that removes almost all of the machinery.
#
# WHY THESE SIGNAL STRENGTHS. s = 4 is the hard regime where the methods can
# separate at all. s = 40 is where real omics data actually sits: breast.TCGA
# block pairs calibrate to s ~ 37-48 (diagnostics/calibrate3.R), where every
# method detects the association every time. The s = 40 row is therefore the
# load-bearing one for any real-data claim -- if naive_cv still overstates
# there, the paper applies to real analyses; if the gap closes, the scope limit
# is worth knowing before a reviewer finds it.
#
# -----------------------------------------------------------------------------
# DESIGN HISTORY -- set DESIGN <- "benchmark" below to restore the previous one
# -----------------------------------------------------------------------------
#
# "benchmark" (previous default) additionally computed:
#
#   c_M     each method's 95th-percentile sim_null statistic, used as its own
#           alpha = 0.05 threshold. Per method, so a method could not win by
#           inflating -- inflation raised its own bar by the same amount.
#   POWER   the share of signal runs clearing c_M, swept over the full
#           SIGNAL_STRENGTHS grid c(1, 2, 3, 4, 6, 8, 12).
#
# Dropped because a demonstration needs no power curve, and because that sweep
# characterises only the hard regime -- all four methods saturate at power 1.00
# where real data lives, so the power axis cannot support a real-data claim.
# Both are still computed and printed under DESIGN = "benchmark".
#
# See cluster/ for the third design: per-dataset permutation testing (283,200
# jobs, ~7,500 CPU-hours), which answers the narrower practitioner question
# "what power does ONE analyst get from method M plus a permutation test?"
# =============================================================================

suppressMessages({
  library(dplyr)
  library(parallel)
})

source("R/generate_data.R")
source("R/cv_approaches.R")
source("cluster/config.R")   # APPROACHES, CV_CONFIG, HP_GRID, N_DATASETS,
                             # SIGNAL_STRENGTHS, BLOCK_STRUCTURE, datasets

# -----------------------------------------------------------------------------
# DESIGN SWITCH -- the one line to change to revert
# -----------------------------------------------------------------------------
#   "demo"      bias table only, s = 0 / 4 / 40        (current)
#   "benchmark" adds c_M thresholds and power curves over the full grid
#
DESIGN <- "demo"

STRENGTHS_DEMO <- c(4, 40)   # hard regime, then where real omics sits
STRENGTHS_FULL <- SIGNAL_STRENGTHS   # c(1, 2, 3, 4, 6, 8, 12), from config.R

SIGNAL_STRENGTHS <- switch(DESIGN,
  demo      = STRENGTHS_DEMO,
  benchmark = STRENGTHS_FULL,
  stop("DESIGN must be \"demo\" or \"benchmark\", got: ", DESIGN)
)

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

  # sim_null is condition s = 0: one axis for every condition, so the estimates
  # table reads straight across from "no signal" to "what real data looks like".
  raw <- raw %>% mutate(s = ifelse(dataset == "sim_null", 0, signal_strength))

  # --- ESTIMATES: the reported statistic, per method per condition ---
  # At s = 0 the truth is 0, so the mean IS the bias. At s > 0 the mean is read
  # ACROSS methods (does naive still sit above caspoc?) and ACROSS conditions
  # (does caspoc rise with s?) -- no threshold is involved either way.
  est <- raw %>%
    group_by(approach, s) %>%
    summarise(
      n         = sum(!is.na(observed_stat)),
      mean_stat = mean(observed_stat, na.rm = TRUE),
      sd_stat   = sd(observed_stat, na.rm = TRUE),
      .groups   = "drop"
    )

  runtime <- raw %>%
    group_by(approach) %>%
    summarise(mean_runtime = mean(runtime_sec, na.rm = TRUE), .groups = "drop")

  out <- list(est = est, runtime = runtime, design = DESIGN)

  # --- DESIGN = "benchmark" only: thresholds and power (see DESIGN HISTORY) ---
  if (identical(DESIGN, "benchmark")) {
    thresholds <- raw %>%
      dplyr::filter(dataset == "sim_null") %>%
      group_by(approach) %>%
      summarise(
        # type = 1 keeps the threshold an actually observed value, so the
        # nominal level is not inflated by interpolating order statistics.
        crit_value = quantile(observed_stat, 1 - ALPHA, na.rm = TRUE, type = 1),
        .groups = "drop"
      )

    out$power <- raw %>%
      dplyr::filter(dataset == "sim_signal") %>%
      left_join(thresholds, by = "approach") %>%
      group_by(approach, s) %>%
      summarise(power = mean(observed_stat > crit_value, na.rm = TRUE),
                .groups = "drop")
    out$thresholds <- thresholds
  }

  out
}

report <- function(a) {
  conds     <- sort(unique(a$est$s))
  approaches <- a$runtime$approach

  cat("\n=======================================================================\n")
  cat("REPORTED ASSOCIATION STRENGTH  (mean of the statistic each method prints)\n")
  cat("=======================================================================\n\n")

  cat(sprintf("%-14s", "Approach"))
  for (s in conds) cat(sprintf("%16s", if (s == 0) "s=0 (null)" else sprintf("s=%g", s)))
  cat(sprintf("%11s\n", "Time(s)"))
  cat(strrep("-", 14 + 16 * length(conds) + 11), "\n")

  for (ap in approaches) {
    cat(sprintf("%-14s", ap))
    for (s in conds) {
      r <- a$est[a$est$approach == ap & a$est$s == s, ]
      if (nrow(r) == 0) cat(sprintf("%16s", "-"))
      else cat(sprintf("%10.3f%6s", r$mean_stat[1], sprintf("(%.2f)", r$sd_stat[1])))
    }
    cat(sprintf("%11.1f\n", a$runtime$mean_runtime[a$runtime$approach == ap]))
  }

  cat("\n  Cells are mean (SD) over", max(a$est$n), "datasets.\n")
  cat("\n  s = 0    truth is 0, so the number printed IS the bias.\n")
  cat("  s > 0    read ACROSS methods (does naive still sit above caspoc?) and\n")
  cat("           ACROSS columns (does caspoc rise with the true signal?).\n")
  cat("           No threshold, no p-value, no rejection rate is involved.\n")
  cat("  s = 40   where real omics data sits (breast.TCGA ~ s 37-48), so this\n")
  cat("           column is what any real-data claim rests on.\n")

  if (identical(a$design, "benchmark")) {
    cat("\n=======================================================================\n")
    cat("POWER  (DESIGN = \"benchmark\": share of signal runs clearing c_M)\n")
    cat("=======================================================================\n\n")
    ps <- sort(unique(a$power$s))
    cat(sprintf("%-14s %10s", "Approach", "c_M"))
    for (s in ps) cat(sprintf("%8s", sprintf("s=%g", s)))
    cat("\n"); cat(strrep("-", 24 + 8 * length(ps)), "\n")
    for (ap in approaches) {
      cat(sprintf("%-14s %10.3f", ap, a$thresholds$crit_value[a$thresholds$approach == ap]))
      row <- a$power %>% dplyr::filter(approach == ap) %>% arrange(s)
      cat(sprintf("%8.2f", row$power), "\n")
    }
    cat("\n  c_M is each method's own 95th-percentile null statistic, so every\n")
    cat("  method is held to a 0.05 false-positive rate and the power numbers\n")
    cat("  are comparable. 0.05 = no better than chance.\n")
  }
  cat("\n")
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
saveRDS(a, file.path(RESULTS_DIR,
                     if (SMOKE) "smoke_experiment_summary.rds" else "experiment_summary.rds"))
report(a)
