#!/usr/bin/env Rscript
# =============================================================================
# run_chunk.R
# Runs a chunk of the job grid identified by SLURM array task ID.
#
# Usage:
#   Rscript cluster/run_chunk.R <ARRAY_TASK_ID> [CHUNK_SIZE]
#
# Each task processes CHUNK_SIZE jobs from the full grid and saves a partial
# results file to cluster/chunks/chunk_<ID>.rds
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript cluster/run_chunk.R <ARRAY_TASK_ID> [CHUNK_SIZE]")

task_id    <- as.integer(args[1])
chunk_size <- if (length(args) >= 2) as.integer(args[2]) else 200

# --- Source helpers ---
source("R/generate_data.R")
source("R/cv_approaches.R")
source("cluster/config.R")

# --- Build job grid and extract this chunk ---
full_grid <- build_job_grid()
n_jobs    <- nrow(full_grid)

start_idx <- (task_id - 1) * chunk_size + 1
end_idx   <- min(task_id * chunk_size, n_jobs)

if (start_idx > n_jobs) {
  message(sprintf("Task %d: no jobs (start=%d > n_jobs=%d). Exiting.", task_id, start_idx, n_jobs))
  quit(save = "no", status = 0)
}

my_jobs <- full_grid[start_idx:end_idx, ]
message(sprintf("Task %d: processing jobs %d-%d (%d jobs)", task_id, start_idx, end_idx, nrow(my_jobs)))

# --- Run each job ---
results <- vector("list", nrow(my_jobs))

for (j in seq_len(nrow(my_jobs))) {
  job <- my_jobs[j, ]

  # Generate data from seed
  ds <- datasets[[job$dataset]]
  data <- ds$generator(seed = job$iteration)

  # Permute Y if this is a permutation run
  if (job$perm_id > 0) {
    set.seed(job$iteration * 1000 + job$perm_id)
    perm_idx <- sample(nrow(data$Y))
    data$Y <- data$Y[perm_idx, , drop = FALSE]
  }

  # Run CV method
  result <- tryCatch(
    run_cv_approach(
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
    ),
    error = function(e) NULL
  )

  results[[j]] <- data.frame(
    dataset        = job$dataset,
    iteration      = job$iteration,
    approach       = job$approach,
    perm_id        = job$perm_id,
    observed_stat  = if (!is.null(result)) result$observed_stat else NA_real_,
    selected_keepX = if (!is.null(result)) result$selected_keepX else NA_real_,
    selected_keepY = if (!is.null(result)) result$selected_keepY else NA_real_,
    runtime_sec    = if (!is.null(result)) result$runtime else NA_real_,
    stringsAsFactors = FALSE
  )

  # Progress every 10 jobs
  if (j %% 10 == 0) {
    message(sprintf("  Task %d: %d/%d jobs done", task_id, j, nrow(my_jobs)))
  }
}

# --- Combine and save ---
chunk_results <- do.call(rbind, results)

out_dir <- "cluster/chunks"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, sprintf("chunk_%04d.rds", task_id))
saveRDS(chunk_results, out_file)
message(sprintf("Task %d: saved %d rows to %s", task_id, nrow(chunk_results), out_file))
