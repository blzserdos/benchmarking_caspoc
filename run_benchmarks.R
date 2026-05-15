# =============================================================================
# run_benchmarks.R
# Main analysis script for CASPOC paper benchmarking
#
# Strategy:
#   1. Build a flat grid of all jobs: dataset × iteration × approach × perm_id
#      (perm_id = 0 is the real data run, 1..N_PERM are permuted runs)
#   2. Each job independently generates data from seed, optionally permutes Y,
#      runs one CV method, and returns a single result row.
#   3. All jobs are dispatched via future_map_dfr for maximum load balancing.
#   4. Post-processing: group by (dataset, iteration, approach), compute
#      permutation p-values from collected null statistics.
#
# This design ensures perfect load balancing — fast naïve CV jobs finish
# quickly, freeing cores for slower nested CV jobs. No core idles while
# work remains.
#
# Usage:
#   Rscript run_benchmarks.R
#   or source("run_benchmarks.R") in an interactive session
# =============================================================================


# --- Setup -------------------------------------------------------------------

source("R/generate_data.R")
source("R/cv_approaches.R")
source("R/evaluate_results.R")

library(future)
library(furrr)


# --- Configuration -----------------------------------------------------------

# Approaches to compare
APPROACHES <- c("naive_cv", "repeated_cv", "nested_cv", "caspoc")

# Shared CV settings
CV_CONFIG <- list(
  ncomp       = 1,
  num_folds   = 10,
  num_repeats = 11,        # for repeated_cv, nested_cv, and caspoc
  num_folds_inner = 5      # for nested_cv
)

# Hyperparameter grid (kept small for initial development; expand later)
HP_GRID <- list(
  keepX_options = c(5, 10, 20, 50),
  keepY_options = c(5, 10, 20)
)

# Number of simulation iterations per dataset
N_ITERATIONS <- 3  # increase for final paper (e.g. 200-500)

# Number of permutations per iteration for significance testing
N_PERM <- 10  # 100 gives 1% resolution on p-values

# Parallelization: number of cores (set to 1 for sequential execution)
N_CORES <- parallelly::availableCores() - 1  # leave one core free

# Chunk size for future_map_dfr: how many jobs per worker dispatch
# Higher = less scheduling overhead, lower = better load balancing
CHUNK_SIZE <- 10

# Output directory
RESULTS_DIR <- "results"


# =============================================================================
# Define datasets
# =============================================================================

datasets <- list(

  # --- Simulation 1: Null data (Type I error) ---
  sim_null = list(
    name = "sim_null",
    type = "simulated",
    generator = function(seed) {
      generate_null_data(
        n = 100, p = 200, q = 50,
        cor_x = 0, cor_y = 0,
        seed = seed
      )
    }
  ),

  # --- Simulation 2: Signal data (Type II error / power) ---
  sim_signal = list(
    name = "sim_signal",
    type = "simulated",
    generator = function(seed) {
      generate_signal_data(
        n = 100, p = 200, q = 50,
        n_comp_true = 1,
        n_relevant_x = 20,
        n_relevant_y = 10,
        signal_strength = 1.0,
        seed = seed
      )
    }
  )#,

  # # --- Real data 1: Breast TCGA (mRNA vs protein) ---
  # real_breast_tcga = list(
  #   name = "real_breast_tcga",
  #   type = "real",
  #   loader = function() {
  #     load_real_dataset("breast_tcga")
  #   }
  # ),

  # # --- Real data 2: Microbiome example ---
  # real_microbiome = list(
  #   name = "real_microbiome",
  #   type = "real",
  #   loader = function() {
  #     load_real_dataset("microbiome_example")
  #   }
  # )
)


# =============================================================================
# Single-job worker: one CV run
# =============================================================================

#' Execute a single CV job
#'
#' Each job generates data from seed, optionally permutes Y, and runs one
#' CV method. This is the atomic unit of parallelism.
#'
#' @param job         A one-row data.frame with columns: dataset, iteration,
#'                    approach, perm_id
#' @param ds_config   Dataset config (generator function, etc.)
#' @param cv_config   CV settings list
#' @param hp_grid     Hyperparameter grid list
#'
#' @return A one-row data.frame with job identifiers + result
run_single_job <- function(job, ds_config, cv_config, hp_grid) {

  # Generate data (deterministic from iteration seed)
  data <- ds_config$generator(seed = job$iteration)

  # If this is a permutation run, shuffle Y
  if (job$perm_id > 0) {
    set.seed(job$iteration * 1000 + job$perm_id)
    perm_idx <- sample(nrow(data$Y))
    data$Y <- data$Y[perm_idx, , drop = FALSE]
  }

  # Run the CV method
  result <- tryCatch(
    run_cv_approach(
      approach        = job$approach,
      X               = data$X,
      Y               = data$Y,
      ncomp           = cv_config$ncomp,
      num_folds       = cv_config$num_folds,
      num_repeats     = cv_config$num_repeats,
      num_folds_inner = cv_config$num_folds_inner,
      keepX_options   = hp_grid$keepX_options,
      keepY_options   = hp_grid$keepY_options,
      seed            = job$iteration
    ),
    error = function(e) NULL
  )

  # Return result row
  data.frame(
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
}


# =============================================================================
# Post-processing: compute permutation p-values
# =============================================================================

#' Compute permutation p-values from flat results
#'
#' Groups by (dataset, iteration, approach), extracts the real run (perm_id=0)
#' and permuted runs (perm_id>0), then computes empirical p-value.
#'
#' @param flat_results  data.frame from the parallel run (one row per job)
#'
#' @return data.frame with one row per (dataset, iteration, approach),
#'         including observed_stat, perm_pvalue, and HP selections
compute_pvalues <- function(flat_results) {

  # Split into real runs and permutation runs to avoid dplyr column ambiguity
  real_runs <- flat_results %>%
    filter(perm_id == 0) %>%
    select(dataset, iteration, approach,
           observed_stat, selected_keepX, selected_keepY, runtime_sec)

  perm_runs <- flat_results %>%
    filter(perm_id > 0) %>%
    group_by(dataset, iteration, approach) %>%
    summarise(
      # Permutation p-value: (n_extreme + 1) / (n_valid + 1)
      # Conservative estimate (Phipson & Smyth, 2010)
      n_perm_completed = sum(!is.na(observed_stat)),
      null_stats       = list(observed_stat[!is.na(observed_stat)]),
      total_perm_runtime = sum(runtime_sec, na.rm = TRUE),
      .groups = "drop"
    )

  # Join and compute p-values
  result <- real_runs %>%
    left_join(perm_runs, by = c("dataset", "iteration", "approach")) %>%
    rowwise() %>%
    mutate(
      perm_pvalue = {
        n_extreme <- sum(null_stats >= observed_stat)
        n_valid   <- length(null_stats)
        (n_extreme + 1) / (n_valid + 1)
      },
      total_runtime_sec = runtime_sec + total_perm_runtime
    ) %>%
    ungroup() %>%
    select(dataset, iteration, approach,
           observed_stat, selected_keepX, selected_keepY,
           perm_pvalue, n_perm_completed, runtime_sec, total_runtime_sec)

  result
}


# =============================================================================
# Main analysis (parallelized with flat job grid)
# =============================================================================

run_all_benchmarks <- function(n_cores = N_CORES, chunk_size = CHUNK_SIZE) {

  # --- Set up parallel backend ---
  if (n_cores > 1) {
    plan(multisession, workers = n_cores)
    message(sprintf("Parallelizing with %d cores", n_cores))
  } else {
    plan(sequential)
    message("Running sequentially (set N_CORES > 1 to parallelize)")
  }
  on.exit(plan(sequential), add = TRUE)

  all_results <- data.frame()

  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    message(sprintf("\n========== Dataset: %s ==========", ds_name))

    if (ds$type == "simulated") {

      # ----- Build flat job grid -----
      job_grid <- expand.grid(
        iteration = seq_len(N_ITERATIONS),
        approach  = APPROACHES,
        perm_id   = 0:N_PERM,  # 0 = real data, 1..N_PERM = permuted
        stringsAsFactors = FALSE
      )
      job_grid$dataset <- ds_name

      n_jobs <- nrow(job_grid)
      message(sprintf("  %d total jobs (%d iterations x %d approaches x %d perm levels)",
                      n_jobs, N_ITERATIONS, length(APPROACHES), N_PERM + 1))
      message(sprintf("  Dispatching with chunk_size = %d across %d cores...",
                      chunk_size, n_cores))

      # ----- Capture closure variables -----
      ds_config_local  <- ds
      cv_config_local  <- CV_CONFIG
      hp_grid_local    <- HP_GRID

      # ----- Run all jobs in parallel -----
      flat_results <- future_map_dfr(
        seq_len(n_jobs),
        function(j) {
          # Source helper files in each worker
          source("R/generate_data.R")
          source("R/cv_approaches.R")

          run_single_job(
            job       = job_grid[j, ],
            ds_config = ds_config_local,
            cv_config = cv_config_local,
            hp_grid   = hp_grid_local
          )
        },
        .options = furrr_options(seed = TRUE, chunk_size = chunk_size),
        .progress = TRUE
      )

      # ----- Save raw flat results -----
      save_results(flat_results,
                   paste0(ds_name, "_flat_results.rds"),
                   RESULTS_DIR)

      # ----- Compute p-values -----
      ds_results <- compute_pvalues(flat_results)
      all_results <- rbind(all_results, ds_results)

      message(sprintf("  Done: %d result rows from %d jobs", nrow(ds_results), n_jobs))

    } else if (ds$type == "real") {
      # ----- Real datasets: single iteration with permutation test -----

      data <- tryCatch(ds$loader(), error = function(e) {
        message(sprintf("  SKIPPING: %s", e$message))
        NULL
      })
      if (is.null(data)) next

      # Build job grid for single iteration
      job_grid <- expand.grid(
        iteration = 42L,
        approach  = APPROACHES,
        perm_id   = 0:N_PERM,
        stringsAsFactors = FALSE
      )
      job_grid$dataset <- ds_name

      n_jobs <- nrow(job_grid)
      message(sprintf("  %d jobs (%d approaches x %d perm levels)",
                      n_jobs, length(APPROACHES), N_PERM + 1))

      # For real data, pass the loaded data directly (no generator)
      # Create a dummy generator that returns the pre-loaded data
      data_local <- data
      ds_config_real <- list(
        generator = function(seed) data_local
      )
      cv_config_local <- CV_CONFIG
      hp_grid_local   <- HP_GRID

      flat_results <- future_map_dfr(
        seq_len(n_jobs),
        function(j) {
          source("R/generate_data.R")
          source("R/cv_approaches.R")

          run_single_job(
            job       = job_grid[j, ],
            ds_config = ds_config_real,
            cv_config = cv_config_local,
            hp_grid   = hp_grid_local
          )
        },
        .options = furrr_options(seed = TRUE, chunk_size = chunk_size),
        .progress = TRUE
      )

      save_results(flat_results,
                   paste0(ds_name, "_flat_results.rds"),
                   RESULTS_DIR)

      ds_results <- compute_pvalues(flat_results)
      all_results <- rbind(all_results, ds_results)

      message(sprintf("  Done: %d result rows from %d jobs", nrow(ds_results), n_jobs))
    }
  }

  # --- Save and summarise ---
  save_results(all_results, "benchmark_results.rds", RESULTS_DIR)

  summary <- compute_summary_stats(all_results)
  save_results(summary, "benchmark_summary.rds", RESULTS_DIR)

  print_comparison_table(summary)

  invisible(all_results)
}


# =============================================================================
# Run
# =============================================================================

if (!interactive()) {
  run_all_benchmarks()
} else {
  cat("Benchmarking script loaded. Call run_all_benchmarks() to execute.\n")
  cat(sprintf("Config: %d iterations, %d permutations, %d approaches, %d datasets\n",
              N_ITERATIONS, N_PERM, length(APPROACHES), length(datasets)))
  cat(sprintf("Total jobs per simulated dataset: %d\n",
              N_ITERATIONS * length(APPROACHES) * (N_PERM + 1)))
  cat(sprintf("Cores: %d, chunk_size: %d\n", N_CORES, CHUNK_SIZE))
}
