# =============================================================================
# run_benchmarks.R
# Main analysis script for CASPOC paper benchmarking
#
# Structure:
#   for each dataset:
#     for each approach:
#       run method -> get observed statistic
#       run permutation test -> get p-value
#       save results
#
# Iterations are parallelized across cores using the future/furrr framework.
# Each iteration is fully independent (different seed, fresh data).
#
# Usage:
#   Rscript run_benchmarks.R
#   or source("run_benchmarks.R") in an interactive session
# =============================================================================


# --- Setup -------------------------------------------------------------------

source("R/generate_data.R")
source("R/cv_approaches.R")
source("R/permutation_test.R")
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
N_ITERATIONS <- 5  # increase for final paper ~100

# Number of permutations per iteration for significance testing
N_PERM <- 100  # inrease for final paper ~1000

# Parallelization: number of cores (set to 1 for sequential execution)
N_CORES <- parallelly::availableCores() - 1  # leave one core free

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
  # # TODO: Replace with actual dataset once selected
  # real_microbiome = list(
  #   name = "real_microbiome",
  #   type = "real",
  #   loader = function() {
  #     load_real_dataset("microbiome_example")
  #   }
  # )
)


# =============================================================================
# Single-iteration worker (called in parallel)
# =============================================================================

#' Run all approaches on one dataset for one iteration, with permutation tests
#'
#' @param iter       Iteration number (used as seed)
#' @param ds_name    Dataset name
#' @param data       Data list with X, Y
#' @param n_perm     Number of permutations for significance testing
#'
#' @return A data.frame of results (one row per approach)
run_single_iteration <- function(iter, ds_name, data,
                                n_perm = N_PERM,
                                approaches = APPROACHES,
                                cv_config = CV_CONFIG,
                                hp_grid = HP_GRID) {

  iter_results <- data.frame()

  for (approach in approaches) {

    # --- Run the method on real data ---
    result <- tryCatch(
      run_cv_approach(
        approach        = approach,
        X               = data$X,
        Y               = data$Y,
        ncomp           = cv_config$ncomp,
        num_folds       = cv_config$num_folds,
        num_repeats     = cv_config$num_repeats,
        num_folds_inner = cv_config$num_folds_inner,
        keepX_options   = hp_grid$keepX_options,
        keepY_options   = hp_grid$keepY_options,
        seed            = iter
      ),
      error = function(e) NULL
    )
    if (is.null(result)) next

    # --- Permutation test for significance ---
    perm_result <- tryCatch(
      compute_permutation_pvalue(
        approach      = approach,
        X             = data$X,
        Y             = data$Y,
        observed_stat = result$observed_stat,
        n_perm        = n_perm,
        seed          = iter,
        # pass through CV settings
        ncomp           = cv_config$ncomp,
        num_folds       = cv_config$num_folds,
        num_repeats     = cv_config$num_repeats,
        num_folds_inner = cv_config$num_folds_inner,
        keepX_options   = hp_grid$keepX_options,
        keepY_options   = hp_grid$keepY_options
      ),
      error = function(e) {
        list(perm_pvalue = NA, observed_stat = result$observed_stat,
             null_distribution = NA, n_perm = 0)
      }
    )

    iter_results <- rbind(iter_results,
                          summarise_single_run(result, perm_result, ds_name, iter))
  }

  iter_results
}


# =============================================================================
# Main analysis loop (parallelized across iterations)
# =============================================================================

run_all_benchmarks <- function(n_cores = N_CORES) {

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
      # ----- Simulated datasets: run N_ITERATIONS in parallel -----

      message(sprintf("  Running %d iterations x %d permutations across %d cores...",
                      N_ITERATIONS, N_PERM, n_cores))

      # Pre-capture variables for worker closures
      generator_fn <- ds$generator
      ds_name_local <- ds_name
      n_perm_local <- N_PERM
      cv_config_local <- CV_CONFIG
      hp_grid_local <- HP_GRID
      approaches_local <- APPROACHES

      iter_results <- future_map_dfr(
        seq_len(N_ITERATIONS),
        function(iter) {
          # Source all helper files inside each worker process
          source("R/generate_data.R")
          source("R/cv_approaches.R")
          source("R/permutation_test.R")
          source("R/evaluate_results.R")

          data <- generator_fn(seed = iter)
          run_single_iteration(iter, ds_name_local, data,
                               n_perm     = n_perm_local,
                               approaches = approaches_local,
                               cv_config  = cv_config_local,
                               hp_grid    = hp_grid_local)
        },
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
      )

      all_results <- rbind(all_results, iter_results)
      message(sprintf("  Done: %d rows collected", nrow(iter_results)))

    } else if (ds$type == "real") {
      # ----- Real datasets: single run with permutation test -----

      data <- tryCatch(ds$loader(), error = function(e) {
        message(sprintf("  SKIPPING: %s", e$message))
        NULL
      })
      if (is.null(data)) next

      message(sprintf("  Running with %d permutations...", N_PERM))
      real_results <- run_single_iteration(42, ds_name, data,
                                           n_perm     = N_PERM,
                                           approaches = APPROACHES,
                                           cv_config  = CV_CONFIG,
                                           hp_grid    = HP_GRID)
      all_results <- rbind(all_results, real_results)
      message(sprintf("  Done: %d rows collected", nrow(real_results)))
    }
  }

  # --- Save and summarise ---
  save_results(all_results, "benchmark_raw_results.rds", RESULTS_DIR)

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
  cat(sprintf("Config: %d iterations, %d permutations, %d approaches, %d datasets, %d cores\n",
              N_ITERATIONS, N_PERM, length(APPROACHES), length(datasets), N_CORES))
}
