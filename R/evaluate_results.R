# =============================================================================
# evaluate_results.R
# Functions to evaluate and summarise benchmarking results
# =============================================================================

library(dplyr)


# =============================================================================
# Summarise a single run
# =============================================================================

#' Extract a one-row summary from a CV result + permutation test
#'
#' @param result        A result list from run_cv_approach()
#' @param perm_result   A result list from compute_permutation_pvalue()
#' @param dataset_name  Name of the dataset
#' @param iteration     Iteration number (for repeated simulation runs)
#'
#' @return A one-row data.frame
summarise_single_run <- function(result, perm_result, dataset_name, iteration = 1) {
  data.frame(
    dataset        = dataset_name,
    iteration      = iteration,
    approach       = result$approach,
    selected_keepX = result$selected_keepX,
    selected_keepY = result$selected_keepY,
    observed_stat  = result$observed_stat,
    perm_pvalue    = perm_result$perm_pvalue,
    runtime_sec    = result$runtime,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Aggregate across iterations (for simulation studies)
# =============================================================================

#' Compute summary statistics across simulation iterations
#'
#' @param results_df  A data.frame of stacked single-run summaries
#'
#' @return A summary data.frame grouped by dataset x approach
compute_summary_stats <- function(results_df) {
  results_df %>%
    group_by(dataset, approach) %>%
    summarise(
      n_iterations     = n(),

      # Test statistic summaries
      mean_stat        = mean(observed_stat, na.rm = TRUE),
      sd_stat          = sd(observed_stat, na.rm = TRUE),

      # Rejection rate at alpha = 0.05 (based on permutation p-values)
      # Interpretation depends on dataset:
      #   null/permuted data -> this is the false positive rate (should be ~0.05)
      #   signal data        -> this is power (higher is better)
      reject_rate_005  = mean(perm_pvalue < 0.05, na.rm = TRUE),

      # Runtime
      mean_runtime     = mean(runtime_sec, na.rm = TRUE),

      .groups = "drop"
    )
}


# =============================================================================
# Save / load results
# =============================================================================

#' Save results to disk
save_results <- function(results_df, filename, results_dir = "results") {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  filepath <- file.path(results_dir, filename)
  saveRDS(results_df, filepath)
  message("Saved results to: ", filepath)
}

#' Load results from disk
load_results <- function(filename, results_dir = "results") {
  filepath <- file.path(results_dir, filename)
  readRDS(filepath)
}


# =============================================================================
# Formatted output for reporting
# =============================================================================

#' Print a comparison table to the console
#'
#' @param summary_df  Output of compute_summary_stats()
print_comparison_table <- function(summary_df) {
  cat("\n=== BENCHMARKING COMPARISON ===\n\n")

  for (ds in unique(summary_df$dataset)) {
    is_null <- grepl("null|permuted", ds, ignore.case = TRUE)
    rate_label <- if (is_null) "FPR" else "Power"

    cat(sprintf("--- Dataset: %s ---\n", ds))
    sub <- summary_df %>% filter(dataset == ds)

    cat(sprintf("%-15s %10s %10s %10s %10s\n",
                "Approach", "Mean_stat", "SD_stat", rate_label, "Time(s)"))
    cat(paste(rep("-", 60), collapse = ""), "\n")

    for (j in seq_len(nrow(sub))) {
      cat(sprintf("%-15s %10.3f %10.3f %10.3f %10.1f\n",
                  sub$approach[j],
                  sub$mean_stat[j],
                  sub$sd_stat[j],
                  sub$reject_rate_005[j],
                  sub$mean_runtime[j]))
    }
    cat("\n")
  }
}
