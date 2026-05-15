#!/usr/bin/env Rscript
# =============================================================================
# collect_results.R
# Collects chunk files, computes permutation p-values, and produces summary.
#
# Run after all SLURM array tasks have completed:
#   cd benchmarking_caspoc
#   Rscript cluster/collect_results.R
# =============================================================================

library(dplyr)
source("R/evaluate_results.R")
source("cluster/config.R")

# --- Collect all chunk files ---
chunk_dir   <- "cluster/chunks"
chunk_files <- list.files(chunk_dir, pattern = "^chunk_.*\\.rds$", full.names = TRUE)

if (length(chunk_files) == 0) {
  stop("No chunk files found in ", chunk_dir, ". Did the SLURM jobs complete?")
}

message(sprintf("Found %d chunk files", length(chunk_files)))

flat_results <- do.call(rbind, lapply(chunk_files, readRDS))
message(sprintf("Total rows: %d", nrow(flat_results)))

# --- Verify completeness ---
expected_grid <- build_job_grid()
n_expected    <- nrow(expected_grid)
n_actual      <- nrow(flat_results)

if (n_actual < n_expected) {
  # Find missing jobs
  flat_results$job_key <- paste(flat_results$dataset, flat_results$iteration,
                                flat_results$approach, flat_results$perm_id, sep = "|")
  expected_grid$job_key <- paste(expected_grid$dataset, expected_grid$iteration,
                                 expected_grid$approach, expected_grid$perm_id, sep = "|")

  missing <- setdiff(expected_grid$job_key, flat_results$job_key)
  warning(sprintf("%d of %d jobs missing (%d completed). Missing jobs:\n%s",
                  length(missing), n_expected, n_actual,
                  paste(head(missing, 20), collapse = "\n")))

  flat_results$job_key <- NULL
  expected_grid$job_key <- NULL
} else {
  message(sprintf("All %d jobs completed", n_expected))
}

# --- Check for duplicates ---
n_unique <- nrow(distinct(flat_results, dataset, iteration, approach, perm_id))
if (n_unique < n_actual) {
  warning(sprintf("%d duplicate rows detected. Keeping first occurrence.",
                  n_actual - n_unique))
  flat_results <- distinct(flat_results, dataset, iteration, approach, perm_id,
                           .keep_all = TRUE)
}

# --- Compute permutation p-values ---
message("Computing permutation p-values...")

# Split into real and permutation runs
real_runs <- flat_results %>%
  filter(perm_id == 0) %>%
  select(dataset, iteration, approach,
         observed_stat, selected_keepX, selected_keepY, runtime_sec)

perm_runs <- flat_results %>%
  filter(perm_id > 0) %>%
  group_by(dataset, iteration, approach) %>%
  summarise(
    n_perm_completed   = sum(!is.na(observed_stat)),
    null_stats         = list(observed_stat[!is.na(observed_stat)]),
    total_perm_runtime = sum(runtime_sec, na.rm = TRUE),
    .groups = "drop"
  )

# Join and compute p-values
all_results <- real_runs %>%
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

# --- Save results ---
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Save flat results (all jobs)
saveRDS(flat_results, file.path(results_dir, "benchmark_flat_results.rds"))
message("Saved: results/benchmark_flat_results.rds")

# Save processed results (one row per iteration x approach)
saveRDS(all_results, file.path(results_dir, "benchmark_results.rds"))
message("Saved: results/benchmark_results.rds")

# Save summary
summary <- compute_summary_stats(all_results)
saveRDS(summary, file.path(results_dir, "benchmark_summary.rds"))
message("Saved: results/benchmark_summary.rds")

# --- Print results ---
print_comparison_table(summary)

# --- Completeness report ---
cat("\n=== COMPLETENESS REPORT ===\n")
completeness <- all_results %>%
  group_by(dataset, approach) %>%
  summarise(
    n_iterations     = n(),
    min_perms        = min(n_perm_completed),
    max_perms        = max(n_perm_completed),
    median_perms     = median(n_perm_completed),
    total_cpu_hours  = sum(total_runtime_sec, na.rm = TRUE) / 3600,
    .groups = "drop"
  )
print(as.data.frame(completeness))
