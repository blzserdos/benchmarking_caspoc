# =============================================================================
# config.R
# Shared configuration for benchmarking (used by both local and cluster runs)
# =============================================================================

# Approaches to compare
APPROACHES <- c("naive_cv", "repeated_cv", "nested_cv", "caspoc")

# Shared CV settings
CV_CONFIG <- list(
  ncomp       = 1,
  num_folds   = 10,
  num_repeats = 11,
  num_folds_inner = 5
)

# Hyperparameter grid
HP_GRID <- list(
  keepX_options = c(5, 10, 20, 50),
  keepY_options = c(5, 10, 20)
)

# Simulation settings
N_ITERATIONS <- 100
N_PERM       <- 100

# Signal strengths to evaluate on signal datasets.
# A single value reproduces the previous behaviour; provide a vector to sweep
# power as a function of signal strength, e.g. c(0.25, 0.5, 1.0, 2.0).
SIGNAL_STRENGTHS <- c(1.0)

# Dataset definitions (generators reference functions from generate_data.R)
#
# Per-dataset fields:
#   run_perm         TRUE/FALSE — whether to run permutation tests for this
#                    dataset. Disabled for sim_null to save compute (FPR is
#                    expected ~0.05 a priori and the main result here is the
#                    mean test statistic).
#   signal_strengths Vector of signal strengths to sweep over (only used by
#                    signal datasets). Use NA_real_ for datasets where it
#                    doesn't apply.
#   generator        function(seed, signal_strength) -> list(X, Y, ...).
#                    signal_strength is ignored by null generators.
datasets <- list(
  sim_null = list(
    name = "sim_null",
    type = "simulated",
    run_perm = FALSE,
    signal_strengths = NA_real_,
    generator = function(seed, signal_strength = NA) {
      generate_null_data(
        n = 100, p = 200, q = 50,
        cor_x = 0, cor_y = 0,
        seed = seed
      )
    }
  ),
  sim_signal = list(
    name = "sim_signal",
    type = "simulated",
    run_perm = TRUE,
    signal_strengths = SIGNAL_STRENGTHS,
    generator = function(seed, signal_strength) {
      generate_signal_data(
        n = 100, p = 200, q = 50,
        n_comp_true = 1,
        n_relevant_x = 20,
        n_relevant_y = 10,
        signal_strength = signal_strength,
        seed = seed
      )
    }
  )
)

# Build the full job grid (one row per atomic CV run)
# The grid is shuffled deterministically so that each chunk gets a mix of
# fast (naive_cv) and slow (nested_cv) jobs, making walltime more uniform.
build_job_grid <- function() {
  grids <- list()
  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    if (ds$type != "simulated") next

    perm_ids <- if (isTRUE(ds$run_perm)) 0:N_PERM else 0L
    strengths <- if (length(ds$signal_strengths) > 0 &&
                     !all(is.na(ds$signal_strengths))) {
      ds$signal_strengths
    } else {
      NA_real_
    }

    g <- expand.grid(
      iteration       = seq_len(N_ITERATIONS),
      approach        = APPROACHES,
      perm_id         = perm_ids,
      signal_strength = strengths,
      stringsAsFactors = FALSE
    )
    g$dataset <- ds_name
    grids[[ds_name]] <- g
  }
  full_grid <- do.call(rbind, grids)

  # Deterministic shuffle so chunk assignment is reproducible
  set.seed(42)
  full_grid <- full_grid[sample(nrow(full_grid)), ]
  rownames(full_grid) <- NULL

  full_grid
}

# Helper to print job count + recommended SLURM array size. Run this on the
# login node after editing config to size submit.sh:
#   Rscript -e 'source("cluster/config.R"); print_grid_summary(200)'
print_grid_summary <- function(chunk_size = 200) {
  g <- build_job_grid()
  n <- nrow(g)
  n_tasks <- ceiling(n / chunk_size)
  cat(sprintf("Total jobs: %d\n", n))
  cat(sprintf("Chunk size: %d -> %d array tasks\n", chunk_size, n_tasks))
  cat("By dataset:\n")
  print(as.data.frame(table(dataset = g$dataset)))
  invisible(list(n_jobs = n, n_tasks = n_tasks))
}
