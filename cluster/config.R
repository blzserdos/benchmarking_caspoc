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
N_ITERATIONS <- 5
N_PERM       <- 5

# Dataset definitions (generators reference functions from generate_data.R)
datasets <- list(
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
  )
)

# Build the full job grid (one row per atomic CV run)
# The grid is shuffled deterministically so that each chunk gets a mix of
# fast (naive_cv) and slow (nested_cv) jobs, making walltime more uniform.
build_job_grid <- function() {
  grids <- list()
  for (ds_name in names(datasets)) {
    ds <- datasets[[ds_name]]
    if (ds$type == "simulated") {
      g <- expand.grid(
        iteration = seq_len(N_ITERATIONS),
        approach  = APPROACHES,
        perm_id   = 0:N_PERM,
        stringsAsFactors = FALSE
      )
      g$dataset <- ds_name
      grids[[ds_name]] <- g
    }
  }
  full_grid <- do.call(rbind, grids)

  # Deterministic shuffle so chunk assignment is reproducible
  set.seed(42)
  full_grid <- full_grid[sample(nrow(full_grid)), ]
  rownames(full_grid) <- NULL

  full_grid
}
