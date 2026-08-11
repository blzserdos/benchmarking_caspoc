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

# Within-block correlation structure: block-diagonal compound symmetry with
# heterogeneous module strengths. Module k has size ~ k^-alpha and within-module
# correlation rho * k^-beta; correlation across modules is zero.
#
# Fitted to the breast.TCGA mRNA block over the first 20 eigenvalue fractions.
# Matches DEFAULT_BLOCK_STRUCTURE in R/generate_data.R, spelled out here so this
# file stands alone (collect_results.R sources it without generate_data.R).
#
# Set to NULL for the original independent-feature simulation. Not a cosmetic
# change: nuisance module factors are independent across blocks but high
# variance within a block, and at n = 100 they shift the detectability threshold
# right by roughly 3-4x relative to independent noise.
#
# beta = 0 (equal-strength modules) is a poor idealisation -- it produces a few
# equal eigenvalue spikes then a cliff, and measured ~2x lower power at
# signal_strength = 4 than the fitted decaying spectrum. See R/generate_data.R.
BLOCK_STRUCTURE <- list(n_modules = 8, alpha = 0.5, rho = 0.8,
                        beta = 0.5, frac = 1.0)

# Where the truly relevant features sit relative to the modules: "spread"
# (distributed) or "aligned" (all inside module 1). This is NOT a minor knob
# once modules differ in strength:
#
#   spread  : oracle cor at s=1 is 0.50 +/- 0.03  |  power at s=4 = 0.38
#   aligned : oracle cor at s=1 is 0.58 +/- 0.15  |  power at s=4 = 0.72
#
# "spread" is the default because it preserves the s/(s+1) identity, so a power
# curve indexed by signal_strength is interpretable. Under "aligned" all 20
# relevant features land in module 1 (rho = 0.8) and much of the variation at
# fixed s is just how the loadings happen to interact with that module factor.
SIGNAL_ALIGNMENT <- "spread"

# Signal strengths to evaluate on signal datasets.
#
# Under BLOCK_STRUCTURE the power transition spans roughly s = 2 to 10 (at
# n = 100, p = 200, q = 50, 20/10 relevant features). Real omics block pairs
# calibrate far above this — breast.TCGA pairs match s ~ 20-70 depending on the
# assumed sparsity — so every method saturates at power 1 on realistic data.
# This grid therefore characterises the hard regime where the methods differ,
# which is the claim the sweep can actually support.
SIGNAL_STRENGTHS <- c(1, 2, 3, 4, 6, 8, 12)

# Optional single-strength mode, for clusters whose MaxArraySize cannot hold the
# full sweep in one array. When STRENGTH_INDEX (0-based) is set in the
# environment, this process sees only that one signal strength, so the sweep can
# be submitted as one array per strength:
#
#   for i in $(seq 0 6); do
#     sbatch --array=1-408 --export=ALL,STRENGTH_INDEX=$i cluster/submit.sh
#   done
#
# sim_null carries no signal strength, so it is attached to index 0 only rather
# than being re-run by every submission.
#
# IMPORTANT: leave STRENGTH_INDEX unset when running collect_results.R — it needs
# the full grid to verify completeness.
STRENGTH_INDEX <- Sys.getenv("STRENGTH_INDEX", unset = "")
if (nzchar(STRENGTH_INDEX)) {
  .i <- suppressWarnings(as.integer(STRENGTH_INDEX)) + 1L
  if (is.na(.i) || .i < 1L || .i > length(SIGNAL_STRENGTHS)) {
    stop(sprintf("STRENGTH_INDEX must be 0..%d, got '%s'",
                 length(SIGNAL_STRENGTHS) - 1L, STRENGTH_INDEX))
  }
  SIGNAL_STRENGTHS <- SIGNAL_STRENGTHS[.i]
  INCLUDE_NULL_DATASET <- (.i == 1L)
  message(sprintf("STRENGTH_INDEX=%s -> signal_strength = %g, sim_null %s",
                  STRENGTH_INDEX, SIGNAL_STRENGTHS,
                  if (INCLUDE_NULL_DATASET) "included" else "skipped"))
} else {
  INCLUDE_NULL_DATASET <- TRUE
}

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
        structure = BLOCK_STRUCTURE,
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
        structure = BLOCK_STRUCTURE,
        signal_alignment = SIGNAL_ALIGNMENT,
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
    if (ds_name == "sim_null" && !INCLUDE_NULL_DATASET) next

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
