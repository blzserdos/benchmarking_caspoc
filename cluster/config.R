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
#
# N_DATASETS counts INDEPENDENT SIMULATED DATASETS (replicate draws from the
# generative model). Not to be confused with CV_CONFIG$num_repeats (repeats of
# the fold partition WITHIN one dataset) or with sPLS's own max.iter/tol (the
# NIPALS convergence loop inside a single fit).
N_DATASETS   <- 100
N_PERM       <- 100

# Within-block correlation structure: OVERLAPPING FACTORS. A feature loads on
# several latent factors rather than belonging to exactly one module, so
# cor(i,j) takes a continuum of values. Spelled out here rather than referencing
# DEFAULT_OVERLAP_STRUCTURE because collect_results.R sources this file without
# generate_data.R.
#
# Set to NULL for the original independent-feature simulation. Not a cosmetic
# change: nuisance factors are independent across blocks but high variance
# within a block, and at n = 100 they shift the detectability threshold right
# substantially relative to independent noise.
#
#   K        number of latent factors
#   gamma    spectrum slope (factor weights w_k ~ k^-gamma)
#   pi_join  membership density; ~1/K would give roughly one factor per feature,
#            i.e. the old disjoint-module structure
#   h2_mean  mean communality -- share of a feature's variance that is shared
#
# The goal is omics-LIKE data for benchmarking CV methods, not a replica of any
# dataset, so these are round numbers. See R/generate_data.R for the rationale.
BLOCK_STRUCTURE <- list(type = "overlap", K = 20, gamma = 1.5,
                        pi_join = 0.5, h2_mean = 0.6)

# The previous structure, disjoint modules, kept for reference and for the
# robustness sweep. Its three qualitative artifacts are why it was replaced:
# 86.3% of feature pairs exactly zero, negative correlations structurally
# impossible, and a spectrum that falls off a cliff after n_modules eigenvalues.
#
#                    mean|c|   RMS   median|c|   %neg   spectrum RMSE
#   real mRNA         0.200   0.249    0.173    44.0%        --
#   disjoint          0.134   0.235    0.066    43.2%       0.91
#   overlap           0.198   0.248    0.162    50.0%       0.67
#
# DISJOINT_STRUCTURE <- list(n_modules = 8, alpha = 0.5, rho = 0.8,
#                            beta = 0.5, frac = 1.0)

# Where the truly relevant features sit relative to the modules: "spread"
# (distributed) or "aligned" (all inside module 1).
#
# INERT under the overlapping structure: factor membership is random per
# feature, so a contiguous block of relevant features is already spread across
# the factor structure and the two modes coincide (verified in
# diagnostics/check_overlap.R). Kept because it still applies if BLOCK_STRUCTURE
# is switched back to disjoint modules, where it mattered a great deal:
#
#   spread  : oracle cor at s=1 is 0.50 +/- 0.03  |  power at s=4 = 0.38
#   aligned : oracle cor at s=1 is 0.58 +/- 0.15  |  power at s=4 = 0.72
SIGNAL_ALIGNMENT <- "spread"

# Signal strengths to evaluate on signal datasets.
#
# Measured under the overlapping structure at n = 100, p = 200, q = 50, 20/10
# relevant features (diagnostics/final_measure.R, 200 reps, proxy statistic):
#
#   s     : 1    1.5  2    2.5  3    4    5    6    8    10   12   16   20
#   power : 0.06 0.03 0.07 0.15 0.28 0.53 0.57 0.75 0.89 0.93 0.97 0.99 1.00
#
# The grid below samples that densely where it is steep (2-6) and anchors both
# ends. Note the transition sits in essentially the same place as under the
# previous disjoint-module structure, so this grid did not need re-placing.
#
# Real omics block pairs calibrate far above it: breast.TCGA pairs match
# s ~ 37-48 (diagnostics/calibrate3.R), where every method saturates at power 1.
# The sweep therefore characterises the HARD REGIME where the methods differ,
# which is the claim it can actually support.
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
      iteration       = seq_len(N_DATASETS),
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
