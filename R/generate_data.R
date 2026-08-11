# =============================================================================
# generate_data.R
# Functions to generate simulated datasets for CASPOC benchmarking
# =============================================================================

library(MASS)  # mvrnorm


# --- Within-block correlation structure --------------------------------------

# Real omics blocks are not collections of independent features: they carry a
# few dominant co-expression / co-occurrence axes. The structure here is
# block-diagonal compound symmetry with HETEROGENEOUS module strengths:
#
#   module k has size ~ k^-alpha  and  within-module correlation rho * k^-beta
#   correlation across modules is zero
#
# Fitted to the breast.TCGA mRNA block by minimising RMSE over the first 20
# eigenvalue fractions at matched n and p:
#   n_modules = 8, alpha = 0.5, rho = 0.8, beta = 0.5   (RMSE 0.83 % points)
#
# beta > 0 matters more than it looks. With equal-strength modules (beta = 0)
# the spectrum is a few equal spikes and then a cliff -- e.g. 3 modules gives
# 19.5/14.7/12.6/1.2/1.2/... against real 18.9/13.1/5.5/5.4/4.2/... That
# idealisation leaves sPLS facing several equally attractive spurious
# directions, and it is measurably pessimistic: power at signal_strength = 4
# was 0.49 under equal modules versus 0.85 under the fitted decaying spectrum.
#
# Module factors are independent ACROSS blocks but high variance WITHIN a
# block, so at n = 100 sPLS is readily distracted by spurious factor-to-factor
# alignments. This is the main reason realistic geometry is harder than
# independent noise.
DEFAULT_BLOCK_STRUCTURE <- list(n_modules = 8, alpha = 0.5, rho = 0.8,
                                beta = 0.5, frac = 1.0)

#' Module sizes: n_modules blocks decaying as k^-alpha over frac * p features
module_sizes <- function(p, n_modules, alpha, frac) {
  w <- (seq_len(n_modules))^(-alpha)
  m <- pmax(2, round(w / sum(w) * frac * p))
  while (sum(m) > p) m[which.max(m)] <- m[which.max(m)] - 1
  m
}

#' Within-module correlations: rho * k^-beta, capped below 1
module_cors <- function(n_modules, rho, beta) {
  pmin(0.95, rho * (seq_len(n_modules))^(-beta))
}

#' Draw an n x p noise matrix with unit marginal variance
#'
#' Three regimes, in precedence order:
#'   structure non-NULL  -> modular correlation (one latent factor per module)
#'   cor_within > 0      -> compound symmetry (all pairs equally correlated)
#'   otherwise           -> independent features
#'
#' The modular case uses a factor representation rather than mvrnorm: it is
#' O(n*p) instead of an eigendecomposition of a p x p matrix, and it keeps each
#' feature's marginal variance at 1 so the signal_strength parametrisation is
#' unchanged by the choice of structure.
#'
#' @param n          Number of samples
#' @param p          Number of features
#' @param cor_within Compound-symmetry correlation (ignored if structure given)
#' @param structure  Optional list(n_modules, alpha, rho, beta, frac). beta may
#'                   be omitted, in which case all modules share correlation rho.
make_noise <- function(n, p, cor_within = 0, structure = NULL) {
  if (!is.null(structure) && structure$rho > 0 && structure$n_modules >= 1) {
    m    <- module_sizes(p, structure$n_modules, structure$alpha, structure$frac)
    beta <- if (is.null(structure$beta)) 0 else structure$beta
    rk   <- module_cors(length(m), structure$rho, beta)
    E <- matrix(rnorm(n * p), nrow = n, ncol = p)
    pos <- 0L
    for (k in seq_along(m)) {
      idx <- (pos + 1L):(pos + m[k]); pos <- pos + m[k]
      E[, idx] <- sqrt(rk[k]) * rnorm(n) + sqrt(1 - rk[k]) * E[, idx]
    }
    return(E)
  }

  if (cor_within > 0) {
    Sigma <- matrix(cor_within, nrow = p, ncol = p)
    diag(Sigma) <- 1
    return(mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma))
  }

  matrix(rnorm(n * p), nrow = n, ncol = p)
}

#' Feature indices carrying true loadings for one component
#'
#' "aligned" places them contiguously, so they all fall inside module 1.
#' "spread" distributes them across modules.
#'
#' This choice is consequential once modules differ in strength. With module 1
#' at rho = 0.8, "aligned" puts every relevant feature in one tightly correlated
#' block, so l'Sigma l swings from draw to draw and signal_strength stops being
#' a tight handle on difficulty: at s = 1 the oracle correlation is 0.58 +/-
#' 0.15 (IQR 0.48-0.69) against the nominal 0.50, and power at s = 4 is 0.72 vs
#' 0.38 for "spread".
#'
#' "spread" is the default: it keeps the s/(s+1) identity (0.50 +/- 0.03), so a
#' power curve indexed by s means what it says. Note "aligned" is not simply the
#' more realistic option — module factors are independent ACROSS blocks, so it
#' models relevant features that happen to share an unrelated block-specific
#' factor, not a shared biological axis driving both blocks. Modelling the
#' latter would mean making the latent Z itself a shared module factor.
relevant_indices <- function(p, n_rel, comp = 1, structure = NULL,
                             alignment = "spread") {
  if (alignment == "spread" && !is.null(structure) && structure$n_modules > 1) {
    m <- module_sizes(p, structure$n_modules, structure$alpha, structure$frac)
    starts <- cumsum(c(0L, head(m, -1)))
    per <- diff(round(seq(0, n_rel, length.out = length(m) + 1)))
    # Offset by component so components use non-overlapping features.
    skip <- (comp - 1) * per
    idx <- unlist(lapply(seq_along(m), function(k) {
      if (per[k] <= 0) return(integer(0))
      take <- skip[k] + seq_len(per[k])
      starts[k] + take[take <= m[k]]
    }))
    return(sort(unique(idx)))
  }

  offset <- (comp - 1) * n_rel
  if (offset + 1 > p) return(integer(0))
  (offset + 1):min(offset + n_rel, p)
}


# --- Null data (Type I error) ------------------------------------------------

#' Generate null data with no association between X and Y
#'
#' X and Y are drawn independently from multivariate normal distributions.
#' Any association found by a method on this data is a false positive.
#'
#' @param n      Number of samples
#' @param p      Number of X features
#' @param q      Number of Y features
#' @param cor_x  Within-block compound-symmetry correlation in X
#' @param cor_y  Within-block compound-symmetry correlation in Y
#' @param structure  Optional list(n_modules, alpha, rho, frac) giving a modular
#'                   within-block correlation; overrides cor_x/cor_y. Use
#'                   DEFAULT_BLOCK_STRUCTURE for the breast.TCGA-calibrated fit.
#'                   NULL reproduces the original independent-feature behaviour.
#' @param seed   Random seed for reproducibility
#'
#' @return A list with components:
#'   \item{X}{n x p matrix}
#'   \item{Y}{n x q matrix}
#'   \item{true_signal}{FALSE (no true association)}
#'   \item{params}{List of generation parameters for logging}
generate_null_data <- function(n, p, q,
                               cor_x = 0,
                               cor_y = 0,
                               structure = NULL,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # X and Y drawn independently of each other; any within-block correlation is
  # applied separately to each, so no cross-block association is induced.
  X <- make_noise(n, p, cor_within = cor_x, structure = structure)
  Y <- make_noise(n, q, cor_within = cor_y, structure = structure)

  colnames(X) <- paste0("X", seq_len(p))
  colnames(Y) <- paste0("Y", seq_len(q))

  list(
    X = X,
    Y = Y,
    true_signal = FALSE,
    params = list(
      n = n, p = p, q = q,
      cor_x = cor_x, cor_y = cor_y,
      structure = structure,
      seed = seed,
      type = "null"
    )
  )
}


# --- Signal data (Type II error / power) -------------------------------------

#' Generate data with a true sparse latent association between X and Y
#'
#' Creates a shared latent variable Z, then constructs:
#'   X = Z %*% t(loadings_x) + noise_x
#'   Y = Z %*% t(loadings_y) + noise_y
#' where loadings are sparse (only a subset of features are nonzero).
#'
#' @param n              Number of samples
#' @param p              Number of X features
#' @param q              Number of Y features
#' @param n_comp_true    Number of true latent components
#' @param n_relevant_x   Number of truly relevant X features per component
#' @param n_relevant_y   Number of truly relevant Y features per component
#' @param signal_strength Variance of the latent signal relative to noise (SNR).
#'   With unit-norm loadings and INDEPENDENT noise the oracle cross-block
#'   correlation is exactly s / (s + 1), and a relevant feature carries
#'   s/n_relevant / (1 + s/n_relevant) of its variance on the shared axis.
#'   Under a modular `structure` the identity becomes s / (s + l'Sigma l): the
#'   noise projected onto the loading direction no longer has unit variance, so
#'   the oracle correlation still centres near s / (s + 1) but varies markedly
#'   from draw to draw (at s = 1: mean 0.53, sd 0.08, against 0.50 / 0.01 for
#'   independent noise).
#' @param cor_x          Background compound-symmetry correlation among X
#' @param cor_y          Background compound-symmetry correlation among Y
#' @param structure      Optional list(n_modules, alpha, rho, frac) giving a
#'   modular within-block correlation; overrides cor_x/cor_y. NULL reproduces
#'   the original independent-feature behaviour.
#' @param signal_alignment "spread" (relevant features distributed across
#'   modules; default, keeps signal_strength interpretable) or "aligned" (all
#'   inside module 1). Only has an effect when structure is supplied. See
#'   relevant_indices() — the choice materially changes power.
#' @param seed           Random seed
#'
#' @return A list with components:
#'   \item{X}{n x p matrix}
#'   \item{Y}{n x q matrix}
#'   \item{true_signal}{TRUE}
#'   \item{true_loadings_x}{p x n_comp_true matrix of true loadings}
#'   \item{true_loadings_y}{q x n_comp_true matrix of true loadings}
#'   \item{true_relevant_x}{Indices of truly relevant X features per component}
#'   \item{true_relevant_y}{Indices of truly relevant Y features per component}
#'   \item{params}{List of generation parameters}
generate_signal_data <- function(n, p, q,
                                 n_comp_true = 1,
                                 n_relevant_x = 10,
                                 n_relevant_y = 5,
                                 signal_strength = 1.0,
                                 cor_x = 0,
                                 cor_y = 0,
                                 structure = NULL,
                                 signal_alignment = c("spread", "aligned"),
                                 seed = NULL) {
  signal_alignment <- match.arg(signal_alignment)
  if (!is.null(seed)) set.seed(seed)

  # --- Generate latent variables ---
  Z <- matrix(rnorm(n * n_comp_true), nrow = n, ncol = n_comp_true)

  # --- Build sparse loadings ---
  # Each component gets its own non-overlapping set of relevant features
  loadings_x <- matrix(0, nrow = p, ncol = n_comp_true)
  loadings_y <- matrix(0, nrow = q, ncol = n_comp_true)

  relevant_x_list <- list()
  relevant_y_list <- list()

  for (comp in seq_len(n_comp_true)) {
    # Select relevant features (non-overlapping across components)
    idx_x <- relevant_indices(p, n_relevant_x, comp, structure, signal_alignment)
    idx_y <- relevant_indices(q, n_relevant_y, comp, structure, signal_alignment)

    # Assign random nonzero loadings to relevant features
    loadings_x[idx_x, comp] <- rnorm(length(idx_x), mean = 0, sd = 1)
    loadings_y[idx_y, comp] <- rnorm(length(idx_y), mean = 0, sd = 1)

    # Normalise loadings to unit length
    loadings_x[, comp] <- loadings_x[, comp] / sqrt(sum(loadings_x[, comp]^2))
    loadings_y[, comp] <- loadings_y[, comp] / sqrt(sum(loadings_y[, comp]^2))

    relevant_x_list[[comp]] <- idx_x
    relevant_y_list[[comp]] <- idx_y
  }

  # --- Build X and Y ---
  # Signal part: Z %*% t(loadings) scaled by signal_strength
  # Noise part: independent, compound symmetry, or modular (see make_noise)
  signal_x <- sqrt(signal_strength) * Z %*% t(loadings_x)
  signal_y <- sqrt(signal_strength) * Z %*% t(loadings_y)

  # Background noise (independent, compound symmetry, or modular)
  noise_x <- make_noise(n, p, cor_within = cor_x, structure = structure)
  noise_y <- make_noise(n, q, cor_within = cor_y, structure = structure)

  X <- signal_x + noise_x
  Y <- signal_y + noise_y

  colnames(X) <- paste0("X", seq_len(p))
  colnames(Y) <- paste0("Y", seq_len(q))

  list(
    X = X,
    Y = Y,
    true_signal = TRUE,
    true_loadings_x = loadings_x,
    true_loadings_y = loadings_y,
    true_relevant_x = relevant_x_list,
    true_relevant_y = relevant_y_list,
    params = list(
      n = n, p = p, q = q,
      n_comp_true = n_comp_true,
      n_relevant_x = n_relevant_x,
      n_relevant_y = n_relevant_y,
      signal_strength = signal_strength,
      cor_x = cor_x, cor_y = cor_y,
      structure = structure,
      signal_alignment = signal_alignment,
      seed = seed,
      type = "signal"
    )
  )
}


# --- Real data loaders -------------------------------------------------------

#' Load and prepare a real dataset for benchmarking
#'
#' @param dataset_name  One of the registered dataset names
#' @param data_dir      Path to raw data directory
#'
#' @return A list with X, Y, true_signal = NA, params
load_real_dataset <- function(dataset_name, data_dir = "data/raw") {

  # TODO: Implement loaders for each real dataset
  # Each loader should return X (n x p) and Y (n x q) matrices

  if (dataset_name == "breast_tcga") {
    # --- Breast TCGA (mRNA vs protein) from mixOmics ---
    # Ships with the mixOmics package
    requireNamespace("mixOmics", quietly = TRUE)
    data("breast.TCGA", package = "mixOmics", envir = environment())
    X <- breast.TCGA$data.train$mrna
    Y <- breast.TCGA$data.train$protein

  } else if (dataset_name == "microbiome_example") {
    # --- Placeholder: microbiome dataset ---
    # TODO: Identify a suitable public microbiome dataset
    #       e.g. from curatedMetagenomicData or HMP
    #       X = microbial abundances, Y = metabolites or host phenotypes
    stop("microbiome_example dataset not yet implemented. ",
         "Please add data to ", data_dir, " and update this loader.")

  } else {
    stop("Unknown dataset: ", dataset_name,
         ". Available: 'breast_tcga', 'microbiome_example'")
  }

  list(
    X = X,
    Y = Y,
    true_signal = NA,  # unknown for real data
    params = list(
      dataset_name = dataset_name,
      n = nrow(X),
      p = ncol(X),
      q = ncol(Y),
      type = "real"
    )
  )
}


# --- Permutation wrapper -----------------------------------------------------

#' Permute Y to destroy true association (for type I error on real data)
#'
#' @param dataset  A dataset list (as returned by generate_* or load_*)
#' @param seed     Random seed for the permutation
#'
#' @return A new dataset list with Y rows shuffled
permute_dataset <- function(dataset, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(dataset$Y)
  perm_idx <- sample(n)

  dataset$Y <- dataset$Y[perm_idx, , drop = FALSE]
  dataset$true_signal <- FALSE
  dataset$params$permuted <- TRUE
  dataset$params$perm_seed <- seed

  dataset
}
