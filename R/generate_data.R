# =============================================================================
# generate_data.R
# Functions to generate simulated datasets for CASPOC benchmarking
# =============================================================================

library(MASS)  # mvrnorm


# --- Null data (Type I error) ------------------------------------------------

#' Generate null data with no association between X and Y
#'
#' X and Y are drawn independently from multivariate normal distributions.
#' Any association found by a method on this data is a false positive.
#'
#' @param n      Number of samples
#' @param p      Number of X features
#' @param q      Number of Y features
#' @param cor_x  Within-block correlation in X (0 = independent features)
#' @param cor_y  Within-block correlation in Y (0 = independent features)
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
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Build covariance matrices
  # TODO: Decide on block structure vs compound symmetry vs Toeplitz

  # For now: compound symmetry (all pairwise correlations equal)
  Sigma_x <- matrix(cor_x, nrow = p, ncol = p)
  diag(Sigma_x) <- 1

  Sigma_y <- matrix(cor_y, nrow = q, ncol = q)
  diag(Sigma_y) <- 1

  # Draw X and Y independently
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_x)
  Y <- mvrnorm(n = n, mu = rep(0, q), Sigma = Sigma_y)

  colnames(X) <- paste0("X", seq_len(p))
  colnames(Y) <- paste0("Y", seq_len(q))

  list(
    X = X,
    Y = Y,
    true_signal = FALSE,
    params = list(
      n = n, p = p, q = q,
      cor_x = cor_x, cor_y = cor_y,
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
#' @param signal_strength Variance of the latent signal relative to noise (SNR)
#' @param cor_x          Background correlation among X features
#' @param cor_y          Background correlation among Y features
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
                                 seed = NULL) {
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
    offset_x <- (comp - 1) * n_relevant_x
    offset_y <- (comp - 1) * n_relevant_y

    idx_x <- (offset_x + 1):min(offset_x + n_relevant_x, p)
    idx_y <- (offset_y + 1):min(offset_y + n_relevant_y, q)

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
  # Noise part: multivariate normal with optional correlation
  signal_x <- sqrt(signal_strength) * Z %*% t(loadings_x)
  signal_y <- sqrt(signal_strength) * Z %*% t(loadings_y)

  # Background noise
  Sigma_x <- matrix(cor_x, nrow = p, ncol = p)
  diag(Sigma_x) <- 1
  Sigma_y <- matrix(cor_y, nrow = q, ncol = q)
  diag(Sigma_y) <- 1

  noise_x <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma_x)
  noise_y <- mvrnorm(n = n, mu = rep(0, q), Sigma = Sigma_y)

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
