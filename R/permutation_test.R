# =============================================================================
# permutation_test.R
# Permutation-based significance testing for CV approaches
#
# Provides a uniform way to determine significance
# for any CV method. The observed test statistic is compared against a
# null distribution built by running the same method on permuted data
# (Y rows shuffled).
# =============================================================================

source("R/generate_data.R")  # for permute_dataset


#' Compute a permutation-based p-value for a CV approach
#'
#' Runs the full CV method on N_perm permuted versions of the dataset
#' and computes the empirical p-value: the fraction of permuted test
#' statistics >= the observed statistic.
#'
#' @param approach       One of: "naive_cv", "repeated_cv", "nested_cv", "caspoc"
#' @param X, Y           Data matrices (original, unpermuted)
#' @param observed_stat  The test statistic from the real data run
#' @param n_perm         Number of permutations (default 100)
#' @param seed           Base seed for permutation reproducibility
#' @param ...            Additional arguments passed to run_cv_approach()
#'
#' @return A list with:
#'   \item{perm_pvalue}{Empirical p-value}
#'   \item{observed_stat}{The observed test statistic (passed through)}
#'   \item{null_distribution}{Numeric vector of permuted test statistics}
#'   \item{n_perm}{Number of permutations used}
compute_permutation_pvalue <- function(approach, X, Y,
                                       observed_stat,
                                       n_perm = 100,
                                       seed = 1,
                                       ...) {

  null_stats <- numeric(n_perm)

  for (p in seq_len(n_perm)) {
    # Permute Y
    set.seed(seed * 1000 + p)  # unique seed per permutation
    perm_idx <- sample(nrow(Y))
    Y_perm <- Y[perm_idx, , drop = FALSE]

    # Run the full method on permuted data
    perm_result <- tryCatch(
      run_cv_approach(
        approach = approach,
        X = X, Y = Y_perm,
        seed = seed,  # same CV seed so fold structure is comparable
        ...
      ),
      error = function(e) NULL
    )

    null_stats[p] <- if (!is.null(perm_result)) perm_result$observed_stat else NA
  }

  # Empirical p-value: fraction of null stats >= observed
  # Add 1 to numerator and denominator for conservative estimate
  # (Phipson & Smyth, 2010)
  n_extreme <- sum(null_stats >= observed_stat, na.rm = TRUE)
  n_valid   <- sum(!is.na(null_stats))
  perm_pvalue <- (n_extreme + 1) / (n_valid + 1)

  list(
    perm_pvalue       = perm_pvalue,
    observed_stat     = observed_stat,
    null_distribution = null_stats,
    n_perm            = n_perm
  )
}
