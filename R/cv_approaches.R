# =============================================================================
# cv_approaches.R
# Uniform wrappers for each cross-validation approach
#
# All wrappers follow the same interface:
#   Input:  X, Y, hyperparameter grid, CV settings
#   Output: A standardised list (see return spec below)
#
# Significance is determined externally via permutation testing
# (see permutation_test.R), NOT by fold-level p-values.
# =============================================================================

library(mixOmics)
library(caret)
library(caspoc)
library(dplyr)
library(tibble)


# =============================================================================
# Standardised output structure
# =============================================================================
#
# Every wrapper returns a list with:
#   $approach        : character — name of the approach
#   $selected_keepX  : numeric — chosen keepX
#   $selected_keepY  : numeric — chosen keepY
#   $observed_stat   : numeric — the method's reported test statistic (correlation)
#                      For naive/repeated CV: best CV correlation (biased)
#                      For nested CV: median outer test correlation (unbiased)
#                      For CASPOC: median test fold correlation (unbiased)
#   $all_results     : data.frame — full results across all HPs and folds
#   $runtime         : numeric — elapsed seconds
#   $details         : list — any approach-specific extras


# =============================================================================
# 1. Naïve k-fold CV
# =============================================================================
#
# Standard k-fold: same folds used for selecting hyperparameters AND reporting
# performance. This is the biased baseline.

run_naive_cv <- function(X, Y,
                         ncomp = 1,
                         num_folds = 10,
                         keepX_options = NULL,
                         keepY_options = NULL,
                         seed = 1,
                         ...) {

  start_time <- proc.time()

  if (is.null(keepX_options)) keepX_options <- make_default_grid(ncol(X))
  if (is.null(keepY_options)) keepY_options <- make_default_grid(ncol(Y))

  set.seed(seed)
  folds <- createFolds(seq_len(nrow(X)), k = num_folds, list = TRUE)

  results <- data.frame()

  for (kx in keepX_options) {
    for (ky in keepY_options) {
      fold_cors <- numeric(num_folds)

      for (i in seq_len(num_folds)) {
        test_idx <- folds[[i]]
        train_idx <- setdiff(seq_len(nrow(X)), test_idx)

        train_X <- scale(X[train_idx, , drop = FALSE])
        train_Y <- scale(Y[train_idx, , drop = FALSE])

        test_X <- scale(X[test_idx, , drop = FALSE],
                        center = attr(train_X, "scaled:center"),
                        scale  = attr(train_X, "scaled:scale"))
        test_Y <- scale(Y[test_idx, , drop = FALSE],
                        center = attr(train_Y, "scaled:center"),
                        scale  = attr(train_Y, "scaled:scale"))

        fit <- tryCatch(
          spls(train_X, train_Y,
               ncomp = ncomp,
               keepX = rep(kx, ncomp),
               keepY = rep(ky, ncomp),
               mode = "regression", scale = TRUE),
          error = function(e) NULL
        )
        if (is.null(fit)) { fold_cors[i] <- NA; next }

        scores_x <- test_X %*% fit$loadings$X[, ncomp, drop = FALSE]
        scores_y <- test_Y %*% fit$loadings$Y[, ncomp, drop = FALSE]

        fold_cors[i] <- cor(scores_x[, 1], scores_y[, 1], method = "spearman")
      }

      results <- rbind(results, data.frame(
        keepX = kx, keepY = ky,
        mean_cor = mean(fold_cors, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  best_idx <- which.max(results$mean_cor)
  elapsed <- (proc.time() - start_time)["elapsed"]

  list(
    approach        = "naive_cv",
    selected_keepX  = results$keepX[best_idx],
    selected_keepY  = results$keepY[best_idx],
    observed_stat   = results$mean_cor[best_idx],
    all_results     = results,
    runtime         = elapsed,
    details         = list(folds = folds)
  )
}


# =============================================================================
# 2. Repeated k-fold CV with median selection
# =============================================================================
#
# Multiple repeats of standard k-fold, select HPs by median correlation.
# Still uses the same folds for tuning and evaluation — biased, but
# the median reduces variance.

run_repeated_cv <- function(X, Y,
                            ncomp = 1,
                            num_folds = 10,
                            num_repeats = 11,
                            keepX_options = NULL,
                            keepY_options = NULL,
                            seed = 1,
                            ...) {

  start_time <- proc.time()

  if (is.null(keepX_options)) keepX_options <- make_default_grid(ncol(X))
  if (is.null(keepY_options)) keepY_options <- make_default_grid(ncol(Y))

  all_results <- data.frame()

  for (rep in seq_len(num_repeats)) {
    set.seed(seed + rep)
    folds <- createFolds(seq_len(nrow(X)), k = num_folds, list = TRUE)

    for (kx in keepX_options) {
      for (ky in keepY_options) {
        fold_cors <- numeric(num_folds)

        for (i in seq_len(num_folds)) {
          test_idx <- folds[[i]]
          train_idx <- setdiff(seq_len(nrow(X)), test_idx)

          train_X <- scale(X[train_idx, , drop = FALSE])
          train_Y <- scale(Y[train_idx, , drop = FALSE])

          test_X <- scale(X[test_idx, , drop = FALSE],
                          center = attr(train_X, "scaled:center"),
                          scale  = attr(train_X, "scaled:scale"))
          test_Y <- scale(Y[test_idx, , drop = FALSE],
                          center = attr(train_Y, "scaled:center"),
                          scale  = attr(train_Y, "scaled:scale"))

          fit <- tryCatch(
            spls(train_X, train_Y,
                 ncomp = ncomp,
                 keepX = rep(kx, ncomp),
                 keepY = rep(ky, ncomp),
                 mode = "regression", scale = TRUE),
            error = function(e) NULL
          )
          if (is.null(fit)) { fold_cors[i] <- NA; next }

          scores_x <- test_X %*% fit$loadings$X[, ncomp, drop = FALSE]
          scores_y <- test_Y %*% fit$loadings$Y[, ncomp, drop = FALSE]

          fold_cors[i] <- cor(scores_x[, 1], scores_y[, 1], method = "spearman")
        }

        all_results <- rbind(all_results, data.frame(
          Repeat = rep, keepX = kx, keepY = ky,
          mean_cor = mean(fold_cors, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  summary_results <- all_results %>%
    group_by(keepX, keepY) %>%
    summarise(median_cor = median(mean_cor, na.rm = TRUE), .groups = "drop")

  best_idx <- which.max(summary_results$median_cor)
  elapsed <- (proc.time() - start_time)["elapsed"]

  list(
    approach        = "repeated_cv",
    selected_keepX  = summary_results$keepX[best_idx],
    selected_keepY  = summary_results$keepY[best_idx],
    observed_stat   = summary_results$median_cor[best_idx],
    all_results     = all_results,
    runtime         = elapsed,
    details         = list(summary = summary_results)
  )
}


# =============================================================================
# 3. Nested CV (repeated)
# =============================================================================
#
# Outer loop: hold out test fold (unbiased evaluation)
# Inner loop: tune hyperparameters on remaining data
# Then retrain on all outer training data with best HPs, evaluate on outer test.
# Repeated R times with different fold partitions for fair comparison with CASPOC.

run_nested_cv <- function(X, Y,
                          ncomp = 1,
                          num_folds = 10,
                          num_folds_inner = 5,
                          num_repeats = 11,
                          keepX_options = NULL,
                          keepY_options = NULL,
                          seed = 1,
                          ...) {

  num_folds_outer <- num_folds

  start_time <- proc.time()

  if (is.null(keepX_options)) keepX_options <- make_default_grid(ncol(X))
  if (is.null(keepY_options)) keepY_options <- make_default_grid(ncol(Y))

  all_outer_results <- data.frame()

  for (rep in seq_len(num_repeats)) {
    set.seed(seed + rep)
    outer_folds <- createFolds(seq_len(nrow(X)), k = num_folds_outer, list = TRUE)

    for (i in seq_len(num_folds_outer)) {
      outer_test_idx  <- outer_folds[[i]]
      outer_train_idx <- setdiff(seq_len(nrow(X)), outer_test_idx)

      X_outer_train <- X[outer_train_idx, , drop = FALSE]
      Y_outer_train <- Y[outer_train_idx, , drop = FALSE]

      # --- Inner loop: tune HPs ---
      inner_folds <- createFolds(seq_len(nrow(X_outer_train)),
                                 k = num_folds_inner, list = TRUE)
      inner_results <- data.frame()

      for (kx in keepX_options) {
        for (ky in keepY_options) {
          inner_cors <- numeric(num_folds_inner)

          for (j in seq_len(num_folds_inner)) {
            inner_test_idx  <- inner_folds[[j]]
            inner_train_idx <- setdiff(seq_len(nrow(X_outer_train)), inner_test_idx)

            train_X <- scale(X_outer_train[inner_train_idx, , drop = FALSE])
            train_Y <- scale(Y_outer_train[inner_train_idx, , drop = FALSE])

            test_X <- scale(X_outer_train[inner_test_idx, , drop = FALSE],
                            center = attr(train_X, "scaled:center"),
                            scale  = attr(train_X, "scaled:scale"))
            test_Y <- scale(Y_outer_train[inner_test_idx, , drop = FALSE],
                            center = attr(train_Y, "scaled:center"),
                            scale  = attr(train_Y, "scaled:scale"))

            fit <- tryCatch(
              spls(train_X, train_Y,
                   ncomp = ncomp,
                   keepX = rep(kx, ncomp),
                   keepY = rep(ky, ncomp),
                   mode = "regression", scale = TRUE),
              error = function(e) NULL
            )
            if (is.null(fit)) { inner_cors[j] <- NA; next }

            scores_x <- test_X %*% fit$loadings$X[, ncomp, drop = FALSE]
            scores_y <- test_Y %*% fit$loadings$Y[, ncomp, drop = FALSE]

            inner_cors[j] <- cor(scores_x[, 1], scores_y[, 1], method = "spearman")
          }

          inner_results <- rbind(inner_results, data.frame(
            keepX = kx, keepY = ky,
            mean_inner_cor = mean(inner_cors, na.rm = TRUE),
            stringsAsFactors = FALSE
          ))
        }
      }

      # Select best HPs from inner loop
      best_inner <- which.max(inner_results$mean_inner_cor)
      best_kx <- inner_results$keepX[best_inner]
      best_ky <- inner_results$keepY[best_inner]

      # --- Retrain on full outer training set with best HPs ---
      train_X <- scale(X_outer_train)
      train_Y <- scale(Y_outer_train)

      test_X <- scale(X[outer_test_idx, , drop = FALSE],
                      center = attr(train_X, "scaled:center"),
                      scale  = attr(train_X, "scaled:scale"))
      test_Y <- scale(Y[outer_test_idx, , drop = FALSE],
                      center = attr(train_Y, "scaled:center"),
                      scale  = attr(train_Y, "scaled:scale"))

      fit_final <- tryCatch(
        spls(train_X, train_Y,
             ncomp = ncomp,
             keepX = rep(best_kx, ncomp),
             keepY = rep(best_ky, ncomp),
             mode = "regression", scale = TRUE),
        error = function(e) NULL
      )

      if (!is.null(fit_final)) {
        scores_x <- test_X %*% fit_final$loadings$X[, ncomp, drop = FALSE]
        scores_y <- test_Y %*% fit_final$loadings$Y[, ncomp, drop = FALSE]
        outer_cor <- cor(scores_x[, 1], scores_y[, 1], method = "spearman")
      } else {
        outer_cor <- NA
      }

      all_outer_results <- rbind(all_outer_results, data.frame(
        Repeat     = rep,
        outer_fold = i,
        best_keepX = best_kx,
        best_keepY = best_ky,
        tune_cor   = inner_results$mean_inner_cor[best_inner],
        test_cor   = outer_cor,
        stringsAsFactors = FALSE
      ))
    }
  }

  elapsed <- (proc.time() - start_time)["elapsed"]

  list(
    approach        = "nested_cv",
    selected_keepX  = median(all_outer_results$best_keepX),
    selected_keepY  = median(all_outer_results$best_keepY),
    observed_stat   = median(all_outer_results$test_cor, na.rm = TRUE),
    all_results     = all_outer_results,
    runtime         = elapsed,
    details         = list()
  )
}


# =============================================================================
# 4. CASPOC
# =============================================================================
#
# Double-split k-fold CV with separate tuning and test folds.
# Wraps the caspoc package.

run_caspoc <- function(X, Y,
                       ncomp = 1,
                       num_folds = 10,
                       num_repeats = 11,
                       keepX_options = NULL,
                       keepY_options = NULL,
                       seed = 1,
                       ...) {

  start_time <- proc.time()

  if (is.null(keepX_options)) keepX_options <- make_default_grid(ncol(X))
  if (is.null(keepY_options)) keepY_options <- make_default_grid(ncol(Y))

  res <- CASPOC(
    X = X, Y = Y,
    ncomp = ncomp,
    numRepeats = num_repeats,
    numFolds = num_folds,
    keepX_options = keepX_options,
    keepY_options = keepY_options,
    base_seed = seed
  )

  # Select best HPs by median tuning correlation
  tune_summary <- res$results_tune_df %>%
    filter(Component == ncomp) %>%
    group_by(KeepX, KeepY) %>%
    summarise(
      median_tune_cor = median(Correlation, na.rm = TRUE),
      .groups = "drop"
    )

  best_idx <- which.max(tune_summary$median_tune_cor)
  best_kx <- tune_summary$KeepX[best_idx]
  best_ky <- tune_summary$KeepY[best_idx]

  # Get test correlation for selected HPs
  test_for_best <- res$results_test_df %>%
    filter(Component == ncomp, KeepX == best_kx, KeepY == best_ky)

  elapsed <- (proc.time() - start_time)["elapsed"]

  list(
    approach        = "caspoc",
    selected_keepX  = best_kx,
    selected_keepY  = best_ky,
    observed_stat   = median(test_for_best$Correlation, na.rm = TRUE),
    all_results     = list(tune = res$results_tune_df, test = res$results_test_df),
    runtime         = elapsed,
    details         = list(
      caspoc_result = res,
      tune_summary  = tune_summary
    )
  )
}


# =============================================================================
# Helper: default hyperparameter grid
# =============================================================================

make_default_grid <- function(num_features) {
  if (num_features <= 10) {
    seq(1, num_features)
  } else {
    step <- ceiling(num_features / 10)
    grid <- seq(1, num_features, by = step)
    if (!(num_features %in% grid)) grid <- c(grid, num_features)
    grid
  }
}


# =============================================================================
# Dispatcher: run any approach by name
# =============================================================================

#' Run a cross-validation approach by name
#'
#' @param approach  One of: "naive_cv", "repeated_cv", "nested_cv", "caspoc"
#' @param X, Y      Data matrices
#' @param ...       Additional arguments passed to the specific wrapper
#'
#' @return Standardised result list (see top of this file)
run_cv_approach <- function(approach, X, Y, ...) {
  fn <- switch(approach,
    naive_cv    = run_naive_cv,
    repeated_cv = run_repeated_cv,
    nested_cv   = run_nested_cv,
    caspoc      = run_caspoc,
    stop("Unknown approach: ", approach)
  )
  fn(X, Y, ...)
}
