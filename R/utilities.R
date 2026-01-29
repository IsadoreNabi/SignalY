#' @title Utility Functions for SignalY
#' @description Internal utility functions for data validation, transformation,
#'   and common operations used across the SignalY package.
#' @name utilities
#' @keywords internal
NULL


#' Validate Input Data Structure
#'
#' @description
#' Validates that input data meets the requirements for SignalY analysis.
#' Performs comprehensive checks on data types, dimensions, missing values,
#' and numeric properties.
#'
#' @param y Numeric vector representing the target signal. Must be a numeric
#'   vector with no infinite values. Missing values (NA) are handled according
#'   to the \code{na_action} parameter.
#' @param X Optional matrix or data frame of candidate predictors. If provided,
#'   must have the same number of rows as length(y).
#' @param time_index Optional vector of time indices. If NULL, sequential
#'   integers will be used.
#' @param na_action Character string specifying how to handle missing values.
#'   Options are "fail" (stop with error), "omit" (remove observations with NA),
#'   or "interpolate" (linear interpolation).
#' @param min_obs Minimum number of observations required. Default is 10.
#' @param verbose Logical indicating whether to print diagnostic messages.
#'
#' @return A list with validated and potentially transformed data:
#' \describe{
#'   \item{y}{Validated numeric vector of target signal}
#'   \item{X}{Validated matrix of predictors (or NULL)}
#'   \item{time_index}{Vector of time indices}
#'   \item{n_obs}{Number of observations}
#'   \item{n_vars}{Number of predictor variables (or 0)}
#'   \item{var_names}{Names of predictor variables (or NULL)}
#'   \item{na_indices}{Indices of removed observations (if any)}
#' }
#'
#' @details
#' This function implements defensive programming principles to ensure data
#' integrity before computationally intensive analyses. The validation process
#' includes:
#'
#' \enumerate{
#'   \item Type checking: Ensures y is numeric and X (if provided) is numeric
#'     matrix/data.frame
#'   \item Dimension checking: Verifies compatible dimensions between y and X
#'   \item Missing value handling: Processes NA values according to specified
#'     action
#'   \item Finiteness checking: Removes or flags infinite values
#'   \item Minimum sample size: Ensures sufficient observations for analysis
#' }
#'
#' @keywords internal
validate_input <- function(y,
                           X = NULL,
                           time_index = NULL,
                           na_action = c("fail", "omit", "interpolate"),
                           min_obs = 10,
                           verbose = FALSE) {

  na_action <- match.arg(na_action)

  if (!is.numeric(y)) {
    stop("Target signal 'y' must be a numeric vector.", call. = FALSE)
  }

  y <- as.numeric(y)
  n_obs <- length(y)

  if (n_obs < min_obs) {
    stop(sprintf(
      "Insufficient observations: %d provided, minimum %d required.",
      n_obs, min_obs
    ), call. = FALSE)
  }

  if (any(!is.finite(y) & !is.na(y))) {
    inf_idx <- which(!is.finite(y) & !is.na(y))
    stop(sprintf(
      "Target signal contains infinite values at indices: %s",
      paste(head(inf_idx, 5), collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.null(X)) {
    if (!is.matrix(X) && !is.data.frame(X)) {
      stop("Predictor matrix 'X' must be a matrix or data.frame.", call. = FALSE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("Predictor matrix 'X' must contain only numeric values.", call. = FALSE)
    }
    if (nrow(X) != n_obs) {
      stop(sprintf(
        "Dimension mismatch: y has %d observations, X has %d rows.",
        n_obs, nrow(X)
      ), call. = FALSE)
    }
    var_names <- colnames(X)
    if (is.null(var_names)) {
      var_names <- paste0("V", seq_len(ncol(X)))
      colnames(X) <- var_names
    }
    n_vars <- ncol(X)
  } else {
    n_vars <- 0
    var_names <- NULL
  }

  if (is.null(time_index)) {
    time_index <- seq_len(n_obs)
  } else {
    if (length(time_index) != n_obs) {
      stop(sprintf(
        "Time index length (%d) must match number of observations (%d).",
        length(time_index), n_obs
      ), call. = FALSE)
    }
  }

  na_y <- is.na(y)
  na_X <- if (!is.null(X)) apply(X, 1, function(row) any(is.na(row))) else rep(FALSE, n_obs)
  na_any <- na_y | na_X

  na_indices <- NULL

  if (any(na_any)) {
    if (na_action == "fail") {
      stop(sprintf(
        "Missing values detected at %d observation(s). Set na_action to 'omit' or 'interpolate'.",
        sum(na_any)
      ), call. = FALSE)
    } else if (na_action == "omit") {
      na_indices <- which(na_any)
      y <- y[!na_any]
      if (!is.null(X)) X <- X[!na_any, , drop = FALSE]
      time_index <- time_index[!na_any]
      n_obs <- length(y)
      if (verbose) {
        message(sprintf("Removed %d observations with missing values.", length(na_indices)))
      }
    } else if (na_action == "interpolate") {
      if (any(na_y)) {
        y <- interpolate_na(y)
      }
      if (!is.null(X) && any(na_X)) {
        for (j in seq_len(ncol(X))) {
          if (any(is.na(X[, j]))) {
            X[, j] <- interpolate_na(X[, j])
          }
        }
      }
      if (verbose) {
        message("Interpolated missing values using linear interpolation.")
      }
    }
  }

  if (n_obs < min_obs) {
    stop(sprintf(
      "After NA handling, only %d observations remain (minimum %d required).",
      n_obs, min_obs
    ), call. = FALSE)
  }

  list(
    y = y,
    X = X,
    time_index = time_index,
    n_obs = n_obs,
    n_vars = n_vars,
    var_names = var_names,
    na_indices = na_indices
  )
}


#' Linear Interpolation for Missing Values
#'
#' @description
#' Performs linear interpolation to fill missing values in a numeric vector.
#' Handles edge cases where NA values occur at the beginning or end of the
#' vector by using nearest non-NA values.
#'
#' @param x Numeric vector potentially containing NA values.
#'
#' @return Numeric vector with NA values replaced by interpolated values.
#'
#' @details
#' The interpolation uses \code{stats::approx} with linear method. For
#' leading NAs, the first non-NA value is used. For trailing NAs, the last
#' non-NA value is used. This ensures no NA values remain in the output.
#'
#' @keywords internal
interpolate_na <- function(x) {
  if (!any(is.na(x))) return(x)

  n <- length(x)
  idx <- seq_len(n)
  non_na <- !is.na(x)

  if (sum(non_na) < 2) {
    stop("Cannot interpolate: fewer than 2 non-NA values.", call. = FALSE)
  }

  interpolated <- stats::approx(
    x = idx[non_na],
    y = x[non_na],
    xout = idx,
    method = "linear",
    rule = 2
  )$y

  interpolated
}


#' Check Stationarity of AR Coefficients
#'
#' @description
#' Verifies whether autoregressive coefficients satisfy stationarity conditions
#' by checking that all roots of the characteristic polynomial lie outside
#' the unit circle.
#'
#' @param phi Numeric vector of AR coefficients (phi_1, phi_2, ..., phi_p).
#'
#' @return A list with:
#' \describe{
#'   \item{is_stationary}{Logical indicating whether the process is stationary}
#'   \item{roots}{Complex roots of the characteristic polynomial}
#'   \item{moduli}{Moduli of the roots}
#' }
#'
#' @details
#' For an AR(p) process, stationarity requires that all roots of the
#' characteristic polynomial \eqn{1 - \phi_1 z - \phi_2 z^2 - ... - \phi_p z^p}
#' have modulus greater than 1. This is equivalent to requiring that all roots
#' of the polynomial lie outside the unit circle in the complex plane.
#'
#' The function constructs the companion matrix and computes its eigenvalues,
#' which are the inverses of the characteristic polynomial roots.
#'
#' @references
#' Hamilton, J. D. (1994). Time Series Analysis. Princeton University Press.
#' Chapter 1.
#'
#' @keywords internal
check_stationarity <- function(phi) {
  p <- length(phi)

  if (p == 0) {
    return(list(is_stationary = TRUE, roots = numeric(0), moduli = numeric(0)))
  }

  coeffs <- c(1, -phi)
  roots <- polyroot(coeffs)
  moduli <- Mod(roots)

  list(
    is_stationary = all(moduli > 1),
    roots = roots,
    moduli = moduli
  )
}


#' Compute Shannon Entropy
#'
#' @description
#' Calculates the Shannon entropy of a probability distribution or, when applied
#' to loadings, the entropy of the squared normalized loadings. High entropy
#' indicates diffuse/uniform distribution (systemic noise), while low entropy
#' indicates concentrated structure.
#'
#' @param x Numeric vector. Will be squared and normalized to form a probability
#'   distribution.
#' @param base Base of the logarithm. Default is 2 (bits).
#' @param normalize Logical. If TRUE, returns normalized entropy (0 to 1 scale).
#'
#' @return Numeric scalar representing entropy value.
#'
#' @details
#' The Shannon entropy is defined as:
#' \deqn{H(p) = -\sum_{i} p_i \log(p_i)}
#'
#' where \eqn{p_i} are the probabilities. For factor loadings, we use squared
#' normalized loadings as the probability distribution:
#' \deqn{p_i = \lambda_i^2 / \sum_j \lambda_j^2}
#'
#' This measures the concentration of explanatory power across variables.
#' Maximum entropy occurs when all loadings are equal (diffuse structure);
#' minimum entropy occurs when a single variable dominates (concentrated
#' structure).
#'
#' @section Interpretation in Signal Analysis:
#' In the context of latent structure extraction:
#' \itemize{
#'   \item **High entropy (near maximum)**: Suggests "maximum entropy systemic
#'     stochasticity" - the component captures diffuse, undifferentiated
#'     movement across all variables (akin to Brownian motion).
#'   \item **Low entropy**: Suggests "differentiated latent structure" - the
#'     component is driven by a subset of variables, indicating meaningful
#'     structural relationships.
#' }
#'
#' @examples
#' uniform_loadings <- rep(1, 10)
#' compute_entropy(uniform_loadings, normalize = TRUE)
#'
#' concentrated_loadings <- c(10, rep(0.1, 9))
#' compute_entropy(concentrated_loadings, normalize = TRUE)
#'
#' @export
compute_entropy <- function(x, base = 2, normalize = FALSE) {
  x <- as.numeric(x)
  x_sq <- x^2
  total <- sum(x_sq)

  if (total == 0) {
    warning("All values are zero; returning NA for entropy.")
    return(NA_real_)
  }

  p <- x_sq / total
  p <- p[p > 0]

  entropy <- -sum(p * log(p, base = base))

  if (normalize) {
    max_entropy <- log(length(x), base = base)
    if (max_entropy > 0) {
      entropy <- entropy / max_entropy
    }
  }

  entropy
}


#' Compute Partial R-squared
#'
#' @description
#' Calculates the partial R-squared for a specific predictor or block of
#' predictors, measuring their unique contribution to explaining variance
#' in the outcome after controlling for other predictors.
#'
#' @param y Numeric vector of the outcome variable.
#' @param X_interest Matrix of predictors of interest.
#' @param X_control Matrix of control predictors. Can be NULL for simple
#'   R-squared.
#' @param weights Optional observation weights.
#'
#' @return Numeric scalar representing partial R-squared (0 to 1).
#'
#' @details
#' Partial R-squared is computed as:
#' \deqn{R^2_{partial} = (SSE_{reduced} - SSE_{full}) / SSE_{reduced}}
#'
#' where SSE is the sum of squared errors. This measures the proportional
#' reduction in unexplained variance achieved by adding the predictors of
#' interest to a model that already contains the control predictors.
#'
#' @references
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). Applied Multiple
#' Regression/Correlation Analysis for the Behavioral Sciences (3rd ed.).
#' Lawrence Erlbaum Associates. Chapter 3.
#'
#' @keywords internal
compute_partial_r2 <- function(y, X_interest, X_control = NULL, weights = NULL) {
  n <- length(y)

  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  if (is.null(X_control)) {
    X_full <- cbind(1, X_interest)
  } else {
    X_reduced <- cbind(1, X_control)
    X_full <- cbind(X_reduced, X_interest)
  }

  fit_full <- stats::lm.wfit(X_full, y, w = weights)
  sse_full <- sum(weights * fit_full$residuals^2)

  if (is.null(X_control)) {
    sse_reduced <- sum(weights * (y - stats::weighted.mean(y, weights))^2)
  } else {
    fit_reduced <- stats::lm.wfit(X_reduced, y, w = weights)
    sse_reduced <- sum(weights * fit_reduced$residuals^2)
  }

  if (sse_reduced == 0) return(NA_real_)

  partial_r2 <- (sse_reduced - sse_full) / sse_reduced
  max(0, min(1, partial_r2))
}


#' Apply Function to Matrix Columns
#'
#' @description
#' Applies a univariate filtering or transformation function to each column
#' of a matrix and returns a consolidated data frame. This utility enables
#' batch processing of panel data where each column represents a different
#' variable or series.
#'
#' @param X Matrix or data frame where each column is a series to process.
#' @param FUN Function to apply to each column. Must accept a numeric vector
#'   and return either a numeric vector of the same length or a list with
#'   a named element (specified by \code{extract}).
#' @param extract Character string specifying which element to extract from
#'   the function output if it returns a list. Default is NULL (use raw output).
#' @param ... Additional arguments passed to FUN.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A data frame with the same number of rows as X, containing the
#'   processed output for each column.
#'
#' @examples
#' X <- matrix(rnorm(200), ncol = 4)
#' colnames(X) <- c("A", "B", "C", "D")
#' result <- apply_to_columns(X, function(x) cumsum(x))
#'
#' @export
apply_to_columns <- function(X, FUN, extract = NULL, ..., verbose = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  col_names <- colnames(X)
  if (is.null(col_names)) col_names <- paste0("V", seq_len(p))

  results <- vector("list", p)

  for (j in seq_len(p)) {
    if (verbose) message(sprintf("Processing column %d/%d: %s", j, p, col_names[j]))

    out <- FUN(X[, j], ...)

    if (!is.null(extract) && is.list(out)) {
      out <- out[[extract]]
    }

    if (length(out) != n) {
      stop(sprintf(
        "Output length (%d) does not match input length (%d) for column '%s'.",
        length(out), n, col_names[j]
      ), call. = FALSE)
    }

    results[[j]] <- out
  }

  result_df <- as.data.frame(results)
  colnames(result_df) <- col_names
  result_df
}


#' Create Block Bootstrap Samples
#'
#' @description
#' Generates block bootstrap samples for time series data, preserving temporal
#' dependence structure. Uses non-overlapping blocks of specified length with
#' random block selection with replacement.
#'
#' @param X Matrix where rows are time points and columns are variables.
#' @param n_boot Number of bootstrap samples to generate.
#' @param block_length Length of each block. If NULL, defaults to
#'   \code{ceiling(sqrt(nrow(X)))}.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A 3-dimensional array of dimension (nrow(X), ncol(X), n_boot)
#'   containing the bootstrap samples.
#'
#' @details
#' Block bootstrapping is essential for time series data to preserve the
#' autocorrelation structure within blocks. The procedure:
#' \enumerate{
#'   \item Divides the original series into blocks of length L
#'   \item Randomly samples blocks with replacement
#'   \item Concatenates sampled blocks to form each bootstrap sample
#' }
#'
#' The default block length of \eqn{\sqrt{T}} balances capturing dependence
#' structure against having sufficient blocks for resampling variation.
#'
#' @section Synchronization:
#' For multivariate series, blocks are sampled synchronously across all
#' variables to preserve cross-sectional dependence. This is critical when
#' bootstrap samples are used to construct confidence intervals for statistics
#' that depend on the joint distribution.
#'
#' @references
#' Kunsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary
#' Observations. The Annals of Statistics, 17(3), 1217-1241.
#'
#' Politis, D. N., & Romano, J. P. (1994). The Stationary Bootstrap. Journal
#' of the American Statistical Association, 89(428), 1303-1313.
#'
#' @keywords internal
block_bootstrap <- function(X, n_boot, block_length = NULL, seed = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(block_length)) {
    block_length <- ceiling(sqrt(n))
  }

  n_blocks <- ceiling(n / block_length)

  if (!is.null(seed)) set.seed(seed)

  boot_array <- array(NA_real_, dim = c(n, p, n_boot))

  for (b in seq_len(n_boot)) {
    selected_blocks <- sample(seq_len(n_blocks), n_blocks, replace = TRUE)

    boot_indices <- integer(0)
    for (block_id in selected_blocks) {
      start_idx <- (block_id - 1) * block_length + 1
      end_idx <- min(block_id * block_length, n)
      boot_indices <- c(boot_indices, start_idx:end_idx)
    }

    boot_indices <- boot_indices[seq_len(n)]
    boot_array[, , b] <- X[boot_indices, , drop = FALSE]
  }

  boot_array
}


#' Format Numeric Values for Display
#'
#' @description
#' Formats numeric values with appropriate precision and scientific notation
#' handling for display in reports and console output.
#'
#' @param x Numeric vector to format.
#' @param digits Number of significant digits. Default is 3.
#' @param scientific Logical. If TRUE, use scientific notation for very large
#'   or small values.
#'
#' @return Character vector of formatted numbers.
#'
#' @keywords internal
format_numeric <- function(x, digits = 3, scientific = FALSE) {
  if (scientific) {
    format(x, digits = digits, scientific = TRUE)
  } else {
    format(round(x, digits), nsmall = digits, scientific = FALSE)
  }
}


#' Safe Division with Zero Handling
#'
#' @description
#' Performs division with protection against division by zero, returning
#' a specified value (default NA) when the denominator is zero or very small.
#'
#' @param numerator Numeric vector of numerators.
#' @param denominator Numeric vector of denominators.
#' @param zero_value Value to return when denominator is effectively zero.
#'   Default is NA.
#' @param tol Tolerance for considering denominator as zero. Default is
#'   .Machine$double.eps.
#'
#' @return Numeric vector of quotients.
#'
#' @keywords internal
safe_divide <- function(numerator, denominator, zero_value = NA_real_,
                        tol = .Machine$double.eps) {
  result <- numerator / denominator
  result[abs(denominator) < tol] <- zero_value
  result
}
