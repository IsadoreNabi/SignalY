#' @title Principal Component Analysis and Dynamic Factor Models
#' @description
#' Implements dimensionality reduction techniques for panel data including
#' PCA with bootstrap significance testing and Dynamic Factor Models (DFM)
#' for extracting common latent factors.
#'
#' @name pca_dfm
NULL


#' Principal Component Analysis with Bootstrap Significance Testing
#'
#' @description
#' Performs PCA on panel data with bootstrap-based significance testing for
#' factor loadings. Identifies which variables load significantly on each
#' principal component using a null distribution constructed via block
#' bootstrapping.
#'
#' @param X Matrix or data frame where rows are observations (time points)
#'   and columns are variables.
#' @param n_components Number of principal components to extract. If NULL,
#'   determined by eigenvalue threshold or explained variance.
#' @param center Logical. Center variables before PCA. Default TRUE.
#' @param scale Logical. Scale variables to unit variance. Default TRUE.
#' @param n_boot Number of bootstrap replications for significance testing.
#'   Default 200.
#' @param block_length Block length for block bootstrap. If NULL, defaults
#'   to \code{ceiling(sqrt(nrow(X)))}.
#' @param alpha Significance level for loading tests. Default 0.05.
#' @param use_fdr Logical. Apply Benjamini-Hochberg FDR correction.
#'   Default FALSE.
#' @param rotation Character string specifying rotation method: "none",
#'   "varimax", or "oblimin". Default "varimax".
#' @param verbose Logical for progress messages.
#'
#' @return A list of class "signaly_pca" containing:
#' \describe{
#'   \item{loadings}{Matrix of factor loadings (rotated if specified)}
#'   \item{scores}{Matrix of component scores}
#'   \item{eigenvalues}{Vector of eigenvalues}
#'   \item{variance_explained}{Proportion of variance explained by each
#'     component}
#'   \item{cumulative_variance}{Cumulative proportion of variance explained}
#'   \item{significant_loadings}{Matrix of logical values indicating
#'     significance}
#'   \item{p_values}{Matrix of bootstrap p-values for loadings}
#'   \item{thresholds}{Cutoff values for significance by component}
#'   \item{entropy}{Shannon entropy of loadings for each component}
#'   \item{summary_by_component}{Data frame summarizing each component}
#'   \item{assignments}{Data frame mapping variables to their dominant
#'     component}
#' }
#'
#' @details
#' The analysis proceeds in several stages:
#'
#' **1. Standard PCA**: Eigendecomposition of the correlation (if scaled) or
#' covariance matrix to extract principal components.
#'
#' **2. Rotation** (optional): Varimax rotation maximizes the variance of
#' squared loadings within components, producing cleaner simple structure.
#' Oblimin allows correlated factors.
#'
#' **3. Bootstrap Significance Testing**: For each bootstrap replicate:
#' \enumerate{
#'   \item Resample rows using block bootstrap (preserving temporal dependence)
#'   \item Perform PCA on resampled data
#'   \item Apply Procrustes rotation to align with original
#'   \item Record absolute loadings
#' }
#' The empirical p-value for each loading is the proportion of bootstrap
#' loadings exceeding the original in absolute value.
#'
#' **4. Entropy Calculation**: Shannon entropy of squared loadings indicates
#' whether explanatory power is concentrated (low entropy) or diffuse (high
#' entropy). High entropy on PC1 suggests systemic co-movement rather than
#' differentiated structure.
#'
#' @section Interpretation in Signal Analysis:
#' \itemize{
#'   \item **High PC1 entropy**: "Maximum entropy systemic stochasticity" -
#'     the dominant factor captures undifferentiated movement, suggesting
#'     noise rather than latent structure.
#'   \item **Low PC1 entropy**: "Differentiated latent structure" - specific
#'     variables dominate, indicating meaningful groupings.
#'   \item **Significant loadings**: Variables with p < alpha after bootstrap
#'     testing reliably load on that component.
#' }
#'
#' @references
#' Jolliffe, I. T. (2002). Principal Component Analysis (2nd ed.). Springer.
#'
#' Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor
#' analysis. Psychometrika, 23(3), 187-200.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n * p), ncol = p)
#' colnames(X) <- paste0("V", 1:p)
#' result <- pca_bootstrap(X, n_components = 3, n_boot = 50)
#' print(result$summary_by_component)
#'
#' @export
pca_bootstrap <- function(X,
                          n_components = NULL,
                          center = TRUE,
                          scale = TRUE,
                          n_boot = 200,
                          block_length = NULL,
                          alpha = 0.05,
                          use_fdr = FALSE,
                          rotation = c("varimax", "none", "oblimin"),
                          verbose = FALSE) {

  rotation <- match.arg(rotation)

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  var_names <- colnames(X)
  if (is.null(var_names)) {
    var_names <- paste0("V", seq_len(p))
    colnames(X) <- var_names
  }

  if (any(is.na(X))) {
    stop("Matrix contains NA values. Please handle missing data first.",
         call. = FALSE)
  }

  if (is.null(block_length)) {
    block_length <- ceiling(sqrt(n))
  }

  if (verbose) message("Performing initial PCA...")

  pca_result <- stats::prcomp(X, center = center, scale. = scale)

  eigenvalues <- pca_result$sdev^2
  var_explained <- eigenvalues / sum(eigenvalues)
  cumvar <- cumsum(var_explained)

  if (is.null(n_components)) {
    n_components <- max(1, sum(eigenvalues > 1))
    n_components <- min(n_components, which(cumvar >= 0.9)[1])
    n_components <- max(1, min(n_components, p - 1, n - 1))
  }
  n_components <- min(n_components, p, n - 1)

  if (verbose) {
    message(sprintf("Extracting %d components (%.1f%% variance explained)",
                    n_components, 100 * cumvar[n_components]))
  }

  loadings <- pca_result$rotation[, seq_len(n_components), drop = FALSE]
  scores <- pca_result$x[, seq_len(n_components), drop = FALSE]

  if (rotation != "none" && n_components > 1) {
    if (verbose) message(sprintf("Applying %s rotation...", rotation))

    if (rotation == "varimax") {
      rot_result <- stats::varimax(loadings)
      loadings <- rot_result$loadings[, , drop = FALSE]
      class(loadings) <- "matrix"
    } else if (rotation == "oblimin") {
      if (requireNamespace("GPArotation", quietly = TRUE)) {
        rot_result <- GPArotation::oblimin(loadings)
        loadings <- rot_result$loadings
      } else {
        warning("GPArotation package not available. Using varimax instead.")
        rot_result <- stats::varimax(loadings)
        loadings <- rot_result$loadings[, , drop = FALSE]
        class(loadings) <- "matrix"
      }
    }
  }

  if (verbose) message(sprintf("Running %d bootstrap replications...", n_boot))

  boot_loadings <- array(NA_real_, dim = c(p, n_components, n_boot))

  for (b in seq_len(n_boot)) {
    if (verbose && b %% 50 == 0) message(sprintf("  Bootstrap %d/%d", b, n_boot))

    n_blocks <- ceiling(n / block_length)
    selected_blocks <- sample(seq_len(n_blocks), n_blocks, replace = TRUE)

    boot_indices <- integer(0)
    for (block_id in selected_blocks) {
      start_idx <- (block_id - 1) * block_length + 1
      end_idx <- min(block_id * block_length, n)
      boot_indices <- c(boot_indices, start_idx:end_idx)
    }
    boot_indices <- boot_indices[seq_len(n)]

    X_boot <- X[boot_indices, , drop = FALSE]

    pca_boot <- tryCatch({
      stats::prcomp(X_boot, center = center, scale. = scale)
    }, error = function(e) NULL)

    if (is.null(pca_boot)) next

    L_boot <- pca_boot$rotation[, seq_len(min(n_components, ncol(pca_boot$rotation))),
                                drop = FALSE]

    if (ncol(L_boot) < n_components) {
      L_boot <- cbind(L_boot, matrix(0, nrow = p, ncol = n_components - ncol(L_boot)))
    }

    if (rotation != "none" && n_components > 1) {
      if (rotation == "varimax") {
        rot_boot <- tryCatch(stats::varimax(L_boot), error = function(e) NULL)
        if (!is.null(rot_boot)) {
          L_boot <- rot_boot$loadings[, , drop = FALSE]
          class(L_boot) <- "matrix"
        }
      }
    }

    procrustes_result <- procrustes_rotation(L_boot, loadings)
    boot_loadings[, , b] <- procrustes_result

  }

  if (verbose) message("Computing bootstrap p-values...")

  p_values <- matrix(NA_real_, nrow = p, ncol = n_components)
  rownames(p_values) <- var_names
  colnames(p_values) <- paste0("PC", seq_len(n_components))

  abs_loadings <- abs(loadings)

  for (j in seq_len(n_components)) {
    for (i in seq_len(p)) {
      boot_vals <- abs(boot_loadings[i, j, ])
      boot_vals <- boot_vals[!is.na(boot_vals)]
      if (length(boot_vals) > 0) {
        p_values[i, j] <- mean(boot_vals >= abs_loadings[i, j])
      }
    }
  }

  if (use_fdr) {
    p_values_adj <- matrix(
      stats::p.adjust(as.vector(p_values), method = "BH"),
      nrow = p, ncol = n_components
    )
    rownames(p_values_adj) <- var_names
    colnames(p_values_adj) <- paste0("PC", seq_len(n_components))
    significant <- p_values_adj < alpha
  } else {
    thresholds <- apply(abs(boot_loadings), 2, function(x) {
      stats::quantile(x, 1 - alpha, na.rm = TRUE)
    })
    significant <- abs_loadings > matrix(thresholds, nrow = p, ncol = n_components,
                                         byrow = TRUE)
  }

  entropy_by_pc <- apply(loadings, 2, function(x) compute_entropy(x, normalize = TRUE))

  summary_df <- data.frame(
    component = paste0("PC", seq_len(n_components)),
    eigenvalue = eigenvalues[seq_len(n_components)],
    var_explained = var_explained[seq_len(n_components)],
    cumulative_var = cumvar[seq_len(n_components)],
    n_significant = colSums(significant),
    entropy = entropy_by_pc,
    stringsAsFactors = FALSE
  )

  assignments <- data.frame(
    variable = var_names,
    dominant_pc = apply(abs_loadings, 1, which.max),
    max_loading = apply(loadings, 1, function(x) x[which.max(abs(x))]),
    abs_max_loading = apply(abs_loadings, 1, max),
    is_significant = apply(significant, 1, any),
    stringsAsFactors = FALSE
  )
  assignments$dominant_pc_name <- paste0("PC", assignments$dominant_pc)

  rownames(loadings) <- var_names
  colnames(loadings) <- paste0("PC", seq_len(n_components))

  result <- list(
    loadings = loadings,
    scores = scores,
    eigenvalues = eigenvalues,
    variance_explained = var_explained,
    cumulative_variance = cumvar,
    significant_loadings = significant,
    p_values = p_values,
    n_components = n_components,
    entropy = entropy_by_pc,
    summary_by_component = summary_df,
    assignments = assignments,
    settings = list(
      center = center,
      scale = scale,
      rotation = rotation,
      n_boot = n_boot,
      block_length = block_length,
      alpha = alpha,
      use_fdr = use_fdr
    )
  )

  class(result) <- c("signaly_pca", "list")

  if (verbose) {
    cat("\n--- PCA SUMMARY ---\n")
    print(summary_df)
    cat(sprintf("\nTotal significant variable-component pairs: %d/%d\n",
                sum(significant), p * n_components))
  }

  result
}


#' Procrustes Rotation for Bootstrap Alignment
#'
#' @description
#' Rotates matrix B to best match target matrix A using orthogonal Procrustes
#' rotation.
#'
#' @param B Matrix to rotate.
#' @param A Target matrix.
#'
#' @return Rotated matrix B.
#'
#' @keywords internal
procrustes_rotation <- function(B, A) {
  svd_result <- svd(t(A) %*% B)
  rotation_matrix <- svd_result$v %*% t(svd_result$u)
  B %*% rotation_matrix
}


#' Dynamic Factor Model Estimation
#'
#' @description
#' Estimates a Dynamic Factor Model (DFM) to extract common latent factors
#' from panel data. Uses principal components as initial estimates and
#' optionally refines via EM algorithm.
#'
#' @param X Matrix or data frame where rows are observations and columns are
#'   variables.
#' @param r Number of factors. If NULL, determined by information criterion.
#' @param p Number of lags in factor VAR dynamics. Default 1.
#' @param ic Character string specifying information criterion for factor
#'   selection: "IC1", "IC2", or "IC3" (Bai & Ng, 2002). Default "IC2".
#' @param max_factors Maximum number of factors to consider. Default
#'   min(10, floor(ncol(X)/2)).
#' @param standardize Logical. Standardize variables before estimation.
#'   Default TRUE.
#' @param verbose Logical for progress messages.
#'
#' @return A list of class "signaly_dfm" containing:
#' \describe{
#'   \item{factors}{Matrix of estimated latent factors (T x r)}
#'   \item{loadings}{Matrix of factor loadings (p x r)}
#'   \item{var_coefficients}{VAR coefficient matrices for factor dynamics}
#'   \item{idiosyncratic_var}{Idiosyncratic variance estimates}
#'   \item{r_selected}{Number of factors selected}
#'   \item{ic_values}{Information criterion values}
#'   \item{fitted_values}{Fitted values from the model}
#'   \item{residuals}{Residuals (idiosyncratic components)}
#' }
#'
#' @details
#' The DFM assumes:
#' \deqn{X_{it} = \lambda_i' F_t + e_{it}}
#'
#' where \eqn{F_t} are common factors, \eqn{\lambda_i} are loadings, and
#' \eqn{e_{it}} are idiosyncratic errors. The factors follow VAR dynamics:
#' \deqn{F_t = A_1 F_{t-1} + ... + A_p F_{t-p} + u_t}
#'
#' Factor selection uses the Bai & Ng (2002) information criteria which
#' penalize over-fitting while consistently estimating the true number of
#' factors.
#'
#' @references
#' Bai, J., & Ng, S. (2002). Determining the number of factors in approximate
#' factor models. Econometrica, 70(1), 191-221.
#'
#' Stock, J. H., & Watson, M. W. (2002). Forecasting using principal components
#' from a large number of predictors. Journal of the American Statistical
#' Association, 97(460), 1167-1179.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), ncol = p)
#' result <- estimate_dfm(X, r = 2)
#' print(dim(result$factors))
#'
#' @export
estimate_dfm <- function(X,
                         r = NULL,
                         p = 1,
                         ic = c("IC2", "IC1", "IC3"),
                         max_factors = NULL,
                         standardize = TRUE,
                         verbose = FALSE) {

  ic <- match.arg(ic)

  X <- as.matrix(X)
  n <- nrow(X)
  k <- ncol(X)

  var_names <- colnames(X)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(k))

  if (standardize) {
    X_std <- scale(X)
    X_mean <- attr(X_std, "scaled:center")
    X_sd <- attr(X_std, "scaled:scale")
    X_work <- X_std
  } else {
    X_work <- X
    X_mean <- rep(0, k)
    X_sd <- rep(1, k)
  }

  if (is.null(max_factors)) {
    max_factors <- min(10, floor(k / 2), floor(n / 2) - 1)
  }

  if (is.null(r)) {
    if (verbose) message("Selecting number of factors via information criteria...")

    ic_values <- numeric(max_factors)

    for (r_test in seq_len(max_factors)) {
      svd_result <- svd(X_work, nu = r_test, nv = r_test)
      F_hat <- sqrt(n) * svd_result$u[, seq_len(r_test), drop = FALSE]
      Lambda_hat <- t(X_work) %*% F_hat / n

      X_hat <- F_hat %*% t(Lambda_hat)
      V_r <- sum((X_work - X_hat)^2) / (n * k)

      if (ic == "IC1") {
        penalty <- r_test * ((n + k) / (n * k)) * log((n * k) / (n + k))
      } else if (ic == "IC2") {
        penalty <- r_test * ((n + k) / (n * k)) * log(min(n, k))
      } else {
        penalty <- r_test * (log(min(n, k)) / min(n, k))
      }

      ic_values[r_test] <- log(V_r) + penalty
    }

    r <- which.min(ic_values)
    if (verbose) {
      message(sprintf("Selected r = %d factors (minimum %s)", r, ic))
    }
  }

  r <- min(r, max_factors, k - 1, n - 1)

  if (verbose) message(sprintf("Estimating DFM with %d factors...", r))

  svd_result <- svd(X_work, nu = r, nv = r)
  F_hat <- sqrt(n) * svd_result$u[, seq_len(r), drop = FALSE]
  Lambda_hat <- t(X_work) %*% F_hat / n

  colnames(F_hat) <- paste0("F", seq_len(r))
  rownames(Lambda_hat) <- var_names
  colnames(Lambda_hat) <- paste0("F", seq_len(r))

  X_fitted <- F_hat %*% t(Lambda_hat)
  residuals <- X_work - X_fitted
  idio_var <- apply(residuals, 2, stats::var)
  names(idio_var) <- var_names

  var_coef <- NULL
  if (p >= 1 && n > r * p + 1) {
    if (verbose) message(sprintf("Estimating VAR(%d) for factor dynamics...", p))

    F_lag_list <- vector("list", p)
    for (lag in seq_len(p)) {
      F_lag_list[[lag]] <- F_hat[seq_len(n - p) + (p - lag), , drop = FALSE]
    }
    F_lags <- do.call(cbind, F_lag_list)
    F_current <- F_hat[(p + 1):n, , drop = FALSE]

    var_fit <- stats::lm.fit(cbind(1, F_lags), F_current)
    var_coef <- var_fit$coefficients
  }

  ic_values_named <- NULL
  if (is.null(r)) {
    ic_values_named <- ic_values
    names(ic_values_named) <- seq_len(max_factors)
  }

  result <- list(
    factors = F_hat,
    loadings = Lambda_hat,
    var_coefficients = var_coef,
    idiosyncratic_var = idio_var,
    r_selected = r,
    ic_values = ic_values_named,
    fitted_values = X_fitted,
    residuals = residuals,
    eigenvalues = svd_result$d^2 / n,
    var_explained = svd_result$d^2 / sum(svd_result$d^2),
    settings = list(
      p = p,
      ic = ic,
      standardize = standardize,
      X_mean = X_mean,
      X_sd = X_sd
    )
  )

  class(result) <- c("signaly_dfm", "list")

  if (verbose) {
    cat(sprintf("\nDFM estimated with r=%d factors, p=%d lags\n", r, p))
    cat(sprintf("Variance explained by factors: %.1f%%\n",
                100 * sum(result$var_explained[seq_len(r)])))
  }

  result
}


#' @export
print.signaly_pca <- function(x, ...) {
  cat("\nSignalY PCA Analysis\n")
  cat("====================\n")
  cat(sprintf("Components extracted: %d\n", x$n_components))
  cat(sprintf("Total variance explained: %.1f%%\n",
              100 * x$cumulative_variance[x$n_components]))
  cat(sprintf("Rotation: %s\n", x$settings$rotation))
  cat("\nSummary by component:\n")
  print(x$summary_by_component, row.names = FALSE)
  invisible(x)
}


#' @export
print.signaly_dfm <- function(x, ...) {
  cat("\nSignalY Dynamic Factor Model\n")
  cat("============================\n")
  cat(sprintf("Factors: %d\n", x$r_selected))
  cat(sprintf("VAR lags: %d\n", x$settings$p))
  cat(sprintf("Variance explained: %.1f%%\n",
              100 * sum(x$var_explained[seq_len(x$r_selected)])))
  invisible(x)
}

#' @export
print.pca_bootstrap <- function(x, ...) {
  
  cat("PCA Bootstrap Results\n")
  cat(sprintf("  Variables: %d, Components: %d\n", nrow(x$loadings), ncol(x$loadings)))
  cat(sprintf("  Variance explained (first %d): %.1f%%\n", 
              min(3, length(x$variance_explained)), 
              100 * sum(x$variance_explained[1:min(3, length(x$variance_explained))])))
  invisible(x)
}

#' @export
print.dfm_result <- function(x, ...) {
  cat("Dynamic Factor Model\n
")
  cat(sprintf("  Factors: %d, Variables: %d, Observations: %d\n",
              ncol(x$factors), nrow(x$loadings), nrow(x$factors)))
  invisible(x)
}