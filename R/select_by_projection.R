#' Select Variables via Predictive Projection
#'
#' @description
#' Variable selection via predictive projection. Finds the minimum-size
#' submodel whose predictive performance is within 1 SE of the full reference
#' model (the fitted Regularized Horseshoe). This method is strictly superior
#' in robustness to marginal criteria (kappa, credible intervals) because it
#' evaluates \emph{joint} predictive performance of growing submodels,
#' captures redundancies among correlated predictors, and is immune to the
#' tau-collapse failure mode of kappa-based selection.
#'
#' Implements the Projection Predictive algorithm (Piironen, Paasiniemi &
#' Vehtari, 2020) directly for the Gaussian family.
#'
#' @param hs_fit Object of class \code{"signaly_horseshoe"} returned by
#'   \code{\link{fit_horseshoe}}.
#' @param cv_k Integer number of folds for cross-validated ELPD evaluation.
#'   Default: 5. Set to \code{NULL} or \code{0} for training ELPD only (faster
#'   but less robust).
#' @param n_draws_max Maximum number of posterior draws to use. If the
#'   posterior has more draws, a random subsample is taken for computational
#'   efficiency. Default: 400.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{selected}{Character vector of selected variable names at the
#'     suggested submodel size.}
#'   \item{suggested_size}{Integer indicating the suggested number of
#'     variables.}
#'   \item{ranking}{Character vector giving the full variable entry order from
#'     the forward search.}
#'   \item{elpd}{Data frame with columns \code{size}, \code{elpd},
#'     \code{elpd_diff} (vs full model), and \code{se_diff}.}
#'   \item{elpd_full}{Scalar ELPD of the full reference model.}
#'   \item{method}{\code{"projection"}}
#' }
#'
#' @details
#' The algorithm proceeds in two phases:
#'
#' \strong{Phase 1 — Forward Search (variable entry order):}
#' At each step \eqn{k}, the variable whose inclusion minimizes the squared
#' projection error \eqn{||\bar{\mu}_{ref} - X_S \hat{\beta}||^2} is added,
#' where \eqn{\bar{\mu}_{ref}} is the posterior mean of the reference
#' predictions and the projection is OLS with intercept.
#'
#' \strong{Phase 2 — ELPD evaluation (submodel quality):}
#' For each submodel size \eqn{k}, the full posterior is projected onto the
#' \eqn{k}-variable submodel. The expected log predictive density (ELPD) is
#' computed either on training data or via K-fold cross-validation. The
#' suggested size is the smallest \eqn{k} where the ELPD difference from the
#' full model is within 1 standard error of zero (the 1-SE rule).
#'
#' For the Gaussian family, the projection of draw \eqn{s} onto submodel
#' variables \eqn{S} is:
#' \deqn{\hat{\beta}_s = (X_S' X_S)^{-1} X_S' \mu_{ref,s}}
#' \deqn{\hat{\mu}_{proj,s} = X_S \hat{\beta}_s}
#'
#' @references
#' Piironen, J., Paasiniemi, M., & Vehtari, A. (2020). Projective inference
#' in high-dimensional problems: Prediction and feature selection. Electronic
#' Journal of Statistics, 14(1), 2155-2197.
#'
#' @examples
#' \dontrun{
#' hs_fit <- fit_horseshoe(y, X, p0 = 5)
#' proj <- select_by_projection(hs_fit, cv_k = 10, verbose = TRUE)
#' print(proj$selected)
#' print(proj$elpd)
#' }
#'
#' @seealso \code{\link{select_by_shrinkage}},
#'   \code{\link{select_by_credible_interval}},
#'   \code{\link{select_by_magnitude}},
#'   \code{\link{fit_horseshoe}}
#'
#' @export
select_by_projection <- function(hs_fit,
                                 cv_k = 5,
                                 n_draws_max = 400,
                                 verbose = TRUE) {
  
  # =========================================================================
  # 1. EXTRACT DATA AND POSTERIOR DRAWS
  # =========================================================================
  validate_horseshoe_fit(hs_fit, "select_by_projection")

  if (is.null(hs_fit$stan_data)) {
    stop(paste0(
      "hs_fit$stan_data not found. This field is required for projection.\n",
      "Please re-fit your model with the current version of fit_horseshoe()."
    ), call. = FALSE)
  }
  
  y <- as.vector(hs_fit$stan_data$y)
  X <- as.matrix(hs_fit$stan_data$X)
  var_names <- hs_fit$var_names
  
  N <- length(y)
  P <- ncol(X)
  
  beta_draws  <- hs_fit$fit$draws("beta",  format = "draws_matrix")  # S_full x P
  alpha_draws <- as.vector(hs_fit$fit$draws("alpha", format = "draws_matrix"))
  sigma_draws <- as.vector(hs_fit$fit$draws("sigma", format = "draws_matrix"))
  
  S_full <- nrow(beta_draws)
  
  # Subsample draws for efficiency
  if (S_full > n_draws_max) {
    idx_draws <- with_seed(42L, {
      sort(sample.int(S_full, n_draws_max))
    })
    beta_draws  <- beta_draws[idx_draws, , drop = FALSE]
    alpha_draws <- alpha_draws[idx_draws]
    sigma_draws <- sigma_draws[idx_draws]
  }
  S <- nrow(beta_draws)
  
  # Reference predictions: S x N matrix
  mu_ref <- sweep(beta_draws %*% t(X), 1, alpha_draws, "+")
  
  if (verbose) {
    message("============================================================")
    message("PROJECTION PREDICTIVE VARIABLE SELECTION")
    message("============================================================")
    message(sprintf("N = %d observations, P = %d predictors, S = %d draws", N, P, S))
  }
  
  # =========================================================================
  # 2. FORWARD SEARCH (determine variable entry order)
  #    At each step, add the variable that minimizes the squared projection
  #    error of the reference mean onto the submodel (OLS with intercept).
  # =========================================================================
  mu_ref_mean <- colMeans(mu_ref)  # N-vector
  
  selected_idx  <- integer(0)
  remaining_idx <- seq_len(P)
  entry_order   <- integer(P)
  
  if (verbose) message("\nPhase 1: Forward search...")
  
  for (step in seq_len(P)) {
    best_j    <- NULL
    best_loss <- Inf
    
    for (j in remaining_idx) {
      S_j   <- c(selected_idx, j)
      X_sub <- cbind(1, X[, S_j, drop = FALSE])  # N x (k+1), with intercept
      
      # OLS projection: minimize ||mu_ref_mean - X_sub %*% beta||^2
      beta_proj <- tryCatch(
        solve(crossprod(X_sub), crossprod(X_sub, mu_ref_mean)),
        error = function(e) {
          XtX <- crossprod(X_sub)
          diag(XtX) <- diag(XtX) + 1e-8
          solve(XtX, crossprod(X_sub, mu_ref_mean))
        }
      )
      residual <- mu_ref_mean - X_sub %*% beta_proj
      loss <- sum(residual^2)
      
      if (loss < best_loss) {
        best_loss <- loss
        best_j    <- j
      }
    }
    
    selected_idx  <- c(selected_idx, best_j)
    remaining_idx <- setdiff(remaining_idx, best_j)
    entry_order[step] <- best_j
    
    if (verbose && (step <= 10 || step %% 10 == 0 || step == P)) {
      message(sprintf("  Step %2d: + %-15s  (residual SS = %.6f)",
                      step, var_names[best_j], best_loss))
    }
  }
  
  if (verbose) message("\nPhase 2: ELPD evaluation...")
  
  # =========================================================================
  # 3. UTILITY: log(mean(exp(x))) with numerical stability
  # =========================================================================
  log_mean_exp <- function(x) {
    mx <- max(x)
    mx + log(mean(exp(x - mx)))
  }
  
  # =========================================================================
  # 4. HELPER: project draws and compute pointwise ELPD
  #    Given design matrices for projection and evaluation, returns
  #    pointwise ELPD for each observation in the evaluation set.
  # =========================================================================
  project_elpd <- function(X_proj, X_eval, y_eval, mu_ref_proj) {
    # X_proj: n_proj x (k+1) — design matrix used for projection (with intercept)
    # X_eval: n_eval x (k+1) — design matrix for evaluation (with intercept)
    # y_eval: n_eval — response for evaluation
    # mu_ref_proj: S x n_proj — reference predictions on projection data
    
    XtX_inv <- tryCatch(
      solve(crossprod(X_proj)),
      error = function(e) {
        XtX <- crossprod(X_proj)
        diag(XtX) <- diag(XtX) + 1e-8
        solve(XtX)
      }
    )
    
    # Project each draw: beta_proj[s,] = XtX_inv %*% X_proj' %*% mu_ref_proj[s,]
    # => beta_proj = mu_ref_proj %*% X_proj %*% XtX_inv   (S x (k+1))
    beta_proj <- mu_ref_proj %*% X_proj %*% XtX_inv
    
    # Predicted on eval data: S x n_eval
    mu_proj_eval <- beta_proj %*% t(X_eval)
    
    # Pointwise ELPD
    n_eval <- length(y_eval)
    elpd_pw <- numeric(n_eval)
    for (i in seq_len(n_eval)) {
      ll <- stats::dnorm(y_eval[i], mu_proj_eval[, i], sigma_draws, log = TRUE)
      elpd_pw[i] <- log_mean_exp(ll)
    }
    elpd_pw
  }
  
  # =========================================================================
  # 5. COMPUTE ELPD AT EACH SUBMODEL SIZE
  # =========================================================================
  use_cv <- !is.null(cv_k) && cv_k > 1
  
  if (use_cv) {
    # -- K-fold cross-validated ELPD --
    if (verbose) message(sprintf("  Using %d-fold cross-validation...", cv_k))
    
    # Create folds (deterministic, balanced)
    fold_ids <- with_seed(123L, {
      ids <- rep(seq_len(cv_k), length.out = N)
      ids[sample.int(N)]
    })
    
    elpd_pointwise_full    <- numeric(N)
    elpd_pointwise_by_size <- matrix(0, nrow = N, ncol = P)
    
    for (f in seq_len(cv_k)) {
      if (verbose) message(sprintf("    Fold %d/%d...", f, cv_k))
      
      test_idx  <- which(fold_ids == f)
      train_idx <- which(fold_ids != f)
      
      # Full model ELPD on test data
      for (i in test_idx) {
        ll <- stats::dnorm(y[i], mu_ref[, i], sigma_draws, log = TRUE)
        elpd_pointwise_full[i] <- log_mean_exp(ll)
      }
      
      # Submodel ELPD at each size
      for (k in seq_len(P)) {
        vars_k      <- entry_order[seq_len(k)]
        X_train_sub <- cbind(1, X[train_idx, vars_k, drop = FALSE])
        X_test_sub  <- cbind(1, X[test_idx,  vars_k, drop = FALSE])
        
        elpd_pw <- project_elpd(
          X_proj      = X_train_sub,
          X_eval      = X_test_sub,
          y_eval      = y[test_idx],
          mu_ref_proj = mu_ref[, train_idx, drop = FALSE]
        )
        
        elpd_pointwise_by_size[test_idx, k] <- elpd_pw
      }
    }
    
  } else {
    # -- Training ELPD (no cross-validation) --
    if (verbose) message("  Using training ELPD (no CV)...")
    
    elpd_pointwise_full    <- numeric(N)
    elpd_pointwise_by_size <- matrix(0, nrow = N, ncol = P)
    
    # Full model
    for (i in seq_len(N)) {
      ll <- stats::dnorm(y[i], mu_ref[, i], sigma_draws, log = TRUE)
      elpd_pointwise_full[i] <- log_mean_exp(ll)
    }
    
    # Each submodel size
    for (k in seq_len(P)) {
      vars_k <- entry_order[seq_len(k)]
      X_sub  <- cbind(1, X[, vars_k, drop = FALSE])
      
      elpd_pw <- project_elpd(
        X_proj      = X_sub,
        X_eval      = X_sub,
        y_eval      = y,
        mu_ref_proj = mu_ref
      )
      
      elpd_pointwise_by_size[, k] <- elpd_pw
    }
  }
  
  elpd_full    <- sum(elpd_pointwise_full)
  elpd_by_size <- colSums(elpd_pointwise_by_size)
  
  # =========================================================================
  # 6. SUGGEST SIZE (1-SE rule)
  #    Smallest k where ELPD difference from full model is within 1 SE.
  # =========================================================================
  elpd_diff <- elpd_full - elpd_by_size
  se_diff   <- numeric(P)
  for (k in seq_len(P)) {
    diff_pw   <- elpd_pointwise_full - elpd_pointwise_by_size[, k]
    se_diff[k] <- sqrt(N * stats::var(diff_pw))
  }
  
  suggested <- P
  for (k in seq_len(P)) {
    if (elpd_diff[k] <= se_diff[k]) {
      suggested <- k
      break
    }
  }
  
  # =========================================================================
  # 7. BUILD RESULTS
  # =========================================================================
  elpd_df <- data.frame(
    size      = seq_len(P),
    elpd      = elpd_by_size,
    elpd_diff = elpd_diff,
    se_diff   = se_diff
  )
  
  selected_vars <- var_names[entry_order[seq_len(suggested)]]
  ranking_vars  <- var_names[entry_order]
  
  if (verbose) {
    message(sprintf("\n--- PROJECTION RESULT ---"))
    message(sprintf("Full model ELPD: %.2f", elpd_full))
    message(sprintf("Suggested size: %d variables (1-SE rule)", suggested))
    message("Selected variables:")
    for (i in seq_along(selected_vars)) {
      message(sprintf("  %2d. %-15s  (ELPD = %.2f, diff = %.2f, SE = %.2f)",
                      i, selected_vars[i],
                      elpd_by_size[i], elpd_diff[i], se_diff[i]))
    }
    message("============================================================")
  }
  
  list(
    selected       = selected_vars,
    suggested_size = suggested,
    ranking        = ranking_vars,
    elpd           = elpd_df,
    elpd_full      = elpd_full,
    method         = "projection"
  )
}