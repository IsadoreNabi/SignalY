#' @title Regularized Horseshoe Regression for Variable Selection
#' @description
#' Implements Bayesian sparse regression using the regularized Horseshoe prior
#' for identifying structurally relevant predictors from high-dimensional
#' candidate variable sets.
#'
#' @name horseshoe
NULL


#' Fit Regularized Horseshoe Regression Model
#'
#' @description
#' Fits a Bayesian linear regression with regularized Horseshoe prior using
#' Stan via cmdstanr. This version includes improved numerical stability
#' and automatic prior calibration.
#'
#' @details
#' The regularized Horseshoe prior (Piironen & Vehtari, 2017) provides adaptive
#' shrinkage that can distinguish between relevant and irrelevant predictors.
#' 
#' \strong{Variable Selection Methods:}
#' 
#' After fitting, variables can be selected using different criteria:
#' \itemize{
#'   \item \code{\link{select_by_credible_interval}}: Selects variables whose
#'     credible interval excludes zero. \strong{Recommended} - most robust method.
#'   \item \code{\link{select_by_shrinkage}}: Selects based on kappa (shrinkage factor).
#'     May underselect when tau is very small.
#'   \item \code{\link{select_by_magnitude}}: Selects based on coefficient magnitude.
#' }
#' 
#' \strong{Note on kappa-based selection:}
#' 
#' The shrinkage factor kappa depends on the global shrinkage parameter tau.
#' In some datasets, the posterior of tau may concentrate near zero, causing
#' all kappa values to be close to 1 even for truly relevant variables.
#' When this happens, the coefficient estimates (beta) remain valid, but
#' kappa-based selection will fail. The function automatically warns when
#' this occurs and recommends using \code{select_by_credible_interval()} instead.
#'
#' @param y Numeric vector of the response variable.
#' @param X Matrix or data frame of predictor variables.
#' @param var_names Optional character vector of variable names.
#' @param p0 Expected number of non-zero coefficients. Default: P/3.
#' @param slab_scale Scale for the regularizing slab. Default: 3.
#' @param slab_df Degrees of freedom for the slab. Default: 4.
#' @param tau_scale Scale multiplier for the global shrinkage prior. 
#'   Default: NULL (auto-calibrated based on data characteristics).
#'   Increase this value (e.g., 10-20) if the model over-shrinks.
#' @param use_qr Use QR decomposition? Default: FALSE.
#' @param standardize Standardize predictors internally? Default: TRUE.
#' @param X_new Optional matrix for out-of-sample prediction.
#' @param iter_warmup Warmup iterations per chain. Default: 1000.
#' @param iter_sampling Sampling iterations per chain. Default: 1000.
#' @param chains Number of MCMC chains. Default: 4.
#' @param adapt_delta Target acceptance probability. Default: 0.95.
#' @param max_treedepth Maximum tree depth. Default: 12.
#' @param seed Random seed.
#' @param verbose Print progress messages?
#'
#' @return A list of class "signaly_horseshoe" with posterior summaries,
#'   diagnostics, and model fit object.
#'
#' @export
fit_horseshoe <- function(y,
                          X,
                          var_names = NULL,
                          p0 = NULL,
                          slab_scale = 3,
                          slab_df = 4,
                          tau_scale = NULL,
                          use_qr = FALSE,
                          standardize = TRUE,
                          X_new = NULL,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          chains = 4,
                          adapt_delta = 0.95,
                          max_treedepth = 12,
                          seed = 123,
                          verbose = TRUE) {
  
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required. Install from https://mc-stan.org/r-packages/",
         call. = FALSE)
  }
  
  cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  if (is.null(cmdstan_path) || nchar(cmdstan_path) == 0) {
    stop("CmdStan not found. Install with cmdstanr::install_cmdstan()", call. = FALSE)
  }
  
  if (use_qr) {
    warning("QR decomposition destroys sparsity structure. Consider use_qr = FALSE for variable selection.")
  }
  
  y <- as.numeric(y)
  X <- as.matrix(X)
  
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must equal length of y.", call. = FALSE)
  }
  
  N <- length(y)
  P <- ncol(X)
  
  if (is.null(var_names)) {
    var_names <- if (!is.null(colnames(X))) colnames(X) else paste0("V", seq_len(P))
  }
  colnames(X) <- var_names
  
  if (is.null(p0)) {
    p0 <- max(1, P / 3)
  }
  
  X_sds <- apply(X, 2, sd)
  X_means <- apply(X, 2, mean)
  data_looks_standardized <- all(abs(X_sds - 1) < 0.3) && all(abs(X_means) < 0.3)
  
  if (is.null(tau_scale)) {
    p0_ratio <- p0 / P
    
    tau_scale <- 1.0 / max(p0_ratio, 0.05)
    tau_scale <- max(5, min(tau_scale, 100))
    
    if (verbose) {
      message("\n+-------------------------------------------------------------+")
      message("| AUTO-CALIBRATING tau_scale                                  |")
      message("+-------------------------------------------------------------+")
      message(sprintf("| Data appears pre-standardized: %-3s                         |",
                      ifelse(data_looks_standardized, "YES", "NO")))
      message(sprintf("| Expected sparsity (p0/P): %.2f                              |", p0_ratio))
      message(sprintf("| Auto-selected tau_scale: %.1f                               |", tau_scale))
      message("|                                                             |")
      message("| If over-shrinking persists, manually increase tau_scale    |")
      message("+-------------------------------------------------------------+\n")
    }
  }
  
  tau_scale <- max(1, min(tau_scale, 200))
  
  if (!is.null(X_new)) {
    X_new_mat <- as.matrix(X_new)
    if (ncol(X_new_mat) != P) {
      stop("X_new must have the same number of columns as X.", call. = FALSE)
    }
    N_new <- nrow(X_new_mat)
  } else {
    N_new <- 0
    X_new_mat <- matrix(0, nrow = 1, ncol = P)
  }
  
  stan_file <- system.file("stan", "regularized_horseshoe.stan", package = "SignalY")
  
  if (stan_file == "" || !file.exists(stan_file)) {
    stan_code <- get_horseshoe_stan_code()
    stan_file <- file.path(tempdir(), "regularized_horseshoe.stan")
    writeLines(stan_code, stan_file)
  }
  
  if (verbose) {
    message("============================================================")
    message("REGULARIZED HORSESHOE REGRESSION")
    message("============================================================")
    message(sprintf("N = %d, P = %d, p0 = %.1f", N, P, p0))
    message(sprintf("tau_scale = %.1f, slab_scale = %.1f, slab_df = %.1f",
                    tau_scale, slab_scale, slab_df))
    message(sprintf("Chains: %d, Warmup: %d, Sampling: %d",
                    chains, iter_warmup, iter_sampling))
  }
  
  stan_data <- list(
    N = N,
    P = P,
    X = X,
    y = y,
    p0 = as.numeric(p0),
    slab_scale = as.numeric(slab_scale),
    slab_df = as.numeric(slab_df),
    tau_scale = as.numeric(tau_scale),
    use_qr = as.integer(use_qr),
    standardize = as.integer(standardize),
    N_new = N_new,
    X_new = X_new_mat
  )
  
  if (verbose) message("\nCompiling Stan model...")
  
  model <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = TRUE),
    stanc_options = list("O1")
  )
  
  if (verbose) message("Sampling from posterior...")
  
  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = min(chains, safe_cores(reserve = 1L)),
    threads_per_chain = 1,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = ifelse(verbose, 500, 0),
    show_messages = verbose
  )
  
  if (verbose) message("\nProcessing results...")
  
  diagnostics <- extract_mcmc_diagnostics(fit, verbose = verbose)
  
  results <- process_horseshoe_results(fit, var_names, verbose = verbose)
  
  if (verbose) message("Computing LOO-CV...")
  loo_result <- compute_horseshoe_loo(fit, verbose = verbose)
  
  if (verbose) message("Generating posterior predictive checks...")
  ppc <- posterior_predictive_check_horseshoe(fit, y, verbose = verbose)
  
  sparsity <- list(
    n_relevant = sum(results$coefficients$kappa_mean < 0.5),
    n_irrelevant = sum(results$coefficients$kappa_mean >= 0.5),
    sparsity_ratio = mean(results$coefficients$kappa_mean >= 0.5),
    m_eff_mean = mean(results$m_eff_draws),
    m_eff_ci = stats::quantile(results$m_eff_draws, c(0.025, 0.975))
  )
  
  output <- list(
    coefficients = results$coefficients,
    hyperparameters = results$hyperparams,
    diagnostics = diagnostics,
    loo = loo_result,
    ppc = ppc,
    posterior_draws = list(
      beta = results$beta_draws,
      kappa = results$kappa_draws,
      tau = results$tau_draws,
      sigma = results$sigma_draws,
      m_eff = results$m_eff_draws
    ),
    fit = fit,
    stan_data = stan_data,
    sparsity = sparsity,
    var_names = var_names,
    settings = list(
      p0 = p0,
      slab_scale = slab_scale,
      slab_df = slab_df,
      tau_scale = tau_scale,
      use_qr = use_qr,
      standardize = standardize,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      chains = chains,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      seed = seed
    )
  )
  
  if (N_new > 0) {
    output$predictions <- fit$draws("y_new_rep", format = "draws_matrix")
  }
  
  class(output) <- c("signaly_horseshoe", "list")
  
  if (verbose) {
    message("\n============================================================")
    message("SUMMARY")
    message("============================================================")
    message(sprintf("Effective non-zeros (m_eff): %.2f [%.2f, %.2f]",
                    sparsity$m_eff_mean, sparsity$m_eff_ci[1], sparsity$m_eff_ci[2]))
    message(sprintf("Variables with kappa < 0.5: %d / %d", 
                    sparsity$n_relevant, P))
    
    message(sprintf("\nPosterior tau: %.3f [%.3f, %.3f]",
                    mean(results$tau_draws),
                    stats::quantile(results$tau_draws, 0.025),
                    stats::quantile(results$tau_draws, 0.975)))
    
    coefs <- results$coefficients
    coefs_sorted <- coefs[order(coefs$kappa_mean), ]
    n_show <- min(10, nrow(coefs_sorted))
    
    message("\nTop variables by relevance (lowest kappa):")
    for (i in seq_len(n_show)) {
      message(sprintf("  %s: beta=%.3f, kappa=%.3f",
                      coefs_sorted$variable[i],
                      coefs_sorted$beta_mean[i],
                      coefs_sorted$kappa_mean[i]))
    }
    
    if (diagnostics$n_divergences > 0) {
      message(sprintf("\n[!] Warning: %d divergences detected", diagnostics$n_divergences))
    }
    
    tau_mean <- mean(results$tau_draws)
    if (tau_mean < 0.05) {
      message("\n[!] Note: Posterior tau is very small (", round(tau_mean, 4), ").")
      message("    The kappa-based selection (select_by_shrinkage) may underselect variables.")
      message("    Recommendation: Use select_by_credible_interval() for variable selection.")
      message("    The coefficient estimates (beta) remain valid and well-calibrated.")
    }
  }
  
  output
}

#' Extract MCMC Diagnostics
#'
#' @description
#' Extracts comprehensive MCMC diagnostics from a cmdstanr fit object.
#'
#' @param fit cmdstanr fit object.
#' @param verbose Logical for messages.
#'
#' @return List of diagnostic metrics.
#'
#' @keywords internal
extract_mcmc_diagnostics <- function(fit, verbose = FALSE) {
  
  sampler_diag <- fit$sampler_diagnostics()
  
  n_divergences <- sum(sampler_diag[, , "divergent__"])
  
  treedepth <- sampler_diag[, , "treedepth__"]
  max_td <- fit$metadata()$max_treedepth
  max_treedepth_hit <- sum(treedepth >= max_td)
  
  energy <- sampler_diag[, , "energy__"]
  energy_bfmi <- apply(energy, 2, function(e) {
    diff_e <- diff(e)
    stats::var(diff_e) / stats::var(e)
  })
  low_bfmi <- sum(energy_bfmi < 0.3)
  
  summ <- fit$summary()
  rhat_vals <- summ$rhat
  rhat_max <- max(rhat_vals, na.rm = TRUE)
  rhat_bad <- sum(rhat_vals > 1.01, na.rm = TRUE)
  
  ess_bulk <- summ$ess_bulk
  ess_tail <- summ$ess_tail
  ess_bulk_min <- min(ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(ess_tail, na.rm = TRUE)
  
  quality <- evaluate_fit_quality_internal(
    n_divergences, max_treedepth_hit, low_bfmi, rhat_max, ess_bulk_min, ess_tail_min
  )
  
  diagnostics <- list(
    n_divergences = n_divergences,
    max_treedepth_hit = max_treedepth_hit,
    energy_bfmi = energy_bfmi,
    low_bfmi = low_bfmi,
    rhat_max = rhat_max,
    rhat_bad = rhat_bad,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    quality = quality
  )
  
  if (verbose) {
    message(sprintf("  Divergences: %d", n_divergences))
    message(sprintf("  Max treedepth reached: %d times", max_treedepth_hit))
    message(sprintf("  Chains with low BFMI: %d", low_bfmi))
    message(sprintf("  Max R-hat: %.4f (parameters with R-hat > 1.01: %d)",
                    rhat_max, rhat_bad))
    message(sprintf("  Min ESS bulk: %.0f, Min ESS tail: %.0f",
                    ess_bulk_min, ess_tail_min))
    message(sprintf("  Quality: %s", quality$status))
  }
  
  diagnostics
}


#' Evaluate Fit Quality
#'
#' @keywords internal
evaluate_fit_quality_internal <- function(n_divergences, max_treedepth_hit,
                                          low_bfmi, rhat_max, ess_bulk_min,
                                          ess_tail_min) {
  issues <- character(0)
  
  if (n_divergences > 0) {
    issues <- c(issues, sprintf("DIVERGENCES: %d (indicates geometric problems)",
                                n_divergences))
  }
  
  if (max_treedepth_hit > 0) {
    issues <- c(issues, sprintf("MAX_TREEDEPTH: reached %d times",
                                max_treedepth_hit))
  }
  
  if (low_bfmi > 0) {
    issues <- c(issues, sprintf("Low BFMI in %d chains", low_bfmi))
  }
  
  if (rhat_max > 1.05) {
    issues <- c(issues, sprintf("High R-hat: %.3f (chains not mixing)", rhat_max))
  } else if (rhat_max > 1.01) {
    issues <- c(issues, sprintf("Marginal R-hat: %.3f", rhat_max))
  }
  
  if (ess_bulk_min < 400) {
    issues <- c(issues, sprintf("Low ESS bulk: %.0f", ess_bulk_min))
  }
  
  if (ess_tail_min < 400) {
    issues <- c(issues, sprintf("Low ESS tail: %.0f", ess_tail_min))
  }
  
  if (length(issues) == 0) {
    status <- "EXCELLENT - All metrics within acceptable ranges"
    passed <- TRUE
  } else if (n_divergences == 0 && rhat_max < 1.05) {
    status <- "ACCEPTABLE - Valid fit with minor warnings"
    passed <- TRUE
  } else {
    status <- "PROBLEMATIC - Review fit before interpretation"
    passed <- FALSE
  }
  
  list(status = status, passed = passed, issues = issues)
}


#' Process Horseshoe Results
#'
#' @keywords internal
process_horseshoe_results <- function(fit, var_names, verbose = FALSE) {
  
  P <- length(var_names)
  
  draws <- fit$draws(format = "draws_df")
  
  beta_vars <- paste0("beta[", seq_len(P), "]")
  beta_draws <- as.matrix(draws[, beta_vars])
  colnames(beta_draws) <- var_names
  
  kappa_vars <- paste0("kappa[", seq_len(P), "]")
  kappa_draws <- as.matrix(draws[, kappa_vars])
  colnames(kappa_draws) <- var_names
  
  m_eff_draws <- draws$m_eff
  tau_draws <- draws$tau
  sigma_draws <- draws$sigma
  alpha_draws <- draws$alpha
  
  coefficients <- data.frame(
    variable = var_names,
    beta_mean = colMeans(beta_draws),
    beta_sd = apply(beta_draws, 2, stats::sd),
    beta_q025 = apply(beta_draws, 2, stats::quantile, 0.025),
    beta_q10 = apply(beta_draws, 2, stats::quantile, 0.10),
    beta_median = apply(beta_draws, 2, stats::median),
    beta_q90 = apply(beta_draws, 2, stats::quantile, 0.90),
    beta_q975 = apply(beta_draws, 2, stats::quantile, 0.975),
    kappa_mean = colMeans(kappa_draws),
    kappa_median = apply(kappa_draws, 2, stats::median),
    stringsAsFactors = FALSE
  )
  
  threshold <- 0.01 * mean(abs(coefficients$beta_mean[coefficients$beta_mean != 0]))
  if (!is.finite(threshold) || threshold == 0) threshold <- 1e-6
  
  coefficients$prob_positive <- colMeans(beta_draws > threshold)
  coefficients$prob_negative <- colMeans(beta_draws < -threshold)
  coefficients$prob_nonzero <- coefficients$prob_positive + coefficients$prob_negative
  
  coefficients$relevance_score <- ifelse(
    coefficients$prob_positive > 0.5,
    (1 - coefficients$kappa_mean) * coefficients$prob_positive,
    ifelse(
      coefficients$prob_negative > 0.5,
      (1 - coefficients$kappa_mean) * coefficients$prob_negative,
      0
    )
  )
  
  coefficients <- coefficients[order(-coefficients$relevance_score), ]
  
  hyperparams <- data.frame(
    parameter = c("tau", "sigma", "alpha", "m_eff"),
    mean = c(mean(tau_draws), mean(sigma_draws), mean(alpha_draws), mean(m_eff_draws)),
    sd = c(stats::sd(tau_draws), stats::sd(sigma_draws), stats::sd(alpha_draws), stats::sd(m_eff_draws)),
    q025 = c(stats::quantile(tau_draws, 0.025), stats::quantile(sigma_draws, 0.025),
             stats::quantile(alpha_draws, 0.025), stats::quantile(m_eff_draws, 0.025)),
    median = c(stats::median(tau_draws), stats::median(sigma_draws),
               stats::median(alpha_draws), stats::median(m_eff_draws)),
    q975 = c(stats::quantile(tau_draws, 0.975), stats::quantile(sigma_draws, 0.975),
             stats::quantile(alpha_draws, 0.975), stats::quantile(m_eff_draws, 0.975)),
    stringsAsFactors = FALSE
  )
  
  list(
    coefficients = coefficients,
    hyperparams = hyperparams,
    beta_draws = beta_draws,
    kappa_draws = kappa_draws,
    tau_draws = tau_draws,
    sigma_draws = sigma_draws,
    m_eff_draws = m_eff_draws
  )
}


#' Compute LOO-CV for Horseshoe Model
#'
#' @keywords internal
compute_horseshoe_loo <- function(fit, verbose = FALSE) {
  
  if (!requireNamespace("loo", quietly = TRUE)) {
    if (verbose) message("Package 'loo' not available. Skipping LOO-CV.")
    return(NULL)
  }
  
  log_lik <- fit$draws("log_lik", format = "draws_matrix")
  
  loo_result <- loo::loo(log_lik, cores = min(2L, safe_cores(1L)))
  
  pareto_k <- loo_result$diagnostics$pareto_k
  k_bad <- sum(pareto_k > 0.7)
  
  if (verbose) {
    message(sprintf("  Pareto k > 0.7: %d/%d observations", k_bad, length(pareto_k)))
  }
  
  loo_result
}


#' Posterior Predictive Check for Horseshoe Model
#'
#' @keywords internal
posterior_predictive_check_horseshoe <- function(fit, y, verbose = FALSE) {
  
  y_rep <- fit$draws("y_rep", format = "draws_matrix")
  
  stats_result <- list(
    mean_obs = mean(y),
    mean_rep = mean(colMeans(y_rep)),
    sd_obs = stats::sd(y),
    sd_rep = mean(apply(y_rep, 1, stats::sd)),
    min_obs = min(y),
    min_rep = mean(apply(y_rep, 1, min)),
    max_obs = max(y),
    max_rep = mean(apply(y_rep, 1, max))
  )
  
  pvals <- list(
    mean = mean(rowMeans(y_rep) >= mean(y)),
    sd = mean(apply(y_rep, 1, stats::sd) >= stats::sd(y)),
    min = mean(apply(y_rep, 1, min) <= min(y)),
    max = mean(apply(y_rep, 1, max) >= max(y))
  )
  
  if (verbose) {
    message(sprintf("  Observed mean: %.2f, Predicted mean: %.2f",
                    stats_result$mean_obs, stats_result$mean_rep))
    message(sprintf("  Observed SD: %.2f, Predicted SD: %.2f",
                    stats_result$sd_obs, stats_result$sd_rep))
  }
  
  list(
    y_rep = y_rep,
    y_obs = y,
    stats = stats_result,
    pvals = pvals
  )
}


#' Print Horseshoe Summary
#'
#' @keywords internal
print_horseshoe_summary <- function(x) {
  
  cat("\n--- DIAGNOSTICS ---\n")
  diag <- x$diagnostics
  cat(sprintf("Divergences: %d\n", diag$n_divergences))
  cat(sprintf("Max R-hat: %.4f\n", diag$rhat_max))
  cat(sprintf("Min ESS (bulk/tail): %.0f / %.0f\n",
              diag$ess_bulk_min, diag$ess_tail_min))
  cat(sprintf("Quality: %s\n", diag$quality$status))
  
  cat("\n--- HYPERPARAMETERS ---\n")
  hp <- x$hyperparameters
  for (i in seq_len(nrow(hp))) {
    cat(sprintf("%s: %.3f (95%% CI: %.3f - %.3f)\n",
                hp$parameter[i], hp$mean[i], hp$q025[i], hp$q975[i]))
  }
  
  cat("\n--- SPARSITY ---\n")
  sp <- x$sparsity
  cat(sprintf("Relevant variables (kappa < 0.5): %d\n", sp$n_relevant))
  cat(sprintf("Irrelevant variables (kappa >= 0.5): %d\n", sp$n_irrelevant))
  cat(sprintf("Sparsity ratio: %.1f%%\n", 100 * sp$sparsity_ratio))
  cat(sprintf("Effective non-zeros (m_eff): %.1f (95%% CI: %.1f - %.1f)\n",
              sp$m_eff_mean, sp$m_eff_ci[1], sp$m_eff_ci[2]))
  
  cat("\n--- TOP RELEVANT VARIABLES (by relevance score) ---\n")
  coef <- x$coefficients
  top <- head(coef[order(-coef$relevance_score), ], 10)
  for (i in seq_len(nrow(top))) {
    cat(sprintf("%2d. %-20s beta=%.3f (kappa=%.2f, relevance=%.2f)\n",
                i, top$variable[i], top$beta_mean[i],
                top$kappa_mean[i], top$relevance_score[i]))
  }
  
  if (!is.null(x$loo)) {
    cat("\n--- LOO-CV ---\n")
    cat(sprintf("elpd_loo: %.1f (SE: %.1f)\n",
                x$loo$estimates["elpd_loo", "Estimate"],
                x$loo$estimates["elpd_loo", "SE"]))
  }
  
  invisible(NULL)
}


#' Get Stan Code for Regularized Horseshoe
#'
#' @description
#' Returns the Stan code for the regularized Horseshoe model with improved
#' numerical stability and prior calibration.
#'
#' @return Character string containing Stan code.
#'
#' @keywords internal
get_horseshoe_stan_code <- function() {
  '
// Regularized Horseshoe Prior - Robust Version
// Based on Piironen & Vehtari (2017) with numerical stability improvements

functions {
  // Compute regularized local scales with numerical safeguards
  vector regularized_hs_scale(vector lambda, real tau, real c2) {
    int P = rows(lambda);
    vector[P] lambda_tilde;
    real tau2 = square(tau) + 1e-10;  // Prevent division issues
    
    for (j in 1:P) {
      real lambda2 = square(lambda[j]) + 1e-10;
      real denom = c2 + tau2 * lambda2;
      lambda_tilde[j] = sqrt(c2 * lambda2 / denom);
    }
    return lambda_tilde;
  }
}

data {
  int<lower=1> N;                    // Number of observations
  int<lower=1> P;                    // Number of predictors
  matrix[N, P] X;                    // Predictor matrix
  vector[N] y;                       // Response vector
  real<lower=0> p0;                  // Prior number of non-zero coefficients
  real<lower=0> slab_scale;          // Scale for the regularizing slab
  real<lower=1> slab_df;             // Degrees of freedom for the slab
  real<lower=0> tau_scale;           // Scale multiplier for global shrinkage
  int<lower=0, upper=1> use_qr;      // Use QR decomposition?
  int<lower=0, upper=1> standardize; // Standardize predictors internally?
  int<lower=0> N_new;                // Number of new observations for prediction
  matrix[N_new > 0 ? N_new : 1, P] X_new;  // New predictor matrix
}

transformed data {
  matrix[N, P] X_work;
  matrix[N, P] Q_ast;
  matrix[P, P] R_ast;
  matrix[P, P] R_ast_inv;
  vector[P] X_mean = rep_vector(0.0, P);
  vector[P] X_sd = rep_vector(1.0, P);
  real y_mean = 0.0;
  real y_sd = 1.0;
  vector[N] y_work;
  real tau0;
  
  // Standardize if requested
  if (standardize == 1) {
    for (j in 1:P) {
      X_mean[j] = mean(X[, j]);
      X_sd[j] = sd(X[, j]);
      if (X_sd[j] < 1e-10) X_sd[j] = 1.0;
    }
    y_mean = mean(y);
    y_sd = sd(y);
    if (y_sd < 1e-10) y_sd = 1.0;
    
    for (j in 1:P) {
      X_work[, j] = (X[, j] - X_mean[j]) / X_sd[j];
    }
    y_work = (y - y_mean) / y_sd;
  } else {
    X_work = X;
    y_work = y;
  }
  
  // Calculate tau0 - the prior scale for global shrinkage
  // Key: tau0 should be O(1) for standardized data to avoid over-shrinkage
  {
    real p0_ratio = p0 / P;
    // Simple formula: tau0 = p0_ratio * tau_scale
    // With tau_scale auto-calibrated in R to give tau0 ~ 0.5-1.0
    tau0 = p0_ratio * tau_scale;
    
    // Bound tau0 to reasonable range
    if (tau0 < 0.1) tau0 = 0.1;
    if (tau0 > 10.0) tau0 = 10.0;
  }
  
  // QR decomposition if requested
  if (use_qr == 1) {
    Q_ast = qr_thin_Q(X_work) * sqrt(N - 1.0);
    R_ast = qr_thin_R(X_work) / sqrt(N - 1.0);
    R_ast_inv = inverse(R_ast);
  } else {
    Q_ast = X_work;
    R_ast = diag_matrix(rep_vector(1.0, P));
    R_ast_inv = diag_matrix(rep_vector(1.0, P));
  }
}

parameters {
  real alpha_std;                    // Intercept (standardized scale)
  vector[P] z;                       // Auxiliary for non-centered parameterization
  real<lower=0.01> tau;              // Global shrinkage - LOWER BOUNDED to prevent collapse
  vector<lower=0>[P] lambda;         // Local shrinkage
  real<lower=0> caux;                // Auxiliary for slab
  real<lower=0> sigma;               // Residual standard deviation
}

transformed parameters {
  real c2;                           // Squared slab scale
  vector[P] lambda_tilde;            // Regularized local scales
  vector[P] beta_tilde;              // Coefficients (working scale)
  vector[N] mu;                      // Linear predictor
  
  // Regularized slab
  c2 = square(slab_scale) * caux;
  
  // Compute regularized local scales
  lambda_tilde = regularized_hs_scale(lambda, tau, c2);
  
  // Non-centered parameterization for coefficients
  for (j in 1:P) {
    beta_tilde[j] = tau * lambda_tilde[j] * z[j];
  }
  
  // Linear predictor
  if (use_qr == 1) {
    mu = rep_vector(alpha_std, N) + Q_ast * beta_tilde;
  } else {
    mu = rep_vector(alpha_std, N) + X_work * beta_tilde;
  }
}

model {
  // Priors
  z ~ std_normal();
  
  // Global shrinkage - half-t with scale tau0
  // Using df=3 for slightly lighter tails than half-Cauchy
  tau ~ student_t(3, 0, tau0);
  
  // Local shrinkage - half-Cauchy(0, 1)
  lambda ~ student_t(1, 0, 1);
  
  // Slab degrees of freedom
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  // Intercept prior
  alpha_std ~ normal(0, 2);
  
  // Residual SD prior
  sigma ~ student_t(3, 0, 1);
  
  // Likelihood
  y_work ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N_new > 0 ? N_new : 1] y_new_rep;
  vector[P] kappa;                   // Shrinkage factors
  real m_eff;                        // Effective number of non-zeros
  vector[P] beta;                    // Coefficients (original scale)
  real alpha;                        // Intercept (original scale)
  real tau0_used = tau0;             // For debugging
  
  // Transform coefficients back to original scale
  if (use_qr == 1) {
    beta = R_ast_inv * beta_tilde;
  } else {
    beta = beta_tilde;
  }
  
  if (standardize == 1) {
    for (j in 1:P) {
      beta[j] = beta[j] * y_sd / X_sd[j];
    }
    alpha = y_mean + alpha_std * y_sd - dot_product(beta, X_mean);
  } else {
    alpha = alpha_std;
  }
  
  // Compute shrinkage factors kappa
  // kappa[j] = 1 / (1 + tau^2 * lambda_tilde[j]^2)
  // kappa near 0 = variable escapes shrinkage (relevant)
  // kappa near 1 = variable is shrunk (irrelevant)
  for (j in 1:P) {
    real tau_lambda_sq = square(tau * lambda_tilde[j]);
    kappa[j] = 1.0 / (1.0 + tau_lambda_sq);
  }
  
  // Effective number of non-zero coefficients
  m_eff = P - sum(kappa);
  
  // Log-likelihood for LOO-CV
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y_work[n] | mu[n], sigma);
  }
  
  // Posterior predictive samples
  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma);
    if (standardize == 1) {
      y_rep[n] = y_rep[n] * y_sd + y_mean;
    }
  }
  
  // Predictions for new data
  if (N_new > 0) {
    matrix[N_new, P] X_new_work;
    vector[N_new] mu_new;
    
    if (standardize == 1) {
      for (j in 1:P) {
        X_new_work[, j] = (X_new[, j] - X_mean[j]) / X_sd[j];
      }
    } else {
      X_new_work = X_new;
    }
    
    if (use_qr == 1) {
      mu_new = rep_vector(alpha_std, N_new) + (X_new_work * R_ast_inv) * beta_tilde;
    } else {
      mu_new = rep_vector(alpha_std, N_new) + X_new_work * beta_tilde;
    }
    
    for (n in 1:N_new) {
      y_new_rep[n] = normal_rng(mu_new[n], sigma);
      if (standardize == 1) {
        y_new_rep[n] = y_new_rep[n] * y_sd + y_mean;
      }
    }
  } else {
    y_new_rep[1] = 0.0;
  }
}
'
}


#' Select Variables Based on Shrinkage
#'
#' @description
#' Variable selection method using shrinkage factors (kappa).
#' Note: This method may underselect when tau collapses to small values.
#' Consider using \code{\link{select_by_credible_interval}} as an alternative.
#'
#' @param hs_fit Object returned by \code{\link{fit_horseshoe}}.
#' @param threshold Kappa threshold. Variables with kappa < threshold are
#'   considered relevant. Default 0.5.
#' @param verbose Logical for messages.
#'
#' @return List with selected variable names and details.
#'
#' @seealso \code{\link{select_by_credible_interval}} for a more robust alternative.
#'
#' @export
select_by_shrinkage <- function(hs_fit, threshold = 0.5, verbose = FALSE) {
  validate_horseshoe_fit(hs_fit, "select_by_shrinkage")

  coef_df <- hs_fit$coefficients
  
  selected <- coef_df[coef_df$kappa_mean < threshold, ]
  selected <- selected[order(selected$kappa_mean), ]
  
  if (verbose) {
    cat(sprintf("\nVariables with kappa < %.2f:\n", threshold))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    if (nrow(selected) > 0) {
      for (i in seq_len(nrow(selected))) {
        cat(sprintf("  %2d. %-20s kappa=%.3f, beta=%.4f, P(nonzero)=%.2f\n",
                    i, selected$variable[i], selected$kappa_mean[i],
                    selected$beta_mean[i], selected$prob_nonzero[i]))
      }
    } else {
      cat("  No variables meet the criterion.\n")
      cat("  Consider using select_by_credible_interval() instead.\n")
    }
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(sprintf("Total: %d relevant variables of %d\n",
                nrow(selected), nrow(coef_df)))
  }
  
  list(
    selected = selected$variable,
    details = selected,
    threshold = threshold,
    method = "kappa"
  )
}


#' Select Variables Based on Credible Intervals
#'
#' @description
#' Robust variable selection method using posterior credible intervals.
#' A variable is selected if its credible interval excludes zero, indicating
#' a statistically meaningful effect. This method is more robust than
#' kappa-based selection when the global shrinkage parameter tau collapses.
#'
#' @param hs_fit Object returned by \code{\link{fit_horseshoe}}.
#' @param prob Probability level for the credible interval. Default 0.95
#'   uses the 95% credible interval (2.5% to 97.5% quantiles).
#' @param verbose Logical for messages.
#'
#' @return List with selected variable names and details.
#'
#' @details
#' This function selects variables whose posterior credible interval does not
#' include zero. This is analogous to checking if a confidence interval
#' excludes zero in frequentist statistics, but with a Bayesian interpretation.
#'
#' The method is particularly useful when:
#' \itemize{
#'   \item The kappa-based selection returns no variables
#'   \item The posterior tau is very small (< 0.05)
#'   \item You want a more interpretable selection criterion
#' }
#'
#' @examples
#' \dontrun{
#' hs_fit <- fit_horseshoe(y, X, p0 = 5)
#' selected <- select_by_credible_interval(hs_fit, prob = 0.95, verbose = TRUE)
#' print(selected$selected)
#' }
#'
#' @seealso \code{\link{select_by_shrinkage}} for kappa-based selection.
#'
#' @export
select_by_credible_interval <- function(hs_fit, prob = 0.95, verbose = FALSE) {
  validate_horseshoe_fit(hs_fit, "select_by_credible_interval")

  coef_df <- hs_fit$coefficients
  
  alpha <- (1 - prob) / 2
  
  if (prob == 0.95) {
    lower_col <- "beta_q025"
    upper_col <- "beta_q975"
  } else if (prob == 0.90) {
    lower_col <- "beta_q10"
    upper_col <- "beta_q90"
  } else {
    if (!is.null(hs_fit$posterior_draws$beta)) {
      beta_draws <- hs_fit$posterior_draws$beta
      lower_q <- apply(beta_draws, 2, stats::quantile, probs = alpha)
      upper_q <- apply(beta_draws, 2, stats::quantile, probs = 1 - alpha)
      coef_df$ci_lower <- lower_q[coef_df$variable]
      coef_df$ci_upper <- upper_q[coef_df$variable]
      lower_col <- "ci_lower"
      upper_col <- "ci_upper"
    } else {
      warning("Custom prob requires posterior draws. Using 95% CI instead.")
      lower_col <- "beta_q025"
      upper_col <- "beta_q975"
      prob <- 0.95
    }
  }
  
  excludes_zero <- coef_df[[lower_col]] > 0 | coef_df[[upper_col]] < 0
  selected <- coef_df[excludes_zero, ]
  
  selected <- selected[order(-abs(selected$beta_mean)), ]
  
  if (verbose) {
    cat(sprintf("\nVariables with %.0f%% CI excluding zero:\n", prob * 100))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    
    if (nrow(selected) > 0) {
      for (i in seq_len(nrow(selected))) {
        cat(sprintf("  %2d. %-15s beta=%7.3f [%7.3f, %7.3f]\n",
                    i, 
                    selected$variable[i], 
                    selected$beta_mean[i],
                    selected[[lower_col]][i],
                    selected[[upper_col]][i]))
      }
    } else {
      cat("  No variables have CI excluding zero.\n")
    }
    cat(paste(rep("-", 60), collapse = ""), "\n")
    cat(sprintf("Total: %d relevant variables of %d\n",
                nrow(selected), nrow(coef_df)))
    
    tau_mean <- mean(hs_fit$posterior_draws$tau)
    if (tau_mean < 0.05) {
      cat(sprintf("\nNote: Posterior tau = %.4f (very small).\n", tau_mean))
      cat("CI-based selection is recommended over kappa-based selection.\n")
    }
  }
  
  list(
    selected = selected$variable,
    details = selected,
    prob = prob,
    method = "credible_interval"
  )
}


#' Select Variables Based on Effect Magnitude
#'
#' @description
#' Simple variable selection based on the magnitude of posterior mean
#' coefficients. Useful as a quick screening method.
#'
#' @param hs_fit Object returned by \code{\link{fit_horseshoe}}.
#' @param threshold Minimum absolute coefficient value. Default 0.1.
#' @param verbose Logical for messages.
#'
#' @return List with selected variable names and details.
#'
#' @export
select_by_magnitude <- function(hs_fit, threshold = 0.1, verbose = FALSE) {
  validate_horseshoe_fit(hs_fit, "select_by_magnitude")

  coef_df <- hs_fit$coefficients
  
  selected <- coef_df[abs(coef_df$beta_mean) > threshold, ]
  selected <- selected[order(-abs(selected$beta_mean)), ]
  
  if (verbose) {
    cat(sprintf("\nVariables with |beta| > %.3f:\n", threshold))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    if (nrow(selected) > 0) {
      for (i in seq_len(nrow(selected))) {
        cat(sprintf("  %2d. %-20s beta=%.4f\n",
                    i, selected$variable[i], selected$beta_mean[i]))
      }
    } else {
      cat("  No variables meet the criterion.\n")
    }
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(sprintf("Total: %d relevant variables of %d\n",
                nrow(selected), nrow(coef_df)))
  }
  
  list(
    selected = selected$variable,
    details = selected,
    threshold = threshold,
    method = "magnitude"
  )
}


#' @export
print.signaly_horseshoe <- function(x, ...) {
  print_horseshoe_summary(x)
}


#' @export
summary.signaly_horseshoe <- function(object, ...) {
  print_horseshoe_summary(object)
  invisible(object)
}

#' @export
print.horseshoe_fit <- function(x, ...) {
  cat("Horseshoe Regression Fit\n")
  cat(sprintf("  Predictors: %d, m_eff: %.1f\n", length(x$kappa), x$m_eff))
  invisible(x)
}