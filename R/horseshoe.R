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
#' Stan via cmdstanr. The Horseshoe prior provides adaptive shrinkage that
#' aggressively shrinks irrelevant coefficients toward zero while allowing
#' truly relevant coefficients to remain large.
#'
#' @param y Numeric vector of the response variable (target signal).
#' @param X Matrix or data frame of predictor variables (candidate signals).
#' @param var_names Optional character vector of variable names. If NULL,
#'   column names of X are used.
#' @param p0 Expected number of non-zero coefficients. If NULL, defaults to
#'   P/3 where P is the number of predictors. This controls the global
#'   shrinkage strength.
#' @param slab_scale Scale parameter for the regularizing slab. Default is 2.
#'   Larger values allow larger coefficients for selected variables.
#' @param slab_df Degrees of freedom for the regularizing slab t-distribution.
#'   Default is 4. Lower values give heavier tails.
#' @param use_qr Logical indicating whether to use QR decomposition for
#'   improved numerical stability with correlated predictors. Default TRUE.
#' @param standardize Logical indicating whether to standardize predictors
#'   internally. Results are returned on original scale. Default TRUE.
#' @param iter_warmup Number of warmup (burn-in) iterations per chain.
#'   Default 2000.
#' @param iter_sampling Number of sampling iterations per chain. Default 2000.
#' @param chains Number of MCMC chains. Default 4.
#' @param adapt_delta Target acceptance probability for HMC. Higher values
#'   reduce divergences but slow sampling. Default 0.99.
#' @param max_treedepth Maximum tree depth for NUTS sampler. Default 15.
#' @param seed Random seed for reproducibility.
#' @param verbose Logical for progress messages.
#'
#' @return A list of class "signaly_horseshoe" containing:
#' \describe{
#'   \item{coefficients}{Data frame with posterior summaries for each
#'     coefficient including mean, SD, credible intervals, shrinkage factor
#'     kappa, and relevance probabilities}
#'   \item{hyperparameters}{Data frame with posterior summaries for
#'     hyperparameters (tau, sigma, alpha, m_eff)}
#'   \item{diagnostics}{MCMC diagnostics including divergences, R-hat, ESS}
#'   \item{loo}{Leave-one-out cross-validation results}
#'   \item{posterior_draws}{Raw posterior draws for all parameters}
#'   \item{fit}{The cmdstanr fit object (if cmdstanr available)}
#'   \item{settings}{Parameters used in the analysis}
#'   \item{sparsity}{Summary of sparsity pattern}
#' }
#'
#' @details
#' The regularized Horseshoe prior (Piironen & Vehtari, 2017) models
#' coefficients as:
#' \deqn{\beta_j | \lambda_j, \tau, c \sim N(0, \tau^2 \tilde\lambda_j^2)}
#'
#' where the regularized local scale is:
#' \deqn{\tilde\lambda_j^2 = \frac{c^2 \lambda_j^2}{c^2 + \tau^2 \lambda_j^2}}
#'
#' This combines:
#' \itemize{
#'   \item **Global shrinkage** \eqn{\tau}: Controls overall sparsity, with
#'     prior calibrated to expected number of non-zero coefficients p0
#'   \item **Local shrinkage** \eqn{\lambda_j}: Half-Cauchy(0,1) allowing
#'     individual coefficients to escape shrinkage
#'   \item **Regularizing slab** c: Prevents coefficients from becoming
#'     unreasonably large for selected variables
#' }
#'
#' @section Shrinkage Factor Interpretation:
#' The shrinkage factor \eqn{\kappa_j} for each coefficient measures how much
#' it is shrunk toward zero:
#' \deqn{\kappa_j \approx \frac{1}{1 + \tau^2 \tilde\lambda_j^2}}
#'
#' \itemize{
#'   \item \eqn{\kappa_j \approx 0}: Coefficient escapes shrinkage (relevant
#'     variable)
#'   \item \eqn{\kappa_j \approx 1}: Coefficient shrunk to zero (irrelevant
#'     variable)
#'   \item \eqn{\kappa_j \approx 0.5}: Boundary case (uncertain relevance)
#' }
#'
#' @section Effective Number of Non-Zero Coefficients:
#' The model estimates m_eff, the effective number of non-zero coefficients:
#' \deqn{m_{eff} = P - \sum_{j=1}^P \kappa_j}
#'
#' This provides a data-driven estimate of the true sparsity level.
#'
#' @section Model Diagnostics:
#' The function performs comprehensive MCMC diagnostics:
#' \itemize{
#'   \item **Divergences**: Indicate geometric problems; should be 0
#'   \item **R-hat**: Chain mixing; should be < 1.01
#'   \item **ESS**: Effective sample size; should be > 400
#'   \item **BFMI**: Bayesian fraction of missing information; should be > 0.3
#' }
#'
#' @references
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of
#' Statistics, 11(2), 5018-5051. \doi{10.1214/17-EJS1337SI}
#'
#' Piironen, J., & Vehtari, A. (2017). On the hyperprior choice for the global
#' shrinkage parameter in the horseshoe prior. Proceedings of Machine Learning
#' Research, 54, 905-913.
#'
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. Biometrika, 97(2), 465-480.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), ncol = p)
#' beta_true <- c(rep(2, 3), rep(0, p - 3))
#' y <- X %*% beta_true + rnorm(n)
#' result <- fit_horseshoe(y, X, iter_warmup = 1000, iter_sampling = 1000)
#' print(result$coefficients)
#' }
#'
#' @seealso \code{\link{select_by_shrinkage}}, \code{\link{signal_analysis}}
#'
#' @export
fit_horseshoe <- function(y,
                          X,
                          var_names = NULL,
                          p0 = NULL,
                          slab_scale = 2,
                          slab_df = 4,
                          use_qr = TRUE,
                          standardize = TRUE,
                          iter_warmup = 2000,
                          iter_sampling = 2000,
                          chains = 4,
                          adapt_delta = 0.99,
                          max_treedepth = 15,
                          seed = 123,
                          verbose = FALSE) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required for Horseshoe regression. ",
      "Please install it from https://mc-stan.org/r-packages/",
      call. = FALSE
    )
  }

  cmdstan_path <- tryCatch(
    cmdstanr::cmdstan_path(),
    error = function(e) NULL
  )

  if (is.null(cmdstan_path) || nchar(cmdstan_path) == 0) {
    stop(
      "CmdStan binary is not found. ",
      "Please install it using cmdstanr::install_cmdstan()",
      call. = FALSE
    )
  }

  y <- as.numeric(y)
  X <- as.matrix(X)

  if (nrow(X) != length(y)) {
    stop("Number of rows in X must equal length of y.", call. = FALSE)
  }

  N <- length(y)
  P <- ncol(X)

  if (is.null(var_names)) {
    if (!is.null(colnames(X))) {
      var_names <- colnames(X)
    } else {
      var_names <- paste0("V", seq_len(P))
    }
  }
  colnames(X) <- var_names

  if (is.null(p0)) {
    p0 <- P / 3
    if (verbose) {
      message(sprintf(
        "Using p0 = %.1f (expecting ~%.0f%% relevant variables)",
        p0, 100 * p0 / P
      ))
    }
  }

  stan_file <- system.file("stan", "regularized_horseshoe.stan",
                           package = "SignalY")

  if (stan_file == "" || !file.exists(stan_file)) {
    stan_code <- get_horseshoe_stan_code()
    stan_file <- file.path(tempdir(), "regularized_horseshoe.stan")
    writeLines(stan_code, stan_file)
  }

  if (verbose) {
    message("============================================================")
    message("REGULARIZED HORSESHOE REGRESSION")
    message("============================================================")
    message(sprintf("N = %d observations, P = %d predictors", N, P))
    message(sprintf("p0 = %.1f, slab_scale = %.1f, slab_df = %.1f",
                    p0, slab_scale, slab_df))
    message(sprintf("QR decomposition: %s, Standardization: %s",
                    ifelse(use_qr, "Yes", "No"),
                    ifelse(standardize, "Yes", "No")))
    message(sprintf("Iterations: %d warmup + %d sampling x %d chains",
                    iter_warmup, iter_sampling, chains))
    message(sprintf("adapt_delta = %.3f, max_treedepth = %d",
                    adapt_delta, max_treedepth))
  }

  scale_global <- (p0 / (P - p0)) / sqrt(N)
  
  stan_data <- list(
    N = N,
    P = P,
    X = X,
    y = y,
    scale_global = scale_global,
    nu_global = 1,
    nu_local = 1,
    slab_scale = slab_scale,
    slab_df = slab_df
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
    parallel_chains = min(chains, parallel::detectCores() - 1),
    threads_per_chain = 1,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = ifelse(verbose, 500, 0),
    show_messages = verbose
  )

  if (verbose) message("\nExtracting diagnostics...")
  diagnostics <- extract_mcmc_diagnostics(fit, verbose = verbose)

  if (verbose) message("Processing results...")
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
    sparsity = sparsity,
    var_names = var_names,
    settings = list(
      p0 = p0,
      slab_scale = slab_scale,
      slab_df = slab_df,
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

  class(output) <- c("signaly_horseshoe", "list")

  if (verbose) {
    message("\n============================================================")
    message("SUMMARY")
    message("============================================================")
    print_horseshoe_summary(output)
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
    mean = colMeans(beta_draws),
    sd = apply(beta_draws, 2, stats::sd),
    q025 = apply(beta_draws, 2, stats::quantile, 0.025),
    q10 = apply(beta_draws, 2, stats::quantile, 0.10),
    median = apply(beta_draws, 2, stats::median),
    q90 = apply(beta_draws, 2, stats::quantile, 0.90),
    q975 = apply(beta_draws, 2, stats::quantile, 0.975),
    kappa_mean = colMeans(kappa_draws),
    kappa_median = apply(kappa_draws, 2, stats::median),
    stringsAsFactors = FALSE
  )

  threshold <- 0.01 * mean(abs(coefficients$mean[coefficients$mean != 0]))
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

  loo_result <- loo::loo(log_lik, cores = min(2, parallel::detectCores() - 1))

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
                i, top$variable[i], top$mean[i],
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
#' Returns the Stan code for the regularized Horseshoe model. This is used
#' when the package is not installed and the inst/stan file is not available.
#'
#' @return Character string containing Stan code.
#'
#' @keywords internal
get_horseshoe_stan_code <- function() {
  '
functions {
  vector regularized_hs_scale(vector lambda, real tau, real c2) {
    int P = rows(lambda);
    vector[P] lambda2 = square(lambda);
    vector[P] tau2_lambda2 = square(tau) * lambda2;
    return sqrt(c2 * tau2_lambda2 ./ (c2 + tau2_lambda2));
  }

  real effective_nonzero(vector lambda_tilde, real tau, real sigma, int n) {
    vector[rows(lambda_tilde)] kappa = 1.0 ./ (1.0 + n * square(lambda_tilde * tau / sigma));
    return sum(1 - kappa);
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  matrix[N, P] X;
  vector[N] y;
  real<lower=0> p0;
  real<lower=0> slab_scale;
  real<lower=1> slab_df;
  int<lower=0, upper=1> use_qr;
  int<lower=0, upper=1> standardize;
  int<lower=0> N_new;
  matrix[N_new > 0 ? N_new : 1, P] X_new;
}

transformed data {
  matrix[N, P] X_work;
  matrix[N, P] Q_ast;
  matrix[P, P] R_ast;
  matrix[P, P] R_ast_inv;
  vector[P] X_mean = rep_vector(0, P);
  vector[P] X_sd = rep_vector(1, P);
  real y_mean = 0;
  real y_sd = 1;
  vector[N] y_work;
  real sigma_y_approx = sd(y);
  real tau0 = (p0 / (P - p0 + 1e-8)) * sigma_y_approx / sqrt(N);

  if (standardize == 1) {
    for (j in 1:P) {
      X_mean[j] = mean(X[, j]);
      X_sd[j] = sd(X[, j]);
      if (X_sd[j] < 1e-10) X_sd[j] = 1;
    }
    y_mean = mean(y);
    y_sd = sd(y);
    if (y_sd < 1e-10) y_sd = 1;
    for (j in 1:P) {
      X_work[, j] = (X[, j] - X_mean[j]) / X_sd[j];
    }
    y_work = (y - y_mean) / y_sd;
  } else {
    X_work = X;
    y_work = y;
  }

  if (use_qr == 1) {
    Q_ast = qr_thin_Q(X_work) * sqrt(N - 1.0);
    R_ast = qr_thin_R(X_work) / sqrt(N - 1.0);
    R_ast_inv = inverse(R_ast);
  } else {
    Q_ast = X_work;
    R_ast = diag_matrix(rep_vector(1, P));
    R_ast_inv = diag_matrix(rep_vector(1, P));
  }
}

parameters {
  real alpha_raw;
  vector[P] beta_raw;
  real<lower=0> tau;
  vector<lower=0>[P] lambda;
  real<lower=0> c2_tilde;
  vector<lower=0>[P] lambda_aux;
  real<lower=0> tau_aux;
  real<lower=0> sigma;
}

transformed parameters {
  vector[P] beta;
  vector[P] lambda_tilde;
  real c2;
  vector[N] mu;
  real alpha;

  c2 = slab_scale^2 * c2_tilde;
  lambda_tilde = regularized_hs_scale(lambda, tau, c2);

  if (use_qr == 1) {
    vector[P] theta = tau * lambda_tilde .* beta_raw;
    beta = R_ast_inv * theta;
  } else {
    beta = tau * lambda_tilde .* beta_raw;
  }

  alpha = alpha_raw * sigma_y_approx;

  if (use_qr == 1) {
    mu = alpha + Q_ast * (tau * lambda_tilde .* beta_raw);
  } else {
    mu = alpha + X_work * beta;
  }
}

model {
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();
  tau_aux ~ cauchy(0, 1);
  tau ~ cauchy(0, tau0);
  lambda_aux ~ cauchy(0, 1);
  lambda ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  sigma ~ cauchy(0, sigma_y_approx);
  y_work ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N_new > 0 ? N_new : 1] y_new_rep;
  vector[P] kappa;
  real m_eff;
  vector[P] beta;
  real alpha;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y_work[n] | mu[n], sigma);
  }

  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }

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
    mu_new = alpha + X_new_work * beta;
    for (n in 1:N_new) {
      y_new_rep[n] = normal_rng(mu_new[n], sigma);
    }
  } else {
    y_new_rep[1] = 0;
  }

  for (j in 1:P) {
    kappa[j] = 1.0 / (1.0 + square(tau * lambda_tilde[j]));
  }

  m_eff = P - sum(kappa);

  if (standardize == 1) {
    for (j in 1:P) {
      beta[j] = beta[j] * y_sd / X_sd[j];
    }
    alpha = y_mean + alpha * y_sd - dot_product(beta, X_mean);
  } else {
    beta = beta;
    alpha = alpha;
  }
}
'
}


#' Select Variables Based on Shrinkage
#'
#' @description
#' Alternative variable selection method using shrinkage factors (kappa)
#' directly. Does not require projpred.
#'
#' @param hs_fit Object returned by \code{\link{fit_horseshoe}}.
#' @param threshold Kappa threshold. Variables with kappa < threshold are
#'   considered relevant. Default 0.5.
#' @param verbose Logical for messages.
#'
#' @return Character vector of selected variable names.
#'
#' @export
select_by_shrinkage <- function(hs_fit, threshold = 0.5, verbose = FALSE) {

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
                    selected$mean[i], selected$prob_nonzero[i]))
      }
    } else {
      cat("  No variables meet the criterion.\n")
      cat("  All variables were shrunk (kappa near 1).\n")
    }
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(sprintf("Total: %d relevant variables of %d\n",
                nrow(selected), nrow(coef_df)))
  }

  list(
    selected = selected$variable,
    details = selected,
    threshold = threshold
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