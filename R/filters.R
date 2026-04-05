#' @title Signal Filtering Methods for Trend Extraction
#' @description
#' This module implements three complementary spectral decomposition methods
#' for extracting latent trend signals from time series: wavelet multiresolution
#' analysis, empirical mode decomposition, and the Grant-Chan embedded
#' Hodrick-Prescott filter.
#'
#' @name filters
NULL


#' Wavelet Multiresolution Analysis Filter
#'
#' @description
#' Performs wavelet-based signal decomposition using the Maximal Overlap
#' Discrete Wavelet Transform (MODWT) to extract trend components at specified
#' frequency scales. This method decomposes the signal into detail coefficients
#' (D1, D2, ..., DJ) capturing progressively lower frequencies and a smooth
#' coefficient (SJ) representing the underlying trend.
#'
#' @param y Numeric vector of the time series to decompose. Length must be at
#'   least \code{2^J}.
#' @param wf Character string specifying the wavelet filter. Options include
#'   "la8" (least asymmetric with 8 vanishing moments, 16 coefficients),
#'   "la16", "la20", "haar", "d4", "d6", "d8", etc. Default is "la8".
#' @param J Integer specifying the decomposition depth (number of levels).
#'   Default is 4, yielding D1-D4 detail levels plus S4 smooth level.
#' @param boundary Character string specifying boundary handling: "periodic"
#'   (default) or "reflection".
#' @param levels_to_combine Integer vector specifying which detail levels to
#'   combine for the trend estimate. Default is \code{c(3, 4)} for D3+D4.
#' @param first_difference Logical. If TRUE, applies wavelet to first
#'   differences and reconstructs via cumulative sum. Default is FALSE.
#' @param verbose Logical indicating whether to print diagnostic messages.
#'
#' @return A list of class "signaly_wavelet" containing:
#' \describe{
#'   \item{trend}{Numeric vector of the extracted trend component}
#'   \item{mra}{Full multiresolution analysis object from waveslim::mra}
#'   \item{detail_levels}{Data frame with all detail level coefficients}
#'   \item{smooth_level}{Vector of the smooth (SJ) coefficients}
#'   \item{combined_levels}{Character string indicating which levels were
#'     combined}
#'   \item{settings}{List of parameters used in the analysis}
#'   \item{diagnostics}{List with variance decomposition and energy
#'     distribution}
#' }
#'
#' @details
#' The MODWT (Maximal Overlap Discrete Wavelet Transform) is preferred over
#' the classical DWT for several reasons relevant to signal extraction:
#'
#' \enumerate{
#'   \item **Translation invariance**: Unlike DWT, MODWT does not depend on
#'     the starting point of the series, producing consistent results regardless
#'     of circular shifts.
#'   \item **Any sample size**: MODWT can be applied to series of any length,
#'     not just powers of 2.
#'   \item **Additive decomposition**: The MRA (multiresolution analysis)
#'     coefficients sum exactly to the original series.
#' }
#'
#' The choice of wavelet filter affects the trade-off between time and frequency
#' localization:
#' \itemize{
#'   \item **la8 (Daubechies least asymmetric, 8 vanishing moments)**: Good
#'     balance of smoothness and localization, recommended for economic data.
#'   \item **Higher order (la16, la20)**: Better frequency resolution at cost
#'     of temporal smearing.
#'   \item **haar**: Maximum time localization but poor frequency resolution.
#' }
#'
#' @section Frequency Interpretation:
#' For a series with unit sampling interval, the detail levels correspond to
#' approximate frequency bands:
#' \itemize{
#'   \item D1: periods 2-4 (highest frequency noise)
#'   \item D2: periods 4-8 (short-term fluctuations)
#'   \item D3: periods 8-16 (medium-term cycles)
#'   \item D4: periods 16-32 (longer cycles)
#'   \item S4: periods > 32 (smooth trend)
#' }
#'
#' For annual economic data, D3+D4 typically captures business cycle dynamics
#' (8-32 year periods), while D1+D2 captures short-term noise.
#'
#' @references
#' Daubechies, I. (1992). Ten Lectures on Wavelets. SIAM.
#'
#' Percival, D. B., & Walden, A. T. (2000). Wavelet Methods for Time Series
#' Analysis. Cambridge University Press.
#'
#' Gencay, R., Selcuk, F., & Whitcher, B. (2002). An Introduction to Wavelets
#' and Other Filtering Methods in Finance and Economics. Academic Press.
#'
#' @examples
#' set.seed(123)
#' y <- cumsum(rnorm(100)) + sin(seq(0, 4*pi, length.out = 100))
#' result <- filter_wavelet(y, wf = "la8", J = 4)
#' plot(y, type = "l", col = "gray")
#' lines(result$trend, col = "red", lwd = 2)
#'
#' @seealso \code{\link[waveslim]{mra}}, \code{\link{filter_emd}},
#'   \code{\link{filter_hpgc}}
#'
#' @export
filter_wavelet <- function(y,
                           wf = "la8",
                           J = 4,
                           boundary = "periodic",
                           levels_to_combine = c(3, 4),
                           first_difference = FALSE,
                           verbose = FALSE) {
  
  if (!requireNamespace("waveslim", quietly = TRUE)) {
    stop("Package 'waveslim' is required. Please install it.", call. = FALSE)
  }
  
  y <- as.numeric(y)
  n <- length(y)
  
  if (n < 2^J) {
    stop(sprintf(
      "Series length (%d) must be at least 2^J = %d for J = %d.",
      n, 2^J, J
    ), call. = FALSE)
  }
  
  if (any(is.na(y))) {
    stop("Input series contains NA values. Please handle missing data first.",
         call. = FALSE)
  }
  
  if (first_difference) {
    y_work <- diff(y)
    if (verbose) message("Applied first differencing before wavelet decomposition.")
  } else {
    y_work <- y
  }
  
  if (verbose) {
    message(sprintf("Applying MODWT with %s filter, J = %d levels...", wf, J))
  }
  
  mra_result <- waveslim::mra(
    y_work,
    wf = wf,
    J = J,
    method = "modwt",
    boundary = boundary
  )
  
  component_names <- names(mra_result)
  detail_names <- grep("^[Dd]\\d+$", component_names, value = TRUE)
  smooth_name <- grep("^[Ss]\\d+$", component_names, value = TRUE)
  
  detail_df <- as.data.frame(mra_result[detail_names])
  names(detail_df) <- toupper(names(detail_df))
  
  smooth_vec <- as.numeric(mra_result[[smooth_name]])
  
  valid_levels <- seq_len(J)
  levels_to_combine <- levels_to_combine[levels_to_combine %in% valid_levels]
  if (length(levels_to_combine) == 0) {
    stop("No valid detail levels specified in 'levels_to_combine'.", call. = FALSE)
  }
  
  combined_names <- paste0("D", levels_to_combine)
  if (!all(combined_names %in% names(detail_df))) {
    available <- names(detail_df)
    stop(sprintf(
      "Requested levels %s not all available. Available: %s",
      paste(combined_names, collapse = ", "),
      paste(available, collapse = ", ")
    ), call. = FALSE)
  }
  
  trend_diff <- rowSums(detail_df[, combined_names, drop = FALSE])
  
  if (first_difference) {
    trend <- c(y[1], y[1] + cumsum(trend_diff))
    if (length(trend) > n) trend <- trend[seq_len(n)]
    if (length(trend) < n) trend <- c(trend, rep(trend[length(trend)], n - length(trend)))
  } else {
    trend <- trend_diff
  }
  
  total_var <- stats::var(y_work)
  var_by_level <- sapply(mra_result, stats::var)
  var_proportion <- var_by_level / sum(var_by_level)
  
  diagnostics <- list(
    total_variance = total_var,
    variance_by_level = var_by_level,
    variance_proportion = var_proportion,
    trend_variance_ratio = stats::var(trend_diff) / total_var
  )
  
  if (verbose) {
    message("Variance decomposition by level:")
    for (i in seq_along(var_proportion)) {
      message(sprintf("  %s: %.1f%%", names(var_proportion)[i],
                      100 * var_proportion[i]))
    }
  }
  
  result <- list(
    trend = trend,
    combined = trend_diff,              
    levels = levels_to_combine,         
    filter = wf,                        
    mra = mra_result,
    detail_levels = detail_df,
    smooth_level = smooth_vec,
    combined_levels = paste(combined_names, collapse = "+"),
    settings = list(
      wf = wf,
      J = J,
      boundary = boundary,
      levels_to_combine = levels_to_combine,
      first_difference = first_difference
    ),
    diagnostics = diagnostics
  )
  
  class(result) <- c("signaly_wavelet", "list")
  result
}


#' Empirical Mode Decomposition Filter
#'
#' @description
#' Applies Empirical Mode Decomposition (EMD) to extract intrinsic mode
#' functions (IMFs) from a time series. Unlike Fourier or wavelet methods,
#' EMD is fully data-adaptive and does not require pre-specified basis
#' functions, making it suitable for non-stationary and non-linear signals.
#'
#' @param y Numeric vector of the time series to decompose.
#' @param boundary Character string specifying boundary handling: "periodic"
#'   (default), "symmetric", "none", or "wave".
#' @param max_imf Maximum number of IMFs to extract. If NULL, extraction
#'   continues until the residue is monotonic.
#' @param stop_rule Character string specifying the stopping criterion for
#'   sifting: "type1" (default), "type2", "type3", "type4", or "type5".
#' @param tol Tolerance for sifting convergence. Default is
#'   \code{sd(y) * 0.1^2}.
#' @param max_sift Maximum number of sifting iterations per IMF. Default is 20.
#' @param verbose Logical indicating whether to print diagnostic messages.
#'
#' @return A list of class "signaly_emd" containing:
#' \describe{
#'   \item{trend}{Numeric vector of the extracted trend (original minus
#'     residue)}
#'   \item{residue}{Numeric vector of the EMD residue (monotonic trend)}
#'   \item{imfs}{Matrix where each column is an IMF, ordered from highest to
#'     lowest frequency}
#'   \item{n_imfs}{Number of IMFs extracted}
#'   \item{original}{Original input series}
#'   \item{settings}{List of parameters used}
#'   \item{diagnostics}{List with IMF statistics}
#' }
#'
#' @details
#' EMD decomposes a signal x(t) into a sum of Intrinsic Mode Functions (IMFs)
#' and a residue:
#' \deqn{x(t) = \sum_{j=1}^{n} c_j(t) + r_n(t)}
#'
#' where each IMF \eqn{c_j(t)} satisfies two conditions:
#' \enumerate{
#'   \item The number of extrema and zero crossings differ by at most one
#'   \item The mean of upper and lower envelopes is zero at each point
#' }
#'
#' The sifting process iteratively extracts IMFs from highest to lowest
#' frequency until the residue becomes monotonic (representing the trend).
#'
#' @section Advantages over Fourier/Wavelet Methods:
#' \itemize{
#'   \item **Adaptive basis**: IMFs are derived from the data itself, not
#'     pre-specified
#'   \item **Handles non-stationarity**: Instantaneous frequency can vary
#'     over time
#'   \item **Handles non-linearity**: No assumption of linear superposition
#'   \item **Preserves local structure**: Better time localization than
#'     Fourier methods
#' }
#'
#' @section Limitations:
#' \itemize{
#'   \item **Mode mixing**: Different scales may appear in the same IMF
#'   \item **End effects**: Boundary conditions can cause artifacts
#'   \item **No formal theory**: Unlike wavelets, lacks rigorous mathematical
#'     foundation
#'   \item **Reproducibility**: Results can vary with stopping criteria
#' }
#'
#' @references
#' Huang, N. E., Shen, Z., Long, S. R., Wu, M. C., Shih, H. H., Zheng, Q.,
#' Yen, N.-C., Tung, C. C., & Liu, H. H. (1998). The empirical mode
#' decomposition and the Hilbert spectrum for nonlinear and non-stationary
#' time series analysis. Proceedings of the Royal Society A, 454(1971), 903-995.
#'
#' Wu, Z., & Huang, N. E. (2009). Ensemble empirical mode decomposition: A
#' noise-assisted data analysis method. Advances in Adaptive Data Analysis,
#' 1(1), 1-41.
#'
#' @examples
#' set.seed(123)
#' t <- seq(0, 10, length.out = 200)
#' y <- sin(2*pi*t) + 0.5*sin(8*pi*t) + 0.1*rnorm(200)
#' result <- filter_emd(y)
#' plot(y, type = "l", col = "gray")
#' lines(result$trend, col = "red", lwd = 2)
#'
#' @seealso \code{\link[EMD]{emd}}, \code{\link{filter_wavelet}},
#'   \code{\link{filter_hpgc}}
#'
#' @export
filter_emd <- function(y,
                       boundary = "periodic",
                       max_imf = NULL,
                       stop_rule = "type1",
                       tol = NULL,
                       max_sift = 20,
                       verbose = FALSE) {
  
  if (!requireNamespace("EMD", quietly = TRUE)) {
    stop("Package 'EMD' is required. Please install it.", call. = FALSE)
  }
  
  y <- as.numeric(y)
  n <- length(y)
  
  if (any(is.na(y))) {
    stop("Input series contains NA values. Please handle missing data first.",
         call. = FALSE)
  }
  
  if (n < 10) {
    stop("Series too short for EMD decomposition (minimum 10 observations).",
         call. = FALSE)
  }
  
  if (is.null(tol)) {
    tol <- stats::sd(y) * 0.01
  }
  
  if (verbose) {
    message(sprintf("Applying EMD with boundary = '%s'...", boundary))
  }
  
  emd_args <- list(
    xt = y,
    boundary = boundary,
    max.sift = max_sift
  )
  
  if (!is.null(max_imf)) {
    emd_args$max.imf <- max_imf
  }
  
  emd_result <- do.call(EMD::emd, emd_args)
  
  if (is.null(emd_result$imf)) {
    warning("EMD did not produce any IMFs. Returning original series as trend.")
    result <- list(
      trend = y,
      residue = y,
      imfs = matrix(nrow = n, ncol = 0),
      n_imfs = 0,
      original = y,
      settings = list(boundary = boundary, max_imf = max_imf, max_sift = max_sift),
      diagnostics = list()
    )
    class(result) <- c("signaly_emd", "list")
    return(result)
  }
  
  imfs <- as.matrix(emd_result$imf)
  residue <- as.numeric(emd_result$residue)
  n_imfs <- ncol(imfs)
  
  trend <- y - residue
  
  imf_vars <- apply(imfs, 2, stats::var)
  imf_energy <- apply(imfs^2, 2, sum)
  total_energy <- sum(y^2)
  
  diagnostics <- list(
    imf_variance = imf_vars,
    imf_energy = imf_energy,
    energy_proportion = imf_energy / total_energy,
    residue_variance = stats::var(residue),
    n_extrema = emd_result$nimf
  )
  
  if (verbose) {
    message(sprintf("Extracted %d IMFs.", n_imfs))
    message("Energy distribution:")
    for (i in seq_len(n_imfs)) {
      message(sprintf("  IMF%d: %.1f%%", i, 100 * diagnostics$energy_proportion[i]))
    }
    message(sprintf("  Residue: %.1f%%", 100 * sum(residue^2) / total_energy))
  }
  
  result <- list(
    trend = trend,
    residue = residue,
    imfs = imfs,
    n_imfs = n_imfs,
    original = y,
    settings = list(
      boundary = boundary,
      max_imf = max_imf,
      stop_rule = stop_rule,
      tol = tol,
      max_sift = max_sift
    ),
    diagnostics = diagnostics
  )
  
  class(result) <- c("signaly_emd", "list")
  result
}


#' Grant-Chan Embedded Hodrick-Prescott Filter
#'
#' @description
#' Implements the Bayesian Hodrick-Prescott filter embedded in an unobserved
#' components model, as developed by Grant and Chan (2017). This approach
#' provides principled uncertainty quantification for the extracted trend
#' through Markov Chain Monte Carlo sampling.
#'
#' @param y Numeric vector of the time series. Will be internally scaled
#'   for numerical stability.
#' @param prior_config Character string or list specifying prior configuration.
#'   Options: "weak" (default), "informative", or "empirical". Alternatively,
#'   a named list with prior parameters (see Details).
#' @param n_chains Integer number of MCMC chains to run. Default is 4.
#' @param iterations Integer total number of MCMC iterations per chain.
#'   Default is 20000.
#' @param burnin Integer number of burn-in iterations to discard. Default is
#'   5000.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @return A list of class "signaly_hpgc" containing:
#' \itemize{
#'   \item \code{trend}: Numeric vector of posterior mean trend
#'   \item \code{trend_lower}: Numeric vector of 2.5 percent posterior quantile
#'   \item \code{trend_upper}: Numeric vector of 97.5 percent posterior quantile
#'   \item \code{cycle}: Numeric vector of posterior mean cycle component
#'   \item \code{cycle_lower}: Numeric vector of 2.5 percent posterior quantile
#'   \item \code{cycle_upper}: Numeric vector of 97.5 percent posterior quantile
#'   \item \code{draws}: List of posterior draws for all parameters
#'   \item \code{diagnostics}: Convergence diagnostics including R-hat and ESS
#'   \item \code{dic}: Deviance Information Criterion
#'   \item \code{settings}: Parameters used in the analysis
#' }
#'
#' @details
#' The Grant-Chan model decomposes the observed series \eqn{y_t} as:
#' \deqn{y_t = \tau_t + c_t}
#'
#' where \eqn{\tau_t} is the trend component and \eqn{c_t} is the cyclical
#' component.
#'
#' **Trend Model (Second-Order Markov Process)**:
#' \deqn{\Delta^2 \tau_t = u_t^\tau, \quad u_t^\tau \sim N(0, \sigma_\tau^2)}
#'
#' This implies the trend growth rate follows a random walk, allowing for
#' time-varying trend growth.
#'
#' **Cycle Model (Stationary AR(2))**:
#' \deqn{c_t = \phi_1 c_{t-1} + \phi_2 c_{t-2} + u_t^c, \quad u_t^c \sim N(0, \sigma_c^2)}
#'
#' with stationarity constraints on \eqn{\phi}.
#'
#' @section Prior Configurations:
#' \describe{
#'   \item{weak}{Diffuse priors allowing data to dominate. Good for initial
#'     exploration.}
#'   \item{informative}{Tighter priors based on typical macroeconomic
#'     dynamics. Suitable when strong smoothness is desired.}
#'   \item{empirical}{Priors calibrated from data moments. Balances
#'     flexibility with data-driven regularization.}
#' }
#'
#' Custom priors can be specified as a list with elements:
#' \itemize{
#'   \item \code{phi_mu}: Mean of phi prior (2-vector)
#'   \item \code{phi_v_i}: Precision matrix for phi prior (2x2)
#'   \item \code{gamma_mu}: Mean of gamma (initial trend growth) prior
#'   \item \code{gamma_v_i}: Precision matrix for gamma prior
#'   \item \code{s_tau}: Upper bound for uniform prior on \eqn{\sigma_\tau^2}
#'   \item \code{s_c_shape}: Shape parameter for inverse-gamma prior on
#'     \eqn{\sigma_c^2}
#'   \item \code{s_c_rate}: Rate parameter for inverse-gamma prior on
#'     \eqn{\sigma_c^2}
#' }
#'
#' @section Relationship to Standard HP Filter:
#' The standard HP filter solves:
#' \deqn{\min_\tau \sum_t (y_t - \tau_t)^2 + \lambda \sum_t (\Delta^2 \tau_t)^2}
#'
#' The Grant-Chan approach embeds this within a probabilistic model where
#' \eqn{\lambda = \sigma_c^2 / \sigma_\tau^2}, allowing this ratio to be
#' estimated from data with full uncertainty quantification.
#'
#' @references
#' Grant, A. L., & Chan, J. C. C. (2017). Reconciling output gaps: Unobserved
#' components model and Hodrick-Prescott filter. Journal of Economic Dynamics
#' and Control, 75, 114-121. \doi{10.1016/j.jedc.2016.12.007}
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). Bayesian
#' Econometric Methods (2nd ed.). Cambridge University Press.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' y <- cumsum(rnorm(100)) + sin(seq(0, 4*pi, length.out = 100))
#' result <- filter_hpgc(y, prior_config = "weak", n_chains = 2,
#'                       iterations = 5000, burnin = 1000)
#' plot(y, type = "l", col = "gray")
#' lines(result$trend, col = "red", lwd = 2)
#' }
#'
#' @seealso \code{\link{filter_wavelet}}, \code{\link{filter_emd}}
#'
#' @export
filter_hpgc <- function(y,
                        prior_config = "weak",
                        n_chains = 4,
                        iterations = 20000,
                        burnin = 5000,
                        verbose = FALSE) {
  
  y_original <- as.numeric(y)
  n <- length(y_original)
  
  if (any(is.na(y_original))) {
    stop("Input series contains NA values. Please handle missing data first.",
         call. = FALSE)
  }
  
  if (n < 15) {
    stop("Series too short for HP-GC filter (minimum 15 observations).",
         call. = FALSE)
  }
  
  y_scaled <- y_original * 100
  p <- 2
  
  if (is.character(prior_config)) {
    prior_config <- match.arg(prior_config, c("weak", "informative", "empirical"))
    
    prior_configs <- list(
      weak = list(
        name = "Weak Priors",
        phi_mu = matrix(c(1.3, -0.7)),
        phi_v_i = diag(10, p),
        gamma_mu = matrix(rep(mean(y_scaled), 2)),
        gamma_v_i = diag(1, p),
        s_tau = 0.1,
        s_c_shape = 1,
        s_c_rate = 1
      ),
      informative = list(
        name = "Informative Priors",
        phi_mu = matrix(c(1.3, -0.7)),
        phi_v_i = diag(1, p),
        gamma_mu = matrix(c(y_scaled[1], y_scaled[1])),
        gamma_v_i = diag(1 / 100, p),
        s_tau = 0.01,
        s_c_shape = 3,
        s_c_rate = 2
      ),
      empirical = list(
        name = "Empirical Priors",
        phi_mu = matrix(c(1.3, -0.7)),
        phi_v_i = diag(2, p),
        gamma_mu = matrix(rep(mean(diff(y_scaled)), 2)),
        gamma_v_i = diag(stats::var(diff(y_scaled)) / 10, p),
        s_tau = stats::sd(y_scaled) / 100,
        s_c_shape = 2,
        s_c_rate = 1.5
      )
    )
    
    prior <- prior_configs[[prior_config]]
  } else if (is.list(prior_config)) {
    required_elements <- c("phi_mu", "phi_v_i", "gamma_mu", "gamma_v_i",
                           "s_tau", "s_c_shape", "s_c_rate")
    missing <- setdiff(required_elements, names(prior_config))
    if (length(missing) > 0) {
      stop(sprintf("Prior config missing elements: %s", paste(missing, collapse = ", ")),
           call. = FALSE)
    }
    prior <- prior_config
    if (is.null(prior$name)) prior$name <- "Custom"
  } else {
    stop("prior_config must be 'weak', 'informative', 'empirical', or a list.",
         call. = FALSE)
  }
  
  if (verbose) {
    message(sprintf("Running HP-GC filter with %s...", prior$name))
    message(sprintf("Chains: %d, Iterations: %d, Burn-in: %d", n_chains, iterations, burnin))
  }
  
  run_single_chain <- function(chain_id) {
    
    x_gamma <- cbind(2:(n + 1), -1:-n)
    h2 <- diag(1, n)
    diag(h2[-1, -n]) <- -2
    diag(h2[-(1:2), -((n - 1):n)]) <- 1
    h2h2 <- crossprod(h2)
    
    h_phi <- diag(1, n)
    phi <- matrix(c(1.34 + stats::rnorm(1, 0, 0.1), -0.7 + stats::rnorm(1, 0, 0.1)))
    
    st_check <- check_stationarity(phi)
    if (!st_check$is_stationary) {
      phi <- matrix(c(1.3, -0.7))
    }
    
    for (i in 1:p) {
      diag(h_phi[-(1:i), -((n - i):n)]) <- -phi[i, ]
    }
    
    s_tau_i <- 1 / stats::runif(1, 0.0001, 0.01)
    s_c_i <- 1 / stats::runif(1, 0.1, 1)
    gamma <- t(rep(y_scaled[1] + stats::rnorm(1, 0, stats::sd(y_scaled) / 10), 2))
    
    n_draws <- iterations - burnin
    draws_tau <- matrix(NA_real_, n, n_draws)
    draws_c <- matrix(NA_real_, n, n_draws)
    draws_phi <- matrix(NA_real_, p, n_draws)
    draws_gamma <- matrix(NA_real_, p, n_draws)
    draws_s_tau_i <- numeric(n_draws)
    draws_s_c_i <- numeric(n_draws)
    
    phi_accepted <- 0
    phi_proposed <- 0
    
    for (draw in seq_len(iterations)) {
      
      alpha <- solve(h2, matrix(c(2 * gamma[1] - gamma[2], -gamma[1], rep(0, n - 2))))
      sh2 <- s_tau_i * h2h2
      shphi <- s_c_i * as.matrix(crossprod(h_phi))
      K_tau <- sh2 + shphi
      mu_tau <- solve(K_tau, sh2 %*% alpha + shphi %*% y_scaled)
      tau <- as.vector(mu_tau + solve(chol(K_tau), stats::rnorm(n)))
      
      cvec <- c(rep(0, p), y_scaled - tau)
      temp <- stats::embed(cvec, 1 + p)
      c_now <- matrix(temp[, 1])
      x_phi <- temp[, -1, drop = FALSE]
      K_phi <- prior$phi_v_i + s_c_i * crossprod(x_phi)
      mu_phi <- solve(K_phi, prior$phi_v_i %*% prior$phi_mu + s_c_i * crossprod(x_phi, c_now))
      phi_can <- mu_phi + solve(chol(K_phi), stats::rnorm(p))
      
      phi_proposed <- phi_proposed + 1
      st <- check_stationarity(phi_can)
      if (st$is_stationary && sum(phi_can) < 0.99 && phi_can[2] - phi_can[1] < 0.99 && phi_can[2] > -0.99) {
        phi <- phi_can
        for (i in 1:p) {
          diag(h_phi[-(1:i), -((n - i):n)]) <- -phi[i, ]
        }
        phi_accepted <- phi_accepted + 1
      }
      
      s_c_i <- stats::rgamma(1, shape = prior$s_c_shape + n / 2,
                             rate = prior$s_c_rate + crossprod(c_now - x_phi %*% phi) / 2)
      
      tausq_sum <- sum(diff(diff(c(gamma[2:1], tau)))^2)
      s_tau_can <- seq(from = max(1e-6, stats::runif(1) / 1000),
                       to = prior$s_tau - stats::runif(1) / 1000,
                       length.out = 400)
      lik <- -n / 2 * log(s_tau_can) - tausq_sum / (2 * s_tau_can)
      lik <- lik - max(lik)
      plik <- exp(lik)
      plik <- plik / sum(plik)
      plik <- cumsum(plik)
      s_tau_i <- 1 / s_tau_can[which(stats::runif(1) < plik)[1]]
      
      sxh2 <- s_tau_i * crossprod(x_gamma, h2h2)
      K_gamma <- as.matrix(prior$gamma_v_i + sxh2 %*% x_gamma)
      mu_gamma <- solve(K_gamma, prior$gamma_v_i %*% prior$gamma_mu + sxh2 %*% tau)
      gamma <- as.vector(mu_gamma + solve(chol(K_gamma), stats::rnorm(2)))
      
      if (draw > burnin) {
        k <- draw - burnin
        draws_tau[, k] <- tau
        draws_c[, k] <- c_now
        draws_phi[, k] <- phi
        draws_gamma[, k] <- gamma
        draws_s_tau_i[k] <- s_tau_i
        draws_s_c_i[k] <- s_c_i
      }
    }
    
    list(
      tau = draws_tau,
      c = draws_c,
      phi = draws_phi,
      gamma = draws_gamma,
      s_tau_i = draws_s_tau_i,
      s_c_i = draws_s_c_i,
      acceptance_rate_phi = phi_accepted / max(1, phi_proposed)
    )
  }
  
  chains <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    if (verbose) message(sprintf("Running chain %d/%d...", ch, n_chains))
    chains[[ch]] <- run_single_chain(ch)
  }
  
  all_tau <- do.call(cbind, lapply(chains, function(x) x$tau))
  all_c <- do.call(cbind, lapply(chains, function(x) x$c))
  all_phi <- do.call(cbind, lapply(chains, function(x) x$phi))
  all_s_c_i <- unlist(lapply(chains, function(x) x$s_c_i))
  
  tau_mean <- rowMeans(all_tau) / 100
  tau_q025 <- apply(all_tau, 1, stats::quantile, 0.025) / 100
  tau_q975 <- apply(all_tau, 1, stats::quantile, 0.975) / 100
  
  c_mean <- rowMeans(all_c) / 100
  c_q025 <- apply(all_c, 1, stats::quantile, 0.025) / 100
  c_q975 <- apply(all_c, 1, stats::quantile, 0.975) / 100
  
  calculate_dic <- function() {
    n_total <- ncol(all_tau)
    dev <- numeric(n_total)
    for (i in seq_len(n_total)) {
      resid <- y_scaled - all_tau[, i] - all_c[, i]
      dev[i] <- -2 * sum(stats::dnorm(resid, 0, sqrt(1 / all_s_c_i[i]), log = TRUE))
    }
    D_bar <- mean(dev)
    theta_bar <- rowMeans(all_tau + all_c)
    rb <- y_scaled - theta_bar
    D_theta_bar <- -2 * sum(stats::dnorm(rb, 0, stats::sd(rb), log = TRUE))
    pD <- D_bar - D_theta_bar
    list(DIC = D_bar + pD, pD = pD, D_bar = D_bar)
  }
  
  dic <- calculate_dic()
  
  acceptance_rates <- sapply(chains, function(x) x$acceptance_rate_phi)
  
  diagnostics <- list(
    acceptance_rates = acceptance_rates,
    mean_acceptance = mean(acceptance_rates),
    n_chains = n_chains,
    n_draws_per_chain = iterations - burnin,
    total_draws = n_chains * (iterations - burnin)
  )
  
  if (verbose) {
    message(sprintf("DIC: %.2f (pD = %.2f)", dic$DIC, dic$pD))
    message(sprintf("Mean phi acceptance rate: %.3f", mean(acceptance_rates)))
  }
  
  result <- list(
    trend = tau_mean,
    trend_lower = tau_q025,
    trend_upper = tau_q975,
    cycle = c_mean,
    cycle_lower = c_q025,
    cycle_upper = c_q975,
    draws = list(
      tau = all_tau / 100,
      c = all_c / 100,
      phi = all_phi,
      s_c_i = all_s_c_i
    ),
    diagnostics = diagnostics,
    dic = dic,
    settings = list(
      prior_name = prior$name,
      n_chains = n_chains,
      iterations = iterations,
      burnin = burnin
    )
  )
  
  class(result) <- c("signaly_hpgc", "list")
  result
}


#' Apply Multiple Filters to a Series
#'
#' @description
#' Convenience function that applies all three filtering methods (wavelet, EMD,
#' HP-GC) to a time series and returns a consolidated comparison of results.
#'
#' @param y Numeric vector of the time series.
#' @param wavelet_wf Wavelet filter for wavelet decomposition. Default "la8".
#' @param wavelet_J Wavelet decomposition depth. Default 4.
#' @param wavelet_levels Levels to combine for wavelet trend. Default c(3, 4).
#' @param hpgc_prior Prior configuration for HP-GC. Default "weak".
#' @param hpgc_chains Number of MCMC chains. Default 4.
#' @param hpgc_iterations MCMC iterations. Default 20000.
#' @param hpgc_burnin MCMC burn-in. Default 5000.
#' @param verbose Logical for progress messages.
#'
#' @return A list of class "signaly_multifilter" containing results from all
#'   three methods and a comparison data frame.
#'
#' @examples
#' \donttest{
#' y <- cumsum(rnorm(100)) + sin(seq(0, 4*pi, length.out = 100))
#' result <- filter_all(y, hpgc_iterations = 5000, hpgc_burnin = 1000)
#' }
#'
#' @export
filter_all <- function(y,
                       wavelet_wf = "la8",
                       wavelet_J = 4,
                       wavelet_levels = c(3, 4),
                       hpgc_prior = "weak",
                       hpgc_chains = 4,
                       hpgc_iterations = 20000,
                       hpgc_burnin = 5000,
                       verbose = FALSE) {
  
  y <- as.numeric(y)
  n <- length(y)
  
  if (verbose) message("Applying wavelet filter...")
  wavelet_result <- filter_wavelet(
    y,
    wf = wavelet_wf,
    J = wavelet_J,
    levels_to_combine = wavelet_levels,
    first_difference = TRUE,
    verbose = verbose
  )
  
  if (verbose) message("Applying EMD filter...")
  emd_result <- filter_emd(y, verbose = verbose)
  
  if (verbose) message("Applying HP-GC filter...")
  hpgc_result <- filter_hpgc(
    y,
    prior_config = hpgc_prior,
    n_chains = hpgc_chains,
    iterations = hpgc_iterations,
    burnin = hpgc_burnin,
    verbose = verbose
  )
  
  comparison <- data.frame(
    original = y,
    wavelet_trend = wavelet_result$trend,
    emd_trend = emd_result$trend,
    hpgc_trend = hpgc_result$trend,
    hpgc_lower = hpgc_result$trend_lower,
    hpgc_upper = hpgc_result$trend_upper
  )
  
  correlations <- stats::cor(comparison[, c("wavelet_trend", "emd_trend", "hpgc_trend")])
  
  result <- list(
    wavelet = wavelet_result,
    emd = emd_result,
    hpgc = hpgc_result,
    comparison = comparison,
    correlations = correlations,
    n_obs = n
  )
  
  class(result) <- c("signaly_multifilter", "list")
  result
}