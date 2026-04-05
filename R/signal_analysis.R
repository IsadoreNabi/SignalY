#' @title Comprehensive Signal Analysis for Panel Data
#' @name signal_analysis
#'
#' @description
#' Master function that orchestrates the complete signal extraction pipeline,
#' integrating spectral decomposition (wavelets, EMD, HP-GC), Bayesian variable' selection (regularized Horseshoe), dimensionality reduction (PCA, DFM), and
#' stationarity testing into a unified analytical framework.
#'
#' The function constructs a target signal Y from candidate variables X in panel
#' data and applies multiple complementary methodologies to extract the latent
#' structure from phenomenological dynamics.
#'
#' @param data A data.frame or matrix containing the panel data. For data.frames,
#'   time should be in rows and variables in columns.
#' @param y_formula Formula specifying how to construct Y from X variables, or
#'   a character string naming the pre-constructed Y column in data.
#' @param time_var Character string naming the time variable (optional, assumes
#'   rows are ordered by time if NULL).
#' @param group_var Character string naming the group/panel variable for panel
#'   data (optional for single time series).
#' @param methods Character vector specifying which methods to apply. Options:
#'   \code{"wavelet"}, \code{"emd"}, \code{"hpgc"}, \code{"horseshoe"},
#'   \code{"pca"}, \code{"dfm"}, \code{"unitroot"}, or \code{"all"} (default).
#' @param filter_config List of configuration options for filtering methods:
#'   \describe{
#'     \item{wavelet_filter}{Wavelet filter type (default: "la8")}
#'     \item{wavelet_levels}{Which detail levels to combine (default: c(3,4))}
#'     \item{emd_max_imf}{Maximum IMFs for EMD (default: 10)}
#'     \item{hpgc_prior}{Prior configuration: "weak", "informative", "empirical" (default: "weak")}
#'     \item{hpgc_chains}{Number of MCMC chains (default: 4)}
#'     \item{hpgc_iterations}{Total iterations per chain (default: 20000)}
#'   }
#' @param horseshoe_config List of configuration for Horseshoe regression:
#'   \describe{
#'     \item{p0}{Expected number of relevant predictors (default: NULL for auto)}
#'     \item{chains}{Number of MCMC chains (default: 4)}
#'     \item{iter_sampling}{Sampling iterations per chain (default: 2000)}
#'     \item{iter_warmup}{Warmup iterations (default: 1000)}
#'     \item{adapt_delta}{Target acceptance rate (default: 0.95)}
#'     \item{use_qr}{Use QR decomposition (default: TRUE)}
#'     \item{selection_method}{Variable selection method: \code{"shrinkage"}
#'       (default, kappa-based), \code{"credible_interval"} (CI excludes zero),
#'       \code{"magnitude"} (|beta| threshold), or \code{"projection"}
#'       (predictive projection — most robust)}
#'     \item{kappa_threshold}{Shrinkage threshold for \code{"shrinkage"} method
#'       (default: 0.5)}
#'     \item{ci_prob}{Probability level for \code{"credible_interval"} method
#'       (default: 0.95)}
#'     \item{magnitude_threshold}{Minimum |beta| for \code{"magnitude"} method
#'       (default: 0.1)}
#'     \item{projection_cv_k}{Number of CV folds for \code{"projection"} method
#'       (default: 5)}
#'   }
#' @param pca_config List of configuration for PCA:
#'   \describe{
#'     \item{n_components}{Number of components (default: NULL for auto)}
#'     \item{rotation}{Rotation method: "none", "varimax", "oblimin" (default: "none")}
#'     \item{n_boot}{Bootstrap replications (default: 1000)}
#'     \item{block_length}{Block length for bootstrap (default: NULL for auto)}
#'     \item{alpha}{Alpha for bootstrap tests (default: 0.05)}
#'   }
#' @param dfm_config List of configuration for Dynamic Factor Models:
#'   \describe{
#'     \item{r}{Number of factors (default: NULL for auto via IC)}
#'     \item{max_factors}{Maximum factors to consider (default: 10)}
#'     \item{p}{VAR lags for factor dynamics (default: 1)}
#'     \item{ic}{Information criterion: "IC1", "IC2", "IC3" (default: "bai_ng_2")}
#'   }
#' @param unitroot_tests Character vector of unit root tests to apply. Options:
#'   \code{"adf"}, \code{"ers"}, \code{"kpss"}, \code{"pp"}, or \code{"all"} (default).
#' @param na_action How to handle missing values: "interpolate", "omit", "fail" (default: "interpolate").
#' @param standardize Logical, whether to standardize variables before analysis (default: TRUE).
#' @param first_difference Logical, whether to first-difference data (default: FALSE).
#' @param verbose Logical, whether to print progress messages (default: TRUE).
#' @param seed Random seed for reproducibility (default: NULL).
#'
#' @return An S3 object of class \code{"signal_analysis"} containing:
#'   \describe{
#'     \item{call}{The matched function call}
#'     \item{data}{Processed input data}
#'     \item{Y}{The constructed target signal}
#'     \item{X}{The predictor matrix}
#'     \item{filters}{Results from spectral decomposition methods}
#'     \item{horseshoe}{Results from Bayesian variable selection}
#'     \item{pca}{Results from PCA with bootstrap}
#'     \item{dfm}{Results from Dynamic Factor Model}
#'     \item{unitroot}{Results from unit root tests}
#'     \item{interpretation}{Automated technical interpretation}
#'     \item{config}{Configuration parameters used}
#'   }
#'
#' @details
#' \strong{Methodological Framework}
#'
#' The signal extraction pipeline distinguishes between latent structure
#' (the underlying data-generating process) and phenomenological dynamics
#' (observed variability). This is achieved through:
#'
#' \enumerate{
#'   \item \strong{Spectral Decomposition}: Separates signal frequencies
#'     \itemize{
#'       \item Wavelets: Multi-resolution analysis via MODWT
#'       \item EMD: Data-adaptive decomposition into intrinsic modes
#'       \item HP-GC: Bayesian unobserved components (trend + cycle)
#'     }
#'   \item \strong{Sparse Regression}: Identifies relevant predictors
#'     \itemize{
#'       \item Regularized Horseshoe: Adaptive shrinkage with slab regularization
#'       \item Shrinkage factors (kappa) quantify predictor relevance
#'     }
#'   \item \strong{Dimensionality Reduction}: Extracts common factors
#'     \itemize{
#'       \item PCA: Static factor structure with bootstrap significance
#'       \item DFM: Dynamic factors with VAR transition dynamics
#'     }
#'   \item \strong{Stationarity Testing}: Characterizes persistence properties
#'     \itemize{
#'       \item Integrated battery of ADF, ERS, KPSS, PP tests
#'       \item Synthesized conclusion on stationarity type
#'     }
#' }
#'
#' \strong{Interpretation Framework}
#'
#' The automated interpretation assesses:
#' \itemize{
#'   \item \strong{Signal Smoothness}: Variance of second differences
#'   \item \strong{Trend Persistence}: Deterministic vs. stochastic via unit roots
#'   \item \strong{Information Topology}: Entropy of PC1 loadings (concentrated vs. diffuse)
#'   \item \strong{Sparsity Ratio}: Proportion of predictors shrunk to zero
#'   \item \strong{Factor Structure}: Number of significant common factors
#' }
#'
#' @examples
#' \donttest{
#' # Generate example panel data
#' set.seed(42)
#' n_time <- 50   
#' n_vars <- 10   
#' 
#' # Create correlated predictors with common factor structure
#' factors <- matrix(rnorm(n_time * 2), n_time, 2)
#' loadings <- matrix(runif(n_vars * 2, -1, 1), n_vars, 2)
#' X <- factors %*% t(loadings) + matrix(rnorm(n_time * n_vars, 0, 0.5), n_time, n_vars)
#' colnames(X) <- paste0("X", 1:n_vars)
#' 
#' # True signal depends on only 3 predictors
#' true_beta <- c(rep(1, 3), rep(0, 7))
#' Y <- X %*% true_beta + rnorm(n_time, 0, 0.5)
#' 
#' # Combine into data frame
#' data <- data.frame(Y = Y, X)
#' 
#' # Run comprehensive analysis
#' # We pass specific configs to make MCMC very fast just for the example
#' result <- signal_analysis(
#'   data = data,
#'   y_formula = "Y",
#'   methods = "all",
#'   verbose = TRUE,
#'   # Configuration for speed (CRAN policy < 5s preferred)
#'   filter_config = list(
#'      hpgc_chains = 1,      
#'      hpgc_iterations = 50, 
#'      hpgc_burnin = 10
#'   ),
#'   horseshoe_config = list(
#'      chains = 1,           
#'      iter_sampling = 50,   
#'      iter_warmup = 10
#'   ),
#'   pca_config = list(
#'      n_boot = 50           
#'   )
#' )
#' 
#' # View interpretation
#' print(result)
#' 
#' # Plot results (requires plotly)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   plot(result)
#' }
#' }
#'
#' @seealso
#' \code{\link{filter_wavelet}}, \code{\link{filter_emd}}, \code{\link{filter_hpgc}},
#' \code{\link{fit_horseshoe}}, \code{\link{select_by_projection}},
#' \code{\link{pca_bootstrap}}, \code{\link{estimate_dfm}},
#' \code{\link{test_unit_root}}
#'
#' @references
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of Statistics,
#' 11(2), 5018-5051. \doi{10.1214/17-EJS1337SI}
#'
#' Bai, J., & Ng, S. (2002). Determining the Number of Factors in Approximate
#' Factor Models. Econometrica, 70(1), 191-221. \doi{10.1111/1468-0262.00273}
#'
#' @export
signal_analysis <- function(data,
                            y_formula,
                            time_var = NULL,
                            group_var = NULL,
                            methods = "all",
                            filter_config = list(),
                            horseshoe_config = list(),
                            pca_config = list(),
                            dfm_config = list(),
                            unitroot_tests = "all",
                            na_action = c("interpolate", "omit", "fail"),
                            standardize = TRUE,
                            first_difference = FALSE,
                            verbose = TRUE,
                            seed = NULL) {
  
  call <- match.call()
  
  if (!is.null(seed)) set.seed(seed)
  
  na_action <- match.arg(na_action)
  
  all_methods <- c("wavelet", "emd", "hpgc", "horseshoe", "pca", "dfm", "unitroot")
  if ("all" %in% methods) {
    methods <- all_methods
  }
  methods <- match.arg(methods, all_methods, several.ok = TRUE)
  
  if (verbose) message("=== SignalY: Comprehensive Signal Analysis ===\n")
  if (verbose) message("[1/7] Preparing data...")
  
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  
  if (is.character(y_formula) && length(y_formula) == 1) {
    if (!y_formula %in% names(data)) {
      stop("Variable '", y_formula, "' not found in data.")
    }
    Y <- data[[y_formula]]
    X_names <- setdiff(names(data), c(y_formula, time_var, group_var))
    X <- as.matrix(data[, X_names, drop = FALSE])
  } else if (inherits(y_formula, "formula")) {
    mf <- model.frame(y_formula, data = data)
    Y <- model.response(mf)
    X <- model.matrix(y_formula, data = data)[, -1, drop = FALSE]  
    X_names <- colnames(X)
  } else {
    stop("y_formula must be a character string or formula.")
  }
  
  n <- length(Y)
  p <- ncol(X)
  
  if (verbose) message("  - Observations: ", n, ", Predictors: ", p)
  
  na_count <- sum(is.na(Y)) + sum(is.na(X))
  if (na_count > 0) {
    if (verbose) message("  - Missing values detected: ", na_count)
    if (na_action == "fail") {
      stop("Missing values present and na_action = 'fail'.")
    } else if (na_action == "omit") {
      complete_cases <- complete.cases(Y, X)
      Y <- Y[complete_cases]
      X <- X[complete_cases, , drop = FALSE]
      n <- length(Y)
      if (verbose) message("  - Removed ", sum(!complete_cases), " incomplete cases.")
    } else {
      Y <- interpolate_na(Y)
      X <- apply(X, 2, interpolate_na)
      if (verbose) message("  - Interpolated missing values.")
    }
  }
  
  if (first_difference) {
    Y <- diff(Y)
    X <- apply(X, 2, diff)
    n <- length(Y)
    if (verbose) message("  - Applied first differencing. New n = ", n)
  }
  
  if (standardize) {
    Y_mean <- mean(Y)
    Y_sd <- sd(Y)
    Y_std <- (Y - Y_mean) / Y_sd
    
    X_means <- colMeans(X)
    X_sds <- apply(X, 2, sd)
    X_std <- scale(X)
    
    if (verbose) message("  - Standardized variables.")
  } else {
    Y_std <- Y
    X_std <- X
    Y_mean <- 0
    Y_sd <- 1
    X_means <- rep(0, p)
    X_sds <- rep(1, p)
  }
  
  results <- list(
    call = call,
    data = list(
      original = data,
      Y = Y,
      X = X,
      Y_std = Y_std,
      X_std = X_std,
      n = n,
      p = p,
      standardization = list(
        Y_mean = Y_mean, Y_sd = Y_sd,
        X_means = X_means, X_sds = X_sds
      )
    ),
    filters = NULL,
    horseshoe = NULL,
    pca = NULL,
    dfm = NULL,
    unitroot = NULL,
    interpretation = NULL,
    config = list(
      methods = methods,
      filter_config = filter_config,
      horseshoe_config = horseshoe_config,
      pca_config = pca_config,
      dfm_config = dfm_config,
      unitroot_verbose = verbose,
      na_action = na_action,
      standardize = standardize,
      first_difference = first_difference,
      seed = seed
    )
  )
  
  filter_methods <- intersect(methods, c("wavelet", "emd", "hpgc"))
  
  if (length(filter_methods) > 0) {
    if (verbose) message("\n[2/7] Applying spectral decomposition...")
    
    results$filters <- list()
    
    fc <- list(
      wavelet_filter = "la8",
      wavelet_levels = c(3, 4),
      emd_max_imf = 10,
      hpgc_prior = "weak",
      hpgc_chains = 4,
      hpgc_iterations = 20000,
      hpgc_burnin = 5000
    )
    fc <- modifyList(fc, filter_config)
    
    if ("wavelet" %in% filter_methods) {
      if (verbose) message("  - Wavelet MRA (", fc$wavelet_filter, ")...")
      results$filters$wavelet <- tryCatch({
        res <- filter_wavelet(
          y = Y_std,
          wf = fc$wavelet_filter,
          levels_to_combine = fc$wavelet_levels,
          first_difference = FALSE
        )
        if (verbose) {
          message("    Done. Combined levels: D",
                  paste(fc$wavelet_levels, collapse = "+D"))
        }
        res
      }, error = function(e) {
        warning("Wavelet filtering failed: ", e$message)
        NULL
      })
    }

    if ("emd" %in% filter_methods) {
      if (verbose) message("  - Empirical Mode Decomposition...")
      results$filters$emd <- tryCatch({
        res <- filter_emd(
          y = Y_std,
          max_imf = fc$emd_max_imf
        )
        if (verbose) {
          message("    Done. Extracted ",
                  res$n_imf, " IMFs.")
        }
        res
      }, error = function(e) {
        warning("EMD failed: ", e$message)
        NULL
      })
    }

    if ("hpgc" %in% filter_methods) {
      if (verbose) message("  - Bayesian HP-GC filter...")
      results$filters$hpgc <- tryCatch({
        res <- filter_hpgc(
          y = Y_std,
          prior_config = fc$hpgc_prior,
          n_chains = fc$hpgc_chains,
          iterations = fc$hpgc_iterations,
          burnin = fc$hpgc_burnin
        )
        if (verbose) {
          message("    Done. Selected prior: ",
                  res$selected_prior)
        }
        res
      }, error = function(e) {
        warning("HP-GC filtering failed: ", e$message)
        NULL
      })
    }
  }
  
  if ("horseshoe" %in% methods) {
    if (verbose) message("\n[3/7] Fitting regularized Horseshoe regression...")
    
    hc <- list(
      p0 = NULL,
      chains = 4,
      iter_sampling = 2000,
      iter_warmup = 1000,
      adapt_delta = 0.95,
      use_qr = TRUE,
      selection_method = "shrinkage",
      kappa_threshold = 0.5,
      ci_prob = 0.95,
      magnitude_threshold = 0.1,
      projection_cv_k = 5
    )
    hc <- modifyList(hc, horseshoe_config)
    
    results$horseshoe <- tryCatch({
      hs_fit <- fit_horseshoe(
        y = Y_std,
        X = X_std,
        p0 = hc$p0,
        chains = hc$chains,
        iter_sampling = hc$iter_sampling,
        iter_warmup = hc$iter_warmup,
        adapt_delta = hc$adapt_delta,
        use_qr = hc$use_qr
      )

      hs_fit$selection <- switch(
        hc$selection_method,
        "shrinkage" = select_by_shrinkage(
          hs_fit, threshold = hc$kappa_threshold
        ),
        "credible_interval" = select_by_credible_interval(
          hs_fit, prob = hc$ci_prob
        ),
        "magnitude" = select_by_magnitude(
          hs_fit, threshold = hc$magnitude_threshold
        ),
        "projection" = select_by_projection(
          hs_fit,
          cv_k = hc$projection_cv_k,
          verbose = verbose
        ),
        stop(sprintf(
          "Unknown selection_method '%s'.",
          hc$selection_method
        ), call. = FALSE)
      )

      if (verbose) {
        n_sel <- length(hs_fit$selection$selected)
        message("  - Selected ", n_sel, " of ", p,
                " predictors (method: ",
                hc$selection_method, ")")
      }
      hs_fit
    }, error = function(e) {
      warning("Horseshoe regression failed: ", e$message)
      NULL
    })
  }
  
  if ("pca" %in% methods) {
    if (verbose) message("\n[4/7] Principal Component Analysis with bootstrap...")
    
    pc <- list(
      n_components = NULL,
      rotation = "none",
      n_boot = 1000,
      block_length = NULL,
      alpha = 0.05
    )
    pc <- modifyList(pc, pca_config)
    
    results$pca <- tryCatch({
      res <- pca_bootstrap(
        X = X_std,
        n_components = pc$n_components,
        rotation = pc$rotation,
        n_boot = pc$n_boot,
        block_length = pc$block_length,
        alpha = pc$alpha
      )
      if (verbose) {
        n_sig <- sum(res$significant_loadings)
        message("  - ", res$n_components,
                " components retained, ",
                n_sig, " significant loadings.")
        message("  - PC1 entropy: ",
                round(res$entropy[1], 3))
      }
      res
    }, error = function(e) {
      warning("PCA failed: ", e$message)
      NULL
    })
  }
  
  if ("dfm" %in% methods) {
    if (verbose) message("\n[5/7] Estimating Dynamic Factor Model...")
    
    dc <- list(
      r = NULL,
      max_factors = 10,
      p = 1,
      ic = "IC2"
    )
    dc <- modifyList(dc, dfm_config)
    
    results$dfm <- tryCatch({
      res <- estimate_dfm(
        X = X_std,
        r = dc$r,
        max_factors = dc$max_factors,
        p = dc$p,
        ic = dc$ic
      )
      if (verbose) {
        message("  - Optimal factors: ",
                res$r_selected, " (", dc$ic, ")")
        message("  - Variance explained: ",
                round(sum(res$variance_explained) * 100,
                      1), "%")
      }
      res
    }, error = function(e) {
      warning("DFM estimation failed: ", e$message)
      NULL
    })
  }
  
  if ("unitroot" %in% methods) {
    if (verbose) message("\n[6/7] Testing for unit roots...")
    
    results$unitroot <- tryCatch({
      res <- test_unit_root(y = Y_std, verbose = verbose)
      if (verbose) {
        message("  - Conclusion: ", res$conclusion)
      }
      res
    }, error = function(e) {
      warning("Unit root testing failed: ", e$message)
      NULL
    })
  }
  
  if (verbose) message("\n[7/7] Generating technical interpretation...")
  
  results$interpretation <- generate_interpretation(results, verbose = verbose)
  
  class(results) <- "signal_analysis"
  
  if (verbose) {
    message("\n=== Analysis Complete ===")
    message("Use print(), summary(), or plot() to examine results.")
  }
  
  return(results)
}


#' @title Generate Automated Technical Interpretation
#' @description Internal function to synthesize results into interpretable narrative.
#' @param results signal_analysis results object
#' @param verbose Logical, print progress
#' @return List containing interpretation components
#' @keywords internal
generate_interpretation <- function(results, verbose = FALSE) {
  
  interpretation <- list(
    signal_characteristics = list(),
    variable_selection = list(),
    factor_structure = list(),
    persistence = list(),
    overall_summary = NULL
  )
  
  if (!is.null(results$filters)) {
    Y <- results$data$Y_std
    n <- length(Y)
    
    d2Y <- diff(diff(Y))
    smoothness <- var(d2Y)
    interpretation$signal_characteristics$smoothness <- smoothness
    interpretation$signal_characteristics$smoothness_interpretation <- 
      if (smoothness < 0.01) "Very smooth (strong trend dominance)"
    else if (smoothness < 0.1) "Moderately smooth"
    else if (smoothness < 0.5) "Moderately volatile"
    else "Highly volatile (noise-dominated)"
    
    if (!is.null(results$filters$hpgc)) {
      trend_var <- var(results$filters$hpgc$trend)
      cycle_var <- var(results$filters$hpgc$cycle)
      trend_share <- trend_var / (trend_var + cycle_var)
      interpretation$signal_characteristics$trend_share <- trend_share
      interpretation$signal_characteristics$trend_interpretation <-
        if (trend_share > 0.7) "Trend-dominated dynamics"
      else if (trend_share > 0.3) "Mixed trend-cycle dynamics"
      else "Cycle-dominated dynamics"
    }
    
    if (!is.null(results$filters$wavelet)) {
      interpretation$signal_characteristics$wavelet_energy <- 
        results$filters$wavelet$level_variance
    }
  }
  
  if (!is.null(results$horseshoe)) {
    hs <- results$horseshoe
    
    kappa_mean <- hs$coefficients$kappa_mean
    n_selected <- length(hs$selection$selected)
    sparsity_ratio <- 1 - n_selected / results$data$p
    
    interpretation$variable_selection$sparsity_ratio <- sparsity_ratio
    interpretation$variable_selection$n_selected <- n_selected
    interpretation$variable_selection$m_eff <- hs$sparsity$m_eff_mean
    interpretation$variable_selection$selected_vars <- hs$selection$selected
    
    interpretation$variable_selection$interpretation <-
      if (sparsity_ratio > 0.9) "Extremely sparse: very few predictors relevant"
    else if (sparsity_ratio > 0.7) "Highly sparse: clear signal concentration"
    else if (sparsity_ratio > 0.5) "Moderately sparse"
    else if (sparsity_ratio > 0.3) "Moderately dense"
    else "Dense: many predictors contribute"
    
    beta_abs <- abs(hs$coefficients$beta_mean)
    top_idx <- order(beta_abs, decreasing = TRUE)[1:min(5, length(beta_abs))]
    interpretation$variable_selection$top_predictors <- data.frame(
      variable = hs$coefficients$variable[top_idx],
      beta = hs$coefficients$beta_mean[top_idx],
      kappa = kappa_mean[top_idx]
    )
  }
  
  if (!is.null(results$pca)) {
    pca <- results$pca
    
    entropy_pc1 <- pca$entropy[1]
    max_entropy <- log(results$data$p)
    normalized_entropy <- entropy_pc1 / max_entropy
    
    interpretation$factor_structure$pc1_entropy <- entropy_pc1
    interpretation$factor_structure$pc1_entropy_normalized <- normalized_entropy
    interpretation$factor_structure$n_significant_components <- 
      sum(pca$variance_explained > 0.05)
    interpretation$factor_structure$total_variance_explained <- 
      sum(pca$variance_explained[1:pca$n_components])
    
    interpretation$factor_structure$topology_interpretation <-
      if (normalized_entropy < 0.3) "Concentrated: few variables dominate PC1"
    else if (normalized_entropy < 0.6) "Moderately distributed loadings"
    else "Diffuse: signal spread across many variables"
    
    loadings_pc1 <- abs(pca$loadings[, 1])
    top_idx <- order(loadings_pc1, decreasing = TRUE)[1:min(5, length(loadings_pc1))]
    interpretation$factor_structure$top_pc1_loadings <- data.frame(
      variable = rownames(pca$loadings)[top_idx],
      loading = pca$loadings[top_idx, 1],
      abs_loading = loadings_pc1[top_idx]
    )
  }
  
  if (!is.null(results$dfm)) {
    interpretation$factor_structure$dfm_n_factors <- results$dfm$n_factors
    interpretation$factor_structure$dfm_variance_explained <- 
      sum(results$dfm$variance_explained)
  }
  
  if (!is.null(results$unitroot)) {
    ur <- results$unitroot
    
    interpretation$persistence$conclusion <- ur$conclusion
    interpretation$persistence$persistence_type <- ur$persistence_type
    
    ur_conclusion <- if (!is.null(ur$conclusion)) {
      if (grepl("STATIONARY.*stationary", ur$conclusion, ignore.case = TRUE)) {
        "stationary"
      } else if (grepl("UNIT ROOT", ur$conclusion, ignore.case = TRUE)) {
        "unit_root"
      } else if (grepl("trend.stationary", ur$conclusion, ignore.case = TRUE)) {
        "trend_stationary"
      } else {
        "inconclusive"
      }
    } else {
      "inconclusive"
    }
    
    interpretation$persistence$interpretation <- switch(
      ur_conclusion,
      "stationary" = "Signal is mean-reverting with transient shocks",
      "trend_stationary" = "Deterministic trend with stationary deviations",
      "difference_stationary" = "Stochastic trend (unit root); shocks have permanent effects",
      "unit_root" = "Stochastic trend (unit root); shocks have permanent effects",
      "inconclusive" = "Mixed evidence; persistence properties unclear",
      "Mixed evidence; persistence properties unclear"
    )
  }
  
  summary_parts <- c()
  
  if (!is.null(interpretation$signal_characteristics$smoothness_interpretation)) {
    summary_parts <- c(summary_parts, 
                       interpretation$signal_characteristics$smoothness_interpretation)
  }
  
  if (!is.null(interpretation$signal_characteristics$trend_interpretation)) {
    summary_parts <- c(summary_parts,
                       interpretation$signal_characteristics$trend_interpretation)
  }
  
  if (!is.null(interpretation$variable_selection$interpretation)) {
    summary_parts <- c(summary_parts,
                       interpretation$variable_selection$interpretation)
  }
  
  if (!is.null(interpretation$factor_structure$topology_interpretation)) {
    summary_parts <- c(summary_parts,
                       interpretation$factor_structure$topology_interpretation)
  }
  
  if (!is.null(interpretation$persistence$interpretation)) {
    summary_parts <- c(summary_parts,
                       interpretation$persistence$interpretation)
  }
  
  interpretation$overall_summary <- paste(summary_parts, collapse = ". ")
  if (nchar(interpretation$overall_summary) > 0) {
    interpretation$overall_summary <- paste0(interpretation$overall_summary, ".")
  }
  
  return(interpretation)
}

#' @title Print Method for signal_analysis Objects
#' @description Print a concise summary of signal analysis results.
#' @param x An object of class \code{signal_analysis}
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @export
print.signal_analysis <- function(x, ...) {
  cat("\n")
  cat("===========================================================================\n")
  cat("                    SignalY: Signal Analysis Results                        \n")
  cat("===========================================================================\n\n")
  
  cat("DATA SUMMARY\n")
  cat("-----------------------------------------------------------------------------\n")
  cat("  Observations:", x$data$n, "\n")
  cat("  Predictors:  ", x$data$p, "\n")
  cat("  Methods:     ", paste(x$config$methods, collapse = ", "), "\n")
  cat("\n")
  
  if (!is.null(x$interpretation$overall_summary) && 
      nchar(x$interpretation$overall_summary) > 0) {
    cat("SYNTHESIS\n")
    cat("-----------------------------------------------------------------------------\n")
    summary_wrapped <- strwrap(x$interpretation$overall_summary, width = 75)
    cat(" ", paste(summary_wrapped, collapse = "\n  "), "\n\n")
  }
  
  cat("KEY METRICS\n")
  cat("-----------------------------------------------------------------------------\n")
  
  if (!is.null(x$horseshoe)) {
    cat("  Horseshoe Regression:\n")
    cat("    - Selected predictors:", length(x$horseshoe$selection$selected), 
        "of", x$data$p, "\n")
    cat("    - Effective non-zeros (m_eff):", 
        round(x$horseshoe$sparsity$m_eff_mean, 2), "\n")
    cat("    - Sparsity ratio:", 
        round(x$horseshoe$sparsity$sparsity_ratio * 100, 1), "%\n")
  }
  
  if (!is.null(x$pca)) {
    cat("  PCA:\n")
    cat("    - Components retained:", x$pca$n_components, "\n")
    cat("    - Variance explained:", 
        round(sum(x$pca$variance_explained[1:x$pca$n_components]) * 100, 1), "%\n")
    cat("    - PC1 entropy:", round(x$pca$entropy[1], 3), "\n")
  }
  
  if (!is.null(x$dfm)) {
    cat("  Dynamic Factor Model:\n")
    cat("    - Optimal factors:", x$dfm$n_factors, "\n")
    cat("    - Variance explained:", 
        round(sum(x$dfm$variance_explained) * 100, 1), "%\n")
  }
  
  if (!is.null(x$unitroot)) {
    cat("  Unit Root Tests:\n")
    cat("    - Conclusion:", x$unitroot$conclusion, "\n")
    cat("    - Type:", x$unitroot$persistence_type, "\n")
  }
  
  cat("\n")
  cat("===========================================================================\n")
  cat("Use summary() for detailed results, plot() for visualizations.\n")
  
  invisible(x)
}


#' @title Summary Method for signal_analysis Objects
#' @description Generate a detailed summary of signal analysis results.
#' @param object An object of class \code{signal_analysis}
#' @param ... Additional arguments (ignored)
#' @return A list containing detailed summaries (invisibly)
#' @export
summary.signal_analysis <- function(object, ...) {
  
  cat("\n")
  cat("+===========================================================================+\n")
  cat("|               SignalY: Detailed Signal Analysis Summary                   |\n")
  cat("+===========================================================================+\n\n")
  
  cat("+-----------------------------------------------------------------------------+\n")
  cat("| 1. DATA OVERVIEW                                                           |\n")
  cat("+-----------------------------------------------------------------------------+\n\n")
  
  cat("  Sample size (n):", object$data$n, "\n")
  cat("  Number of predictors (p):", object$data$p, "\n")
  cat("  Standardized:", object$config$standardize, "\n")
  cat("  First-differenced:", object$config$first_difference, "\n")
  cat("  NA handling:", object$config$na_action, "\n\n")
  
  if (!is.null(object$filters)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 2. SPECTRAL DECOMPOSITION                                                  |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    if (!is.null(object$filters$wavelet)) {
      cat("  WAVELET MULTI-RESOLUTION ANALYSIS\n")
      cat("  -------------------------------------\n")
      wf_name <- object$filters$wavelet$filter %||% object$config$filter_config$wavelet_filter %||% "N/A"
      cat("  Filter:", wf_name, "\n")
      cat("  Levels combined:", paste(object$filters$wavelet$levels, collapse = ", "), "\n")
      comb_var <- tryCatch(var(object$filters$wavelet$combined), error = function(e) NA)
      cat("  Combined signal variance:", 
          if (!is.na(comb_var)) round(comb_var, 4) else "N/A", "\n\n")
    }
    
    if (!is.null(object$filters$emd)) {
      cat("  EMPIRICAL MODE DECOMPOSITION\n")
      cat("  -------------------------------------\n")
      cat("  Number of IMFs:", if (!is.null(object$filters$emd$n_imf)) object$filters$emd$n_imf else "N/A", "\n")
      res_var <- tryCatch(var(object$filters$emd$residue), error = function(e) NA)
      cat("  Residue variance:", if (!is.na(res_var)) round(res_var, 4) else "N/A", "\n\n")
    }
    
    if (!is.null(object$filters$hpgc)) {
      cat("  BAYESIAN HP-GC FILTER\n")
      cat("  -------------------------------------\n")
      cat("  Selected prior:", object$filters$hpgc$selected_prior %||% "N/A", "\n")
      trend_var <- tryCatch(var(object$filters$hpgc$trend), error = function(e) NA)
      cycle_var <- tryCatch(var(object$filters$hpgc$cycle), error = function(e) NA)
      cat("  Trend variance:", if (!is.na(trend_var)) round(trend_var, 4) else "N/A", "\n")
      cat("  Cycle variance:", if (!is.na(cycle_var)) round(cycle_var, 4) else "N/A", "\n")
      if (!is.null(object$filters$hpgc$dic) && is.numeric(object$filters$hpgc$dic)) {
        cat("  DIC:", round(object$filters$hpgc$dic, 2), "\n")
      }
      cat("\n")
    }
    
    if (!is.null(object$interpretation$signal_characteristics)) {
      sc <- object$interpretation$signal_characteristics
      cat("  SIGNAL CHARACTERISTICS\n")
      cat("  -------------------------------------\n")
      if (!is.null(sc$smoothness) && is.numeric(sc$smoothness)) {
        cat("  Smoothness (var of Delta^2Y):", round(sc$smoothness, 4), "\n")
        cat("  Interpretation:", sc$smoothness_interpretation %||% "N/A", "\n")
      }
      if (!is.null(sc$trend_share) && is.numeric(sc$trend_share)) {
        cat("  Trend share:", round(sc$trend_share * 100, 1), "%\n")
        cat("  Interpretation:", sc$trend_interpretation %||% "N/A", "\n")
      }
      cat("\n")
    }
  }
  
  if (!is.null(object$horseshoe)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 3. REGULARIZED HORSESHOE REGRESSION                                        |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    hs <- object$horseshoe
    
    cat("  MCMC DIAGNOSTICS\n")
    cat("  -------------------------------------\n")
    if (!is.null(hs$diagnostics)) {
      n_div <- hs$diagnostics$n_divergences
      rhat <- hs$diagnostics$rhat_max
      ess_bulk <- hs$diagnostics$ess_bulk_min
      ess_tail <- hs$diagnostics$ess_tail_min
      
      cat("  Divergences:", if (!is.null(n_div)) n_div else "N/A", "\n")
      cat("  Max R-hat:", if (!is.null(rhat) && is.numeric(rhat)) round(rhat, 3) else "N/A", "\n")
      cat("  Min bulk ESS:", if (!is.null(ess_bulk) && is.numeric(ess_bulk)) round(ess_bulk, 0) else "N/A", "\n")
      cat("  Min tail ESS:", if (!is.null(ess_tail) && is.numeric(ess_tail)) round(ess_tail, 0) else "N/A", "\n")
    }
    cat("\n")
    
    cat("  SHRINKAGE SUMMARY\n")
    cat("  -------------------------------------\n")
    tau_val <- NA
    if (!is.null(hs$hyperparameters) && is.data.frame(hs$hyperparameters)) {
      tau_row <- hs$hyperparameters[hs$hyperparameters$parameter == "tau", ]
      if (nrow(tau_row) > 0 && is.numeric(tau_row$mean)) {
        tau_val <- tau_row$mean[1]
      }
    }
    if (!is.na(tau_val)) {
      cat("  Global tau (mean):", round(tau_val, 4), "\n")
    } else {
      cat("  Global tau (mean): N/A\n")
    }
    
    m_eff_val <- if (!is.null(hs$sparsity$m_eff_mean) && is.numeric(hs$sparsity$m_eff_mean)) {
      hs$sparsity$m_eff_mean
    } else NA
    if (!is.na(m_eff_val)) {
      cat("  Effective non-zeros (m_eff):", round(m_eff_val, 2), "\n")
    } else {
      cat("  Effective non-zeros (m_eff): N/A\n")
    }
    
    cat("  Selected predictors:", 
        if (!is.null(hs$selection$selected)) length(hs$selection$selected) else 0, 
        "of", object$data$p, "\n")
    cat("  Kappa threshold:", hs$selection$threshold %||% 0.5, "\n\n")
    
    cat("  TOP PREDICTORS (by |beta|)\n")
    cat("  -------------------------------------\n")
    if (!is.null(object$interpretation$variable_selection$top_predictors)) {
      tp <- object$interpretation$variable_selection$top_predictors
      if (is.data.frame(tp) && nrow(tp) > 0) {
        threshold <- hs$selection$threshold %||% 0.5
        for (i in 1:nrow(tp)) {
          var_name <- if (!is.null(tp$variable[i])) tp$variable[i] else paste0("V", i)
          beta_val <- if (!is.null(tp$beta[i]) && is.numeric(tp$beta[i])) tp$beta[i] else NA
          kappa_val <- if (!is.null(tp$kappa[i]) && is.numeric(tp$kappa[i])) tp$kappa[i] else NA
          
          if (!is.na(beta_val) && !is.na(kappa_val)) {
            cat(sprintf("  %2d. %-20s beta = %7.4f  kappa = %.3f %s\n",
                        i, var_name, beta_val, kappa_val,
                        ifelse(kappa_val < threshold, "*", "")))
          }
        }
      }
    }
    cat("\n")
    
    if (!is.null(hs$loo) && !is.null(hs$loo$estimates)) {
      cat("  MODEL FIT (LOO-CV)\n")
      cat("  -------------------------------------\n")
      
      elpd_est <- tryCatch(hs$loo$estimates["elpd_loo", "Estimate"], error = function(e) NA)
      elpd_se <- tryCatch(hs$loo$estimates["elpd_loo", "SE"], error = function(e) NA)
      looic <- tryCatch(hs$loo$estimates["looic", "Estimate"], error = function(e) NA)
      p_loo <- tryCatch(hs$loo$estimates["p_loo", "Estimate"], error = function(e) NA)
      
      if (!is.na(elpd_est) && !is.na(elpd_se)) {
        cat("  ELPD:", round(elpd_est, 2), "+/-", round(elpd_se, 2), "\n")
      }
      if (!is.na(looic)) {
        cat("  LOOIC:", round(looic, 2), "\n")
      }
      if (!is.na(p_loo)) {
        cat("  p_loo:", round(p_loo, 2), "\n\n")
      }
    }
  }
  
  if (!is.null(object$pca)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 4. PRINCIPAL COMPONENT ANALYSIS                                            |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    pca <- object$pca
    
    cat("  VARIANCE EXPLAINED\n")
    cat("  -------------------------------------\n")
    n_comp <- if (!is.null(pca$n_components) && is.numeric(pca$n_components)) pca$n_components else 0
    if (n_comp > 0 && !is.null(pca$variance_explained) && length(pca$variance_explained) >= n_comp) {
      for (i in 1:min(5, n_comp)) {
        ve_i <- pca$variance_explained[i]
        ve_cum <- sum(pca$variance_explained[1:i])
        if (is.numeric(ve_i) && is.numeric(ve_cum)) {
          cat(sprintf("  PC%d: %.1f%% (cumulative: %.1f%%)\n", i, ve_i * 100, ve_cum * 100))
        }
      }
    }
    cat("\n")
    
    cat("  ENTROPY ANALYSIS\n")
    cat("  -------------------------------------\n")
    if (!is.null(pca$entropy) && length(pca$entropy) > 0 && is.numeric(pca$entropy[1])) {
      cat("  PC1 entropy:", round(pca$entropy[1], 3), "\n")
    }
    if (!is.null(object$interpretation$factor_structure$pc1_entropy_normalized) &&
        is.numeric(object$interpretation$factor_structure$pc1_entropy_normalized)) {
      cat("  Normalized entropy:", 
          round(object$interpretation$factor_structure$pc1_entropy_normalized, 3), "\n")
      cat("  Interpretation:", 
          object$interpretation$factor_structure$topology_interpretation, "\n")
    }
    cat("\n")
    
    cat("  TOP PC1 LOADINGS\n")
    cat("  -------------------------------------\n")
    if (!is.null(object$interpretation$factor_structure$top_pc1_loadings)) {
      tl <- object$interpretation$factor_structure$top_pc1_loadings
      if (is.data.frame(tl) && nrow(tl) > 0) {
        for (i in 1:nrow(tl)) {
          cat(sprintf("  %2d. %-20s loading = %7.4f\n",
                      i, tl$variable[i], tl$loading[i]))
        }
      }
    }
    cat("\n")
    
    cat("  BOOTSTRAP SIGNIFICANCE\n")
    cat("  -------------------------------------\n")
    cat("  Bootstrap replications:", if (!is.null(pca$n_boot)) pca$n_boot else "N/A", "\n")
    cat("  Significance level:", if (!is.null(pca$significance_level)) pca$significance_level else "N/A", "\n")
    if (!is.null(pca$loadings_significant) && is.logical(pca$loadings_significant)) {
      cat("  Significant loadings:", sum(pca$loadings_significant), 
          "of", length(pca$loadings_significant), "\n\n")
    } else {
      cat("  Significant loadings: N/A\n\n")
    }
  }
  
  if (!is.null(object$dfm)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 5. DYNAMIC FACTOR MODEL                                                    |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    dfm <- object$dfm
    
    cat("  FACTOR SELECTION\n")
    cat("  -------------------------------------\n")
    cat("  Information criterion:", object$config$dfm_config$ic %||% "bai_ng_2", "\n")
    cat("  Optimal factors:", if (!is.null(dfm$n_factors)) dfm$n_factors else "N/A", "\n")
    
    var_expl <- if (!is.null(dfm$variance_explained) && is.numeric(dfm$variance_explained) && length(dfm$variance_explained) > 0) {
      round(sum(dfm$variance_explained) * 100, 1)
    } else "N/A"
    cat("  Total variance explained:", var_expl, if (is.numeric(var_expl)) "%" else "", "\n\n")
    
    if (!is.null(dfm$ic_values) && length(dfm$ic_values) > 0) {
      cat("  INFORMATION CRITERIA\n")
      cat("  -------------------------------------\n")
      for (i in 1:min(5, length(dfm$ic_values))) {
        if (is.numeric(dfm$ic_values[i]) && !is.na(dfm$ic_values[i])) {
          cat(sprintf("  %d factors: IC = %.4f\n", i, dfm$ic_values[i]))
        }
      }
      cat("\n")
    }
  }
  
  if (!is.null(object$unitroot)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 6. UNIT ROOT TESTS                                                         |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    ur <- object$unitroot
    
    cat("  SYNTHESIS\n")
    cat("  -------------------------------------\n")
    cat("  Conclusion:", ur$conclusion %||% "N/A", "\n")
    cat("  Type:", ur$persistence_type %||% "N/A", "\n")
    cat("  Interpretation:", object$interpretation$persistence$interpretation %||% "N/A", "\n\n")
    
    cat("  TEST RESULTS SUMMARY\n")
    cat("  -------------------------------------\n")
    
    if (!is.null(ur$tests) && length(ur$tests) > 0) {
      for (test_name in names(ur$tests)) {
        test <- ur$tests[[test_name]]
        stat_val <- test$statistic
        p_val <- test$p_value
        if (!is.null(stat_val) && !is.null(p_val) && 
            is.numeric(stat_val) && is.numeric(p_val)) {
          cat(sprintf("  %-20s stat = %7.3f  p = %.4f\n",
                      test_name, stat_val, p_val))
        }
      }
    }
    cat("\n")
  }
  
  cat("+-----------------------------------------------------------------------------+\n")
  cat("| 7. OVERALL SYNTHESIS                                                        |\n")
  cat("+-----------------------------------------------------------------------------+\n\n")
  
  if (!is.null(object$interpretation$overall_summary) &&
      nchar(object$interpretation$overall_summary) > 0) {
    summary_wrapped <- strwrap(object$interpretation$overall_summary, width = 75)
    cat(" ", paste(summary_wrapped, collapse = "\n  "), "\n\n")
  }
  
  class(object) <- c("summary.signal_analysis", class(object))
  
  invisible(object)
}

#' @title Plot Method for signal_analysis Objects
#' @description
#' Generate interactive diagnostic visualizations for signal
#' analysis results using plotly. Produces a unified HTML widget
#' dashboard covering filters, horseshoe regression, PCA, and
#' DFM results.
#'
#' @param x An object of class \code{signal_analysis}.
#' @param which Character vector specifying which panels to
#'   include. Options: \code{"all"} (default), \code{"filters"},
#'   \code{"horseshoe"}, \code{"pca"}, \code{"dfm"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A plotly htmlwidget (invisibly). Prints/displays
#'   automatically in interactive sessions.
#' @export
plot.signal_analysis <- function(x, which = "all", ...) {

  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for plot(). ",
      "Install with: install.packages('plotly')",
      call. = FALSE
    )
  }

  all_panels <- c("filters", "horseshoe", "pca", "dfm")
  if ("all" %in% which) which <- all_panels
  which <- intersect(which, all_panels)

  panels <- list()

  # ----------------------------------------------------------
  # FILTERS panel
  # ----------------------------------------------------------
  if ("filters" %in% which && !is.null(x$filters)) {
    Y <- x$data$Y_std
    n <- length(Y)
    ti <- seq_len(n)

    p_base <- plotly::plot_ly(
      x = ti, y = Y,
      name = "Original", type = "scatter",
      mode = "lines",
      line = list(color = "black", width = 1,
                  opacity = 0.3)
    )

    if (!is.null(x$filters$wavelet)) {
      pw <- plotly::add_trace(
        p_base, y = x$filters$wavelet$combined,
        name = "Wavelet Trend",
        line = list(color = "#00cc96", width = 2)
      )
      pw <- plotly::layout(
        pw, yaxis = list(title = "Wavelet")
      )
      panels[["wavelet"]] <- pw
    }

    if (!is.null(x$filters$hpgc)) {
      ph <- plotly::add_trace(
        p_base, y = x$filters$hpgc$trend,
        name = "HP-GC Trend",
        line = list(color = "#636efa", width = 2)
      )
      ph <- plotly::layout(
        ph, yaxis = list(title = "HP-GC Trend")
      )
      panels[["hpgc_trend"]] <- ph

      pc <- plotly::plot_ly(
        x = ti, y = x$filters$hpgc$cycle,
        name = "HP-GC Cycle", type = "scatter",
        mode = "lines",
        line = list(color = "#ef553b", width = 1.5)
      )
      pc <- plotly::layout(
        pc, yaxis = list(title = "HP-GC Cycle")
      )
      pc <- plotly::add_lines(
        pc, x = c(1, n), y = c(0, 0),
        line = list(color = "gray", dash = "dash",
                    width = 1),
        showlegend = FALSE, inherit = FALSE
      )
      panels[["hpgc_cycle"]] <- pc
    }

    if (!is.null(x$filters$emd)) {
      pe <- plotly::add_trace(
        p_base, y = x$filters$emd$trend,
        name = "EMD Trend",
        line = list(color = "#ab63fa", width = 2)
      )
      pe <- plotly::layout(
        pe, yaxis = list(title = "EMD")
      )
      panels[["emd"]] <- pe
    }
  }

  # ----------------------------------------------------------
  # HORSESHOE panel
  # ----------------------------------------------------------
  if ("horseshoe" %in% which && !is.null(x$horseshoe)) {
    hs <- x$horseshoe
    coef_df <- hs$coefficients
    sel_vars <- hs$selection$selected %||% character(0)
    threshold <- hs$selection$threshold %||% 0.5

    ord <- order(abs(coef_df$beta_mean),
                 decreasing = TRUE)
    top_n <- min(20L, nrow(coef_df))
    ord <- ord[seq_len(top_n)]
    df <- coef_df[ord, ]
    df$is_selected <- df$variable %in% sel_vars
    df$color <- ifelse(df$is_selected,
                       "#636efa", "#aaaaaa")
    # order top to bottom
    df$variable <- factor(
      df$variable,
      levels = rev(df$variable)
    )

    # Coefficient forest plot
    p_coef <- plotly::plot_ly(
      data = df, y = ~variable, x = ~beta_mean,
      type = "scatter", mode = "markers",
      marker = list(color = ~color, size = 8),
      error_x = list(
        type = "data",
        symmetric = FALSE,
        arrayminus = ~(beta_mean - beta_q025),
        array = ~(beta_q975 - beta_mean),
        color = ~color
      ),
      name = "beta (95% CI)",
      showlegend = FALSE
    )
    p_coef <- plotly::layout(
      p_coef,
      shapes = list(list(
        type = "line", x0 = 0, x1 = 0,
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "red", dash = "dash")
      )),
      xaxis = list(title = "beta"),
      yaxis = list(title = ""),
      title = "Coefficients (95% CI)"
    )
    panels[["hs_coef"]] <- p_coef

    # Kappa bar plot
    p_kap <- plotly::plot_ly(
      data = df, y = ~variable, x = ~kappa_mean,
      type = "bar", orientation = "h",
      marker = list(
        color = ifelse(df$kappa_mean < threshold,
                       "#636efa", "#aaaaaa")
      ),
      name = "kappa", showlegend = FALSE
    )
    p_kap <- plotly::layout(
      p_kap,
      shapes = list(list(
        type = "line",
        x0 = threshold, x1 = threshold,
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "red", dash = "dash",
                    width = 2)
      )),
      xaxis = list(title = "kappa"),
      yaxis = list(title = ""),
      title = "Shrinkage (kappa)"
    )
    panels[["hs_kappa"]] <- p_kap

    # kappa vs |beta| scatter
    p_kb <- plotly::plot_ly(
      data = coef_df,
      x = ~kappa_mean, y = ~abs(beta_mean),
      type = "scatter", mode = "markers",
      text = ~variable, hoverinfo = "text+x+y",
      marker = list(
        color = ifelse(
          coef_df$variable %in% sel_vars,
          "#636efa", "#aaaaaa"
        ),
        size = 8
      ),
      showlegend = FALSE
    )
    p_kb <- plotly::layout(
      p_kb,
      shapes = list(list(
        type = "line",
        x0 = threshold, x1 = threshold,
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "red", dash = "dash")
      )),
      xaxis = list(title = "kappa"),
      yaxis = list(title = "|beta|"),
      title = "Shrinkage vs |beta|"
    )
    panels[["hs_scatter"]] <- p_kb

    # PPC if available
    if (!is.null(hs$ppc) && !is.null(hs$ppc$y_rep)) {
      y_obs <- x$data$Y_std
      y_rep_mean <- colMeans(hs$ppc$y_rep)
      rng <- range(c(y_obs, y_rep_mean), na.rm = TRUE)
      p_ppc <- plotly::plot_ly(
        x = y_obs, y = y_rep_mean,
        type = "scatter", mode = "markers",
        marker = list(color = "#636efa", size = 5,
                      opacity = 0.6),
        name = "Obs vs Pred", showlegend = FALSE
      )
      p_ppc <- plotly::add_lines(
        p_ppc, x = rng, y = rng,
        line = list(color = "red", width = 2,
                    dash = "dash"),
        showlegend = FALSE, inherit = FALSE
      )
      p_ppc <- plotly::layout(
        p_ppc,
        xaxis = list(title = "Observed Y"),
        yaxis = list(title = "Predicted Y (mean)"),
        title = "Posterior Predictive Check"
      )
      panels[["hs_ppc"]] <- p_ppc
    }
  }

  # ----------------------------------------------------------
  # PCA panel
  # ----------------------------------------------------------
  if ("pca" %in% which && !is.null(x$pca)) {
    pca <- x$pca
    p_dim <- x$data$p
    vn <- colnames(x$data$X)
    if (is.null(vn)) vn <- paste0("X", seq_len(p_dim))

    n_show <- min(10L, length(pca$variance_explained))
    ve <- pca$variance_explained[seq_len(n_show)]
    pc_names <- paste0("PC", seq_len(n_show))
    col_scree <- ifelse(
      seq_len(n_show) <= pca$n_components,
      "#636efa", "#aaaaaa"
    )

    # Scree plot
    p_scree <- plotly::plot_ly(
      x = pc_names, y = ve * 100, type = "bar",
      marker = list(color = col_scree),
      name = "Var %", showlegend = FALSE
    )
    p_scree <- plotly::add_trace(
      p_scree,
      x = pc_names, y = cumsum(ve) * 100,
      type = "scatter", mode = "lines+markers",
      line = list(color = "red"),
      marker = list(color = "red", size = 6),
      name = "Cumulative", showlegend = FALSE,
      yaxis = "y"
    )
    p_scree <- plotly::layout(
      p_scree,
      xaxis = list(title = "Principal Component",
                    categoryorder = "array",
                    categoryarray = pc_names),
      yaxis = list(title = "Variance Explained (%)"),
      title = "Scree Plot"
    )
    panels[["pca_scree"]] <- p_scree

    # PC1 loadings
    l1 <- pca$loadings[, 1]
    l_ord <- order(abs(l1), decreasing = TRUE)
    l_top <- l_ord[seq_len(min(20L, length(l1)))]
    load_df <- data.frame(
      var = factor(vn[l_top], levels = rev(vn[l_top])),
      loading = l1[l_top],
      stringsAsFactors = FALSE
    )
    p_load <- plotly::plot_ly(
      data = load_df, y = ~var, x = ~loading,
      type = "bar", orientation = "h",
      marker = list(
        color = ifelse(load_df$loading > 0,
                       "#636efa", "#ef553b")
      ),
      showlegend = FALSE
    )
    p_load <- plotly::layout(
      p_load,
      shapes = list(list(
        type = "line", x0 = 0, x1 = 0,
        y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "gray", dash = "dash")
      )),
      xaxis = list(title = "Loading"),
      yaxis = list(title = ""),
      title = "PC1 Loadings (Top 20)"
    )
    panels[["pca_loadings"]] <- p_load

    # Score plot PC1 vs PC2
    if (pca$n_components >= 2) {
      p_score <- plotly::plot_ly(
        x = pca$scores[, 1],
        y = pca$scores[, 2],
        type = "scatter", mode = "markers",
        marker = list(color = "#636efa", size = 5,
                      opacity = 0.6),
        showlegend = FALSE
      )
      p_score <- plotly::layout(
        p_score,
        xaxis = list(title = sprintf(
          "PC1 (%.1f%%)",
          pca$variance_explained[1] * 100
        )),
        yaxis = list(title = sprintf(
          "PC2 (%.1f%%)",
          pca$variance_explained[2] * 100
        )),
        title = "Score Plot (PC1 vs PC2)"
      )
      panels[["pca_scores"]] <- p_score
    }

    # Entropy
    n_ent <- min(5L, length(pca$entropy))
    ent_names <- paste0("PC", seq_len(n_ent))
    p_ent <- plotly::plot_ly(
      x = ent_names,
      y = pca$entropy[seq_len(n_ent)],
      type = "bar",
      marker = list(color = "#00cc96"),
      showlegend = FALSE
    )
    p_ent <- plotly::layout(
      p_ent,
      shapes = list(list(
        type = "line", x0 = 0, x1 = 1,
        xref = "paper",
        y0 = log(p_dim), y1 = log(p_dim),
        line = list(color = "red", dash = "dash")
      )),
      xaxis = list(title = "Principal Component",
                    categoryorder = "array",
                    categoryarray = ent_names),
      yaxis = list(title = "Shannon Entropy"),
      title = "Loading Entropy by Component"
    )
    panels[["pca_entropy"]] <- p_ent
  }

  # ----------------------------------------------------------
  # DFM panel
  # ----------------------------------------------------------
  if ("dfm" %in% which && !is.null(x$dfm)) {
    dfm <- x$dfm
    n <- x$data$n
    ti <- seq_len(n)

    # IC plot
    if (!is.null(dfm$ic_values)) {
      n_ic <- length(dfm$ic_values)
      ic_x <- seq_len(n_ic)
      p_ic <- plotly::plot_ly(
        x = ic_x, y = dfm$ic_values,
        type = "scatter", mode = "lines+markers",
        marker = list(
          color = ifelse(
            ic_x == dfm$r_selected, "red",
            "#636efa"
          ),
          size = ifelse(
            ic_x == dfm$r_selected, 12, 6
          )
        ),
        line = list(color = "#636efa"),
        showlegend = FALSE
      )
      p_ic <- plotly::layout(
        p_ic,
        xaxis = list(title = "Number of Factors"),
        yaxis = list(title = "IC Value"),
        title = "Information Criterion"
      )
      panels[["dfm_ic"]] <- p_ic
    }

    # Factor time series
    if (!is.null(dfm$factors)) {
      n_plot <- min(3L, ncol(dfm$factors))
      for (fi in seq_len(n_plot)) {
        pf <- plotly::plot_ly(
          x = ti, y = dfm$factors[, fi],
          type = "scatter", mode = "lines",
          line = list(color = "#636efa",
                      width = 1.5),
          showlegend = FALSE
        )
        pf <- plotly::add_lines(
          pf, x = c(1, n), y = c(0, 0),
          line = list(color = "gray", dash = "dash",
                      width = 1),
          showlegend = FALSE, inherit = FALSE
        )
        pf <- plotly::layout(
          pf,
          xaxis = list(title = "Time"),
          yaxis = list(title = paste0("F", fi)),
          title = paste("Factor", fi)
        )
        panels[[paste0("dfm_f", fi)]] <- pf
      }
    }

    # Variance explained
    if (!is.null(dfm$variance_explained)) {
      nv <- length(dfm$variance_explained)
      f_names <- paste0("F", seq_len(nv))
      p_ve <- plotly::plot_ly(
        x = f_names,
        y = dfm$variance_explained * 100,
        type = "bar",
        marker = list(color = "#636efa"),
        showlegend = FALSE
      )
      p_ve <- plotly::layout(
        p_ve,
        xaxis = list(title = "Factor",
                      categoryorder = "array",
                      categoryarray = f_names),
        yaxis = list(title = "Var Explained (%)"),
        title = "Variance Explained by Factor"
      )
      panels[["dfm_ve"]] <- p_ve
    }
  }

  # ----------------------------------------------------------
  # ASSEMBLE
  # ----------------------------------------------------------
  if (length(panels) == 0L) {
    message("No results available to plot.")
    return(invisible(NULL))
  }

  n_panels <- length(panels)
  final <- plotly::subplot(
    panels,
    nrows = ceiling(n_panels / 2),
    shareX = FALSE,
    shareY = FALSE,
    titleX = TRUE,
    titleY = TRUE,
    margin = 0.06
  )

  final <- plotly::layout(
    final,
    title = list(
      text = "SignalY: Interactive Analysis Dashboard",
      x = 0.5
    ),
    showlegend = FALSE
  )

  final
}


#' @title Interactive Plot for Signal Analysis (Alias)
#' @description
#' Alias for \code{plot(x)} kept for backwards compatibility.
#' @param x An object of class \code{signal_analysis}.
#' @param ... Passed to \code{plot.signal_analysis}.
#' @return A plotly htmlwidget.
#' @export
iplot <- function(x, ...) {
  plot.signal_analysis(x, ...)
}