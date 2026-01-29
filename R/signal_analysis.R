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
#'     \item{kappa_threshold}{Shrinkage threshold for selection (default: 0.5)}
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
#' # Plot results
#' plot(result)
#' }
#'
#' @seealso
#' \code{\link{filter_wavelet}}, \code{\link{filter_emd}}, \code{\link{filter_hpgc}},
#' \code{\link{fit_horseshoe}}, \code{\link{pca_bootstrap}}, \code{\link{estimate_dfm}},
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
  
  # ===========================================================================
  # Data Preparation
  # ===========================================================================
  
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
  
  # ===========================================================================
  # Initialize Results
  # ===========================================================================
  
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
  
  # ===========================================================================
  # Spectral Decomposition (Filters)
  # ===========================================================================
  
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
      tryCatch({
        results$filters$wavelet <- filter_wavelet(
          y = Y_std,
          wf = fc$wavelet_filter,
          levels_to_combine = fc$wavelet_levels,
          first_difference = FALSE
        )
        if (verbose) message("    Done. Combined levels: D", 
                             paste(fc$wavelet_levels, collapse = "+D"))
      }, error = function(e) {
        warning("Wavelet filtering failed: ", e$message)
        results$filters$wavelet <<- NULL
      })
    }
    
    if ("emd" %in% filter_methods) {
      if (verbose) message("  - Empirical Mode Decomposition...")
      tryCatch({
        results$filters$emd <- filter_emd(
          y = Y_std,
          max_imf = fc$emd_max_imf
        )
        if (verbose) message("    Done. Extracted ", results$filters$emd$n_imf, " IMFs.")
      }, error = function(e) {
        warning("EMD failed: ", e$message)
        results$filters$emd <<- NULL
      })
    }
    
    if ("hpgc" %in% filter_methods) {
      if (verbose) message("  - Bayesian HP-GC filter...")
      tryCatch({
        results$filters$hpgc <- filter_hpgc(
          y = Y_std,
          prior_config = fc$hpgc_prior,
          n_chains = fc$hpgc_chains,
          iterations = fc$hpgc_iterations,
          burnin = fc$hpgc_burnin
        )
        if (verbose) message("    Done. Selected prior: ", results$filters$hpgc$selected_prior)
      }, error = function(e) {
        warning("HP-GC filtering failed: ", e$message)
        results$filters$hpgc <<- NULL
      })
    }
  }
  
  # ===========================================================================
  # Bayesian Variable Selection (Horseshoe)
  # ===========================================================================
  
  if ("horseshoe" %in% methods) {
    if (verbose) message("\n[3/7] Fitting regularized Horseshoe regression...")
    
    hc <- list(
      p0 = NULL,
      chains = 4,
      iter_sampling = 2000,
      iter_warmup = 1000,
      adapt_delta = 0.95,
      use_qr = TRUE,
      kappa_threshold = 0.5
    )
    hc <- modifyList(hc, horseshoe_config)
    
    tryCatch({
      results$horseshoe <- fit_horseshoe(
        y = Y_std,
        X = X_std,
        p0 = hc$p0,
        chains = hc$chains,
        iter_sampling = hc$iter_sampling,
        iter_warmup = hc$iter_warmup,
        adapt_delta = hc$adapt_delta,
        use_qr = hc$use_qr
      )
      
      results$horseshoe$selection <- select_by_shrinkage(
        results$horseshoe,
        threshold = hc$kappa_threshold
      )
      
      if (verbose) {
        n_selected <- length(results$horseshoe$selection$selected)
        message("  - Selected ", n_selected, " of ", p, " predictors (kappa < ", 
                hc$kappa_threshold, ")")
      }
    }, error = function(e) {
      warning("Horseshoe regression failed: ", e$message)
      results$horseshoe <<- NULL
    })
  }
  
  # ===========================================================================
  # PCA with Bootstrap
  # ===========================================================================
  
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
    
    tryCatch({
      results$pca <- pca_bootstrap(
        X = X_std,
        n_components = pc$n_components,
        rotation = pc$rotation,
        n_boot = pc$n_boot,
        block_length = pc$block_length,
        alpha = pc$alpha
      )
      
      if (verbose) {
        n_sig <- sum(results$pca$loadings_significant)
        message("  - ", results$pca$n_components, " components retained, ",
                n_sig, " significant loadings.")
        message("  - PC1 entropy: ", round(results$pca$entropy[1], 3))
      }
    }, error = function(e) {
      warning("PCA failed: ", e$message)
      results$pca <<- NULL
    })
  }
  
  # ===========================================================================
  # Dynamic Factor Model
  # ===========================================================================
  
  if ("dfm" %in% methods) {
    if (verbose) message("\n[5/7] Estimating Dynamic Factor Model...")
    
    dc <- list(
      r = NULL,
      max_factors = 10,
      p = 1,
      ic = "IC2"
    )
    dc <- modifyList(dc, dfm_config)
    
    tryCatch({
      results$dfm <- estimate_dfm(
        X = X_std,
        r = dc$r,
        max_factors = dc$max_factors,
        p = dc$p,
        ic = dc$ic
      )
      
      if (verbose) {
        message("  - Optimal factors: ", results$dfm$n_factors, 
                " (", dc$ic, ")")
        message("  - Variance explained: ", 
                round(sum(results$dfm$variance_explained) * 100, 1), "%")
      }
    }, error = function(e) {
      warning("DFM estimation failed: ", e$message)
      results$dfm <<- NULL
    })
  }
  
  # ===========================================================================
  # Unit Root Tests
  # ===========================================================================
  
  if ("unitroot" %in% methods) {
    if (verbose) message("\n[6/7] Testing for unit roots...")
    
    tryCatch({
      results$unitroot <- test_unit_root(
        y = Y_std,
        verbose = verbose
      )
      
      if (verbose) {
        message("  - Synthesis: ", results$unitroot$synthesis$conclusion)
      }
    }, error = function(e) {
      warning("Unit root testing failed: ", e$message)
      results$unitroot <<- NULL
    })
  }
  
  # ===========================================================================
  # Automated Technical Interpretation
  # ===========================================================================
  
  if (verbose) message("\n[7/7] Generating technical interpretation...")
  
  results$interpretation <- generate_interpretation(results, verbose = verbose)
  
  # ===========================================================================
  # Finalize
  # ===========================================================================
  
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
  
  # ---------------------------------------------------------------------------
  # Signal Characteristics (from filters)
  # ---------------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------------
  # Variable Selection (from Horseshoe)
  # ---------------------------------------------------------------------------
  
  if (!is.null(results$horseshoe)) {
    hs <- results$horseshoe
    
    # Usar la estructura correcta: $coefficients y $sparsity
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
    
    # Usar $coefficients$mean en lugar de $summary$beta_mean
    beta_abs <- abs(hs$coefficients$mean)
    top_idx <- order(beta_abs, decreasing = TRUE)[1:min(5, length(beta_abs))]
    interpretation$variable_selection$top_predictors <- data.frame(
      variable = hs$coefficients$variable[top_idx],
      beta = hs$coefficients$mean[top_idx],
      kappa = kappa_mean[top_idx]
    )
  }
  
  # ---------------------------------------------------------------------------
  # Factor Structure (from PCA/DFM)
  # ---------------------------------------------------------------------------
  
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
  
  # ---------------------------------------------------------------------------
  # Persistence Properties (from Unit Root)
  # ---------------------------------------------------------------------------
  
  if (!is.null(results$unitroot)) {
    ur <- results$unitroot
    
    interpretation$persistence$conclusion <- ur$synthesis$conclusion
    interpretation$persistence$confidence <- ur$synthesis$confidence
    interpretation$persistence$evidence <- ur$synthesis$evidence
    
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
  
  # ---------------------------------------------------------------------------
  # Overall Summary
  # ---------------------------------------------------------------------------
  
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


# =============================================================================
# S3 Methods: print, summary, plot
# =============================================================================

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
  
  # Data summary
  cat("DATA SUMMARY\n")
  cat("-----------------------------------------------------------------------------\n")
  cat("  Observations:", x$data$n, "\n")
  cat("  Predictors:  ", x$data$p, "\n")
  cat("  Methods:     ", paste(x$config$methods, collapse = ", "), "\n")
  cat("\n")
  
  # Interpretation
  if (!is.null(x$interpretation$overall_summary) && 
      nchar(x$interpretation$overall_summary) > 0) {
    cat("SYNTHESIS\n")
    cat("-----------------------------------------------------------------------------\n")
    # Word wrap the summary
    summary_wrapped <- strwrap(x$interpretation$overall_summary, width = 75)
    cat(" ", paste(summary_wrapped, collapse = "\n  "), "\n\n")
  }
  
  # Key metrics
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
    cat("    - Conclusion:", x$unitroot$synthesis$conclusion, "\n")
    cat("    - Confidence:", x$unitroot$synthesis$confidence, "\n")
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
  
  # ===========================================================================
  # 1. Data Overview
  # ===========================================================================
  
  cat("+-----------------------------------------------------------------------------+\n")
  cat("| 1. DATA OVERVIEW                                                           |\n")
  cat("+-----------------------------------------------------------------------------+\n\n")
  
  cat("  Sample size (n):", object$data$n, "\n")
  cat("  Number of predictors (p):", object$data$p, "\n")
  cat("  Standardized:", object$config$standardize, "\n")
  cat("  First-differenced:", object$config$first_difference, "\n")
  cat("  NA handling:", object$config$na_action, "\n\n")
  
  # ===========================================================================
  # 2. Spectral Decomposition
  # ===========================================================================
  
  if (!is.null(object$filters)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 2. SPECTRAL DECOMPOSITION                                                  |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    if (!is.null(object$filters$wavelet)) {
      cat("  WAVELET MULTI-RESOLUTION ANALYSIS\n")
      cat("  -------------------------------------\n")
      # CORREGIDO: Fallback robusto para el nombre del filtro
      wf_name <- object$filters$wavelet$filter %||% object$config$filter_config$wavelet_filter %||% "N/A"
      cat("  Filter:", wf_name, "\n")
      cat("  Levels combined:", paste(object$filters$wavelet$levels, collapse = ", "), "\n")
      # Robusto para combined signal variance
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
    
    # Signal characteristics
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
  
  # ===========================================================================
  # 3. Horseshoe Regression
  # ===========================================================================
  
  if (!is.null(object$horseshoe)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 3. REGULARIZED HORSESHOE REGRESSION                                        |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    hs <- object$horseshoe
    
    cat("  MCMC DIAGNOSTICS\n")
    cat("  -------------------------------------\n")
    if (!is.null(hs$diagnostics)) {
      # Acceso robusto a diagnósticos
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
    # Usar la estructura correcta: $hyperparameters es un data frame
    # Acceso robusto a tau
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
    
    # Acceso robusto a m_eff
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
      
      # Acceso robusto a las estimaciones de LOO
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
  
  # ===========================================================================
  # 4. PCA
  # ===========================================================================
  
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
  
  # ===========================================================================
  # 5. Dynamic Factor Model
  # ===========================================================================
  
  if (!is.null(object$dfm)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 5. DYNAMIC FACTOR MODEL                                                    |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    dfm <- object$dfm
    
    cat("  FACTOR SELECTION\n")
    cat("  -------------------------------------\n")
    # CORREGIDO: usar ic en lugar de ic_criterion
    cat("  Information criterion:", object$config$dfm_config$ic %||% "bai_ng_2", "\n")
    cat("  Optimal factors:", if (!is.null(dfm$n_factors)) dfm$n_factors else "N/A", "\n")
    
    # Acceso robusto a variance_explained
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
  
  # ===========================================================================
  # 6. Unit Root Tests
  # ===========================================================================
  
  if (!is.null(object$unitroot)) {
    cat("+-----------------------------------------------------------------------------+\n")
    cat("| 6. UNIT ROOT TESTS                                                         |\n")
    cat("+-----------------------------------------------------------------------------+\n\n")
    
    ur <- object$unitroot
    
    cat("  SYNTHESIS\n")
    cat("  -------------------------------------\n")
    cat("  Conclusion:", ur$synthesis$conclusion %||% "N/A", "\n")
    cat("  Confidence:", ur$synthesis$confidence %||% "N/A", "\n")
    cat("  Interpretation:", object$interpretation$persistence$interpretation %||% "N/A", "\n\n")
    
    cat("  TEST RESULTS SUMMARY\n")
    cat("  -------------------------------------\n")
    
    # Print summary table of tests
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
  
  # ===========================================================================
  # 7. Overall Summary
  # ===========================================================================
  
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
#' @description Generate diagnostic visualizations for signal analysis results.
#' @param x An object of class \code{signal_analysis}
#' @param which Character vector specifying which plots to create. Options:
#'   \code{"all"}, \code{"filters"}, \code{"horseshoe"}, \code{"pca"}, \code{"dfm"},
#'   \code{"unitroot"}. Default is \code{"all"}.
#' @param ask Logical, whether to prompt before each plot (default: TRUE in interactive mode)
#' @param ... Additional arguments passed to plotting functions
#' @return Invisibly returns the input object
#' @export
plot.signal_analysis <- function(x, which = "all", ask = NULL, ...) {
  
  # Determine if interactive
  if (is.null(ask)) {
    ask <- interactive() && length(which) > 1
  }
  
  # Expand "all"
  all_plots <- c("filters", "horseshoe", "pca", "dfm")
  if ("all" %in% which) {
    which <- all_plots
  }
  
  # Set up plotting
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  # ===========================================================================
  # Filter Plots
  # ===========================================================================
  
  if ("filters" %in% which && !is.null(x$filters)) {
    
    # Time series and decomposition
    graphics::par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))
    
    Y <- x$data$Y_std
    n <- length(Y)
    time_idx <- 1:n
    
    # Original series
    graphics::plot(time_idx, Y, type = "l", col = "black", lwd = 1.5,
                   main = "Original Series (Standardized)",
                   xlab = "Time", ylab = "Y")
    graphics::grid()
    
    # HP-GC decomposition if available
    if (!is.null(x$filters$hpgc)) {
      graphics::plot(time_idx, Y, type = "l", col = "gray60", lwd = 1,
                     main = "HP-GC Trend-Cycle Decomposition",
                     xlab = "Time", ylab = "Y")
      graphics::lines(time_idx, x$filters$hpgc$trend, col = "blue", lwd = 2)
      graphics::legend("topright", legend = c("Original", "Trend"),
                       col = c("gray60", "blue"), lwd = c(1, 2), bty = "n")
      graphics::grid()
      
      graphics::plot(time_idx, x$filters$hpgc$cycle, type = "l", col = "red", lwd = 1.5,
                     main = "Cyclical Component",
                     xlab = "Time", ylab = "Cycle")
      graphics::abline(h = 0, lty = 2, col = "gray40")
      graphics::grid()
    } else if (!is.null(x$filters$wavelet)) {
      # Wavelet combined signal
      graphics::plot(time_idx, x$filters$wavelet$combined, type = "l", 
                     col = "darkgreen", lwd = 1.5,
                     main = paste("Wavelet Combined Signal (D", 
                                  paste(x$filters$wavelet$levels, collapse = "+D"), ")", sep = ""),
                     xlab = "Time", ylab = "Combined")
      graphics::grid()
      
      # Smooth component
      graphics::plot(time_idx, x$filters$wavelet$smooth, type = "l",
                     col = "purple", lwd = 1.5,
                     main = "Wavelet Smooth Component",
                     xlab = "Time", ylab = "Smooth")
      graphics::grid()
    }
    
    # EMD if available
    if (!is.null(x$filters$emd) && x$filters$emd$n_imf >= 2) {
      graphics::par(mfrow = c(min(4, x$filters$emd$n_imf + 1), 1), 
                    mar = c(3, 4, 2, 2))
      
      graphics::plot(time_idx, Y, type = "l", col = "black", lwd = 1,
                     main = "EMD Decomposition", xlab = "", ylab = "Original")
      
      for (i in 1:min(3, x$filters$emd$n_imf)) {
        graphics::plot(time_idx, x$filters$emd$imf[, i], type = "l",
                       col = "steelblue", lwd = 1,
                       main = "", xlab = "", ylab = paste("IMF", i))
      }
    }
  }
  
  # ===========================================================================
  # Horseshoe Plots
  # ===========================================================================
  
  if ("horseshoe" %in% which && !is.null(x$horseshoe)) {
    
    hs <- x$horseshoe
    p <- x$data$p
    var_names <- colnames(x$data$X)
    if (is.null(var_names)) var_names <- paste0("X", 1:p)
    
    graphics::par(mfrow = c(2, 2), mar = c(5, 6, 4, 2))
    
    # Usar la estructura correcta: $coefficients es un data frame
    # Nota: coefficients está ordenado por relevance_score, necesitamos reordenar
    coef_df <- hs$coefficients
    
    # Crear vectores con el orden original de variables
    beta_mean <- coef_df$mean
    beta_lower <- coef_df$q025
    beta_upper <- coef_df$q975
    kappa_mean <- coef_df$kappa_mean
    coef_var_names <- coef_df$variable
    
    # Determinar cuáles están seleccionados (es un vector de nombres)
    is_selected <- coef_var_names %in% hs$selection$selected
    
    # Ordenar por |beta| para el plot
    ord <- order(abs(beta_mean), decreasing = TRUE)
    top_n <- min(20, length(beta_mean))
    ord <- ord[1:top_n]
    
    y_pos <- top_n:1
    xlim <- range(c(beta_lower[ord], beta_upper[ord]), na.rm = TRUE)
    
    graphics::plot(beta_mean[ord], y_pos, xlim = xlim, yaxt = "n",
                   pch = 19, col = ifelse(is_selected[ord], "blue", "gray60"),
                   main = "Coefficient Estimates (95% CI)",
                   xlab = expression(beta), ylab = "")
    graphics::segments(beta_lower[ord], y_pos, beta_upper[ord], y_pos,
                       col = ifelse(is_selected[ord], "blue", "gray60"))
    graphics::abline(v = 0, lty = 2, col = "red")
    graphics::axis(2, at = y_pos, labels = coef_var_names[ord], las = 1, cex.axis = 0.7)
    
    # Shrinkage factors (kappa)
    graphics::barplot(kappa_mean[ord], horiz = TRUE, names.arg = coef_var_names[ord],
                      col = ifelse(kappa_mean[ord] < hs$selection$threshold, 
                                   "steelblue", "gray70"),
                      main = "Shrinkage Factors (kappa)",
                      xlab = expression(kappa), las = 1, cex.names = 0.7)
    graphics::abline(v = hs$selection$threshold, lty = 2, col = "red", lwd = 2)
    
    # Kappa vs |beta|
    graphics::plot(kappa_mean, abs(beta_mean), pch = 19,
                   col = ifelse(is_selected, "blue", "gray60"),
                   main = expression("Shrinkage vs |"*beta*"|"),
                   xlab = expression(kappa), ylab = expression("|"*beta*"|"))
    graphics::abline(v = hs$selection$threshold, lty = 2, col = "red")
    graphics::grid()
    
    # Posterior predictive check
    if (!is.null(hs$ppc) && !is.null(hs$ppc$y_rep_mean)) {
      y_rep_mean <- hs$ppc$y_rep_mean
      y_obs <- x$data$Y_std
      
      graphics::plot(y_obs, y_rep_mean, pch = 19, col = "steelblue", cex = 0.7,
                     main = "Posterior Predictive Check",
                     xlab = "Observed Y", ylab = "Predicted Y (mean)")
      graphics::abline(0, 1, col = "red", lwd = 2)
      graphics::grid()
    }
  }
  
  # ===========================================================================
  # PCA Plots
  # ===========================================================================
  
  if ("pca" %in% which && !is.null(x$pca)) {
    
    pca <- x$pca
    p <- x$data$p
    var_names <- colnames(x$data$X)
    if (is.null(var_names)) var_names <- paste0("X", 1:p)
    
    graphics::par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))
    
    # Scree plot
    n_show <- min(10, length(pca$variance_explained))
    ve <- pca$variance_explained[1:n_show]
    
    graphics::barplot(ve * 100, names.arg = paste0("PC", 1:n_show),
                      col = c(rep("steelblue", pca$n_components), 
                              rep("gray70", n_show - pca$n_components)),
                      main = "Scree Plot",
                      xlab = "Principal Component", ylab = "Variance Explained (%)")
    graphics::lines(1:n_show - 0.5, cumsum(ve) * 100, type = "b", pch = 19, col = "red")
    graphics::axis(4, at = pretty(cumsum(ve) * 100), col.axis = "red")
    
    # PC1 loadings
    loadings_pc1 <- pca$loadings[, 1]
    ord <- order(abs(loadings_pc1), decreasing = TRUE)
    top_n <- min(20, p)
    ord <- ord[1:top_n]
    
    graphics::barplot(loadings_pc1[ord], names.arg = var_names[ord], horiz = TRUE,
                      col = ifelse(loadings_pc1[ord] > 0, "steelblue", "coral"),
                      main = "PC1 Loadings (Top 20)",
                      xlab = "Loading", las = 1, cex.names = 0.7)
    graphics::abline(v = 0, lty = 2)
    
    # PC1 vs PC2 scores
    if (pca$n_components >= 2) {
      graphics::plot(pca$scores[, 1], pca$scores[, 2], pch = 19, col = "steelblue",
                     main = "Score Plot (PC1 vs PC2)",
                     xlab = paste0("PC1 (", round(pca$variance_explained[1] * 100, 1), "%)"),
                     ylab = paste0("PC2 (", round(pca$variance_explained[2] * 100, 1), "%)"))
      graphics::grid()
    }
    
    # Entropy by component
    n_ent <- min(5, length(pca$entropy))
    graphics::barplot(pca$entropy[1:n_ent], names.arg = paste0("PC", 1:n_ent),
                      col = "darkgreen",
                      main = "Loading Entropy by Component",
                      xlab = "Principal Component", ylab = "Shannon Entropy")
    graphics::abline(h = log(p), lty = 2, col = "red")  # Max entropy line
    graphics::text(n_ent, log(p), "Max entropy", pos = 1, col = "red", cex = 0.8)
  }
  
  # ===========================================================================
  # DFM Plots
  # ===========================================================================
  
  if ("dfm" %in% which && !is.null(x$dfm)) {
    
    dfm <- x$dfm
    n <- x$data$n
    time_idx <- 1:n
    
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
    
    # Information criterion
    if (!is.null(dfm$ic_values)) {
      n_ic <- length(dfm$ic_values)
      graphics::plot(1:n_ic, dfm$ic_values, type = "b", pch = 19,
                     col = ifelse(1:n_ic == dfm$n_factors, "red", "steelblue"),
                     main = "Information Criterion",
                     xlab = "Number of Factors", ylab = "IC Value")
      graphics::points(dfm$n_factors, dfm$ic_values[dfm$n_factors], 
                       pch = 19, col = "red", cex = 2)
      graphics::grid()
    }
    
    # Factor time series
    if (!is.null(dfm$factors)) {
      n_plot <- min(3, ncol(dfm$factors))
      for (i in 1:n_plot) {
        graphics::plot(time_idx, dfm$factors[, i], type = "l", col = "steelblue", lwd = 1.5,
                       main = paste("Factor", i),
                       xlab = "Time", ylab = paste0("F", i))
        graphics::abline(h = 0, lty = 2, col = "gray40")
        graphics::grid()
      }
    }
    
    # Variance explained
    if (!is.null(dfm$variance_explained)) {
      n_show <- length(dfm$variance_explained)
      graphics::barplot(dfm$variance_explained * 100, 
                        names.arg = paste0("F", 1:n_show),
                        col = "steelblue",
                        main = "Variance Explained by Factor",
                        xlab = "Factor", ylab = "Variance Explained (%)")
    }
  }
  
  # Reset par
  graphics::par(mfrow = c(1, 1))
  
  invisible(x)
}


# =============================================================================
# Helper: Null-coalescing operator
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @title Interactive Plot for Signal Analysis
#' @description
#' Generates an interactive dashboard using plotly to explore the results of
#' the signal analysis. Allows zooming, panning, and toggling traces.
#'
#' @param x An object of class \code{signal_analysis}
#'
#' @return A plotly object (HTML widget).
#' @export
iplot <- function(x) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this function. Please install it.")
  }
  
  Y <- x$data$Y_std
  n <- length(Y)
  time_idx <- 1:n
  
  plot_list <- list()
  
  p_base <- plotly::plot_ly(x = time_idx, y = Y, name = "Original Series", 
                            type = 'scatter', mode = 'lines',
                            line = list(color = 'black', width = 1, opacity = 0.3))
  
  if (!is.null(x$filters$wavelet)) {
    wavelet_trend <- x$filters$wavelet$combined 
    
    p_wavelet <- plotly::add_trace(p_base, y = wavelet_trend, name = "Wavelet Trend",
                                   line = list(color = '#00cc96', width = 2, opacity = 1))
    p_wavelet <- plotly::layout(p_wavelet, yaxis = list(title = "Wavelet"))
    
    plot_list[[length(plot_list) + 1]] <- p_wavelet
  }
  
  if (!is.null(x$filters$hpgc)) {
    p_hpgc_trend <- plotly::add_trace(p_base, y = x$filters$hpgc$trend, name = "HP-GC Trend",
                                      line = list(color = '#636efa', width = 2, opacity = 1))
    p_hpgc_trend <- plotly::layout(p_hpgc_trend, yaxis = list(title = "HP-GC Trend"))
    
    plot_list[[length(plot_list) + 1]] <- p_hpgc_trend
    
    p_hpgc_cycle <- plotly::plot_ly(x = time_idx, y = x$filters$hpgc$cycle, 
                                    name = "HP-GC Cycle", type = 'scatter', mode = 'lines',
                                    line = list(color = '#ef553b', width = 1.5))
    p_hpgc_cycle <- plotly::layout(p_hpgc_cycle, yaxis = list(title = "HP-GC Cycle"))
    
    p_hpgc_cycle <- plotly::add_lines(p_hpgc_cycle, x = c(1, n), y = c(0, 0), 
                                      line = list(color = "gray", dash = "dash", width = 1),
                                      showlegend = FALSE, inherit = FALSE)
    
    plot_list[[length(plot_list) + 1]] <- p_hpgc_cycle
  }
  
  if (!is.null(x$filters$emd)) {
    p_emd <- plotly::add_trace(p_base, y = x$filters$emd$trend, name = "EMD Trend",
                               line = list(color = '#ab63fa', width = 2, opacity = 1))
    p_emd <- plotly::layout(p_emd, yaxis = list(title = "EMD"))
    
    plot_list[[length(plot_list) + 1]] <- p_emd
  }
  
  if (length(plot_list) > 0) {
    
    final_plot <- plotly::subplot(plot_list, nrows = length(plot_list), 
                                  shareX = TRUE, titleY = TRUE)
    
    final_plot <- plotly::layout(final_plot,
                                 title = "SignalY: Interactive Analysis Dashboard",
                                 xaxis = list(
                                   title = "Time Index",
                                   rangeslider = list(visible = TRUE)
                                 ),
                                 hovermode = "x unified",
                                 legend = list(orientation = "h", x = 0, y = 1.1)
    )
    
    return(final_plot)
    
  } else {
    message("No filtering results found to plot.")
    return(invisible(NULL))
  }
}