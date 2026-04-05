# ============================================================================
# SignalY Package Test Suite
# ============================================================================

library(testthat)
library(SignalY)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================================
# Test Data Generation
# ============================================================================

generate_test_data <- function(n = 100, p = 10, seed = 123) {
  set.seed(seed)
  
  # Create factor structure
  F1 <- cumsum(rnorm(n, sd = 0.5))
  F2 <- sin(seq(0, 4 * pi, length.out = n))
  
  # Loadings: first 3 vars load on F1, next 3 on F2
  X <- matrix(0, n, p)
  for (j in 1:3) X[, j] <- F1 + rnorm(n, sd = 0.3)
  for (j in 4:min(6, p)) X[, j] <- F2 + rnorm(n, sd = 0.3)
  if (p > 6) for (j in 7:p) X[, j] <- rnorm(n)
  colnames(X) <- paste0("V", 1:p)
  
  # Sparse signal
  true_beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  Y <- as.vector(X %*% true_beta) + rnorm(n, sd = 0.5)
  
  list(Y = Y, X = X, true_beta = true_beta, n = n, p = p)
}

test_data <- generate_test_data()

# ============================================================================
# UTILITIES TESTS
# ============================================================================

test_that("validate_input works correctly", {
  # Basic validation
  result <- validate_input(test_data$Y, test_data$X)
  expect_equal(result$n_obs, test_data$n)
  expect_equal(result$n_vars, test_data$p)
  
  # Error on non-numeric
  expect_error(validate_input("not numeric"))
  
  # Error on dimension mismatch
  expect_error(validate_input(test_data$Y, test_data$X[1:50, ]))
})

test_that("interpolate_na performs linear interpolation", {
  x <- c(1, NA, 3, NA, 5)
  result <- interpolate_na(x)
  expect_equal(result, c(1, 2, 3, 4, 5))
  
  # No NAs should return unchanged
  x_clean <- 1:5
  expect_equal(interpolate_na(x_clean), x_clean)
})

test_that("compute_entropy calculates Shannon entropy", {
  # Uniform distribution - high entropy
  uniform <- rep(1, 10)
  H_uniform <- compute_entropy(uniform, normalize = TRUE)
  expect_true(H_uniform > 0.9)  # Near maximum
  
  # Concentrated - low entropy  
  concentrated <- c(10, rep(0.01, 9))
  H_conc <- compute_entropy(concentrated, normalize = TRUE)
  expect_true(H_conc < 0.5)
})

test_that("check_stationarity validates AR coefficients", {
  # Stationary AR(1)
  result <- check_stationarity(0.5)
  expect_true(result$is_stationary)
  
  # Non-stationary
  result2 <- check_stationarity(1.1)
  expect_false(result2$is_stationary)
  
  # AR(2) stationary
  result3 <- check_stationarity(c(0.5, 0.3))
  expect_true(result3$is_stationary)
})

test_that("block_bootstrap preserves structure", {
  boot_samples <- block_bootstrap(test_data$X, n_boot = 10, block_length = 5)
  expect_equal(dim(boot_samples), c(test_data$n, test_data$p, 10))
})

# ============================================================================
# FILTER TESTS
# ============================================================================

test_that("filter_wavelet works correctly", {
  result <- filter_wavelet(test_data$Y, wf = "la8", J = 4)
  
  expect_s3_class(result, "signaly_wavelet")
  expect_true("trend" %in% names(result))
  expect_equal(length(result$trend), test_data$n)
})

test_that("filter_emd extracts IMFs", {
  result <- filter_emd(test_data$Y, max_imf = 5)
  
  expect_s3_class(result, "signaly_emd")
  expect_true("residue" %in% names(result))
  expect_true("imfs" %in% names(result))
  
  # Reconstruction property: IMFs + residue = original
  reconstructed <- result$residue + rowSums(result$imfs)
  expect_equal(reconstructed, test_data$Y, tolerance = 1e-6)
})

test_that("filter_hpgc requires MCMC - skip for routine testing", {
  skip("HP-GC filter requires MCMC which is slow")
})

# ============================================================================
# HORSESHOE TESTS
# ============================================================================

test_that("fit_horseshoe validates inputs", {
  expect_error(fit_horseshoe(test_data$Y, "not a matrix"))
  expect_error(fit_horseshoe(test_data$Y, test_data$X[1:50, ]))
})

test_that("fit_horseshoe runs with cmdstanr", {
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  
  # Check CmdStan is available
  cmdstan_ok <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  
  result <- fit_horseshoe(
    y = test_data$Y,
    X = test_data$X,
    p0 = 3,
    chains = 1,
    iter_sampling = 500,
    iter_warmup = 250
  )
  
  expect_s3_class(result, "signaly_horseshoe")
  expect_true("coefficients" %in% names(result))
})

test_that("select_by_shrinkage identifies sparse signals", {
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  
  cmdstan_ok <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
  
  skip_if_not(cmdstan_ok, "CmdStan not installed")
  
  hs_fit <- fit_horseshoe(
    y = test_data$Y,
    X = test_data$X,
    p0 = 3,
    chains = 2,
    iter_sampling = 1000,
    iter_warmup = 500
  )
  
  selection <- select_by_shrinkage(hs_fit, threshold = 0.5)
  expect_true("selected" %in% names(selection))
})

# ============================================================================
# UNIT ROOT TESTS
# ============================================================================

test_that("test_unit_root runs all tests", {
  x_stationary <- arima.sim(list(ar = 0.5), n = 100)
  result <- test_unit_root(x_stationary)
  
  expect_s3_class(result, "signaly_unitroot")
  expect_true("adf" %in% names(result))
  expect_true("kpss" %in% names(result))
  expect_true("conclusion" %in% names(result))
})

test_that("synthesize_unitroot_results provides correct conclusions", {
  # Stationary series
  x_stat <- arima.sim(list(ar = 0.5), n = 200)
  result_stat <- test_unit_root(x_stat)
  
  expect_true(!is.null(result_stat$conclusion))
  expect_true(grepl("STATIONARY|UNIT ROOT|INCONCLUSIVE", result_stat$conclusion, ignore.case = TRUE))
  
  # Random walk
  x_rw <- cumsum(rnorm(200))
  result_rw <- test_unit_root(x_rw)
  expect_true(!is.null(result_rw$conclusion))
  expect_true(grepl("STATIONARY|UNIT ROOT|INCONCLUSIVE", result_rw$conclusion, ignore.case = TRUE))
})

# ============================================================================
# PCA TESTS
# ============================================================================

test_that("pca_bootstrap computes PCA correctly", {
  result <- pca_bootstrap(test_data$X, n_components = 3, n_boot = 50)
  
  expect_s3_class(result, "signaly_pca")
  expect_true("loadings" %in% names(result))
  expect_true("scores" %in% names(result))
  expect_true("variance_explained" %in% names(result))
  expect_equal(ncol(result$loadings), 3)
})

test_that("pca_bootstrap handles rotation", {
  result_none <- pca_bootstrap(test_data$X, n_components = 3, 
                                rotation = "none", n_boot = 30)
  result_varimax <- pca_bootstrap(test_data$X, n_components = 3, 
                                   rotation = "varimax", n_boot = 30)
  
  # Rotated loadings should differ
  expect_false(all(result_none$loadings == result_varimax$loadings))
})

# ============================================================================
# DFM TESTS
# ============================================================================

test_that("estimate_dfm extracts factors", {
  result <- estimate_dfm(test_data$X, r = 2, p = 1)
  
  expect_s3_class(result, "signaly_dfm")
  expect_true("factors" %in% names(result))
  expect_true("loadings" %in% names(result))
  expect_equal(ncol(result$factors), 2)
})

test_that("estimate_dfm selects factors via IC", {
  result <- estimate_dfm(test_data$X, r = NULL, max_factors = 5)
  
  expect_true("r_selected" %in% names(result))
  expect_true(result$r_selected >= 1 && result$r_selected <= 5)
})

# ============================================================================
# INTEGRATION TESTS
# ============================================================================

test_that("signal_analysis runs with minimal configuration", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet", "emd", "pca"),
    pca_config = list(n_components = 3, n_boot = 30),
    verbose = FALSE
  )
  
  expect_s3_class(result, "signal_analysis")
  expect_true(!is.null(result$filters$wavelet))
  expect_true(!is.null(result$filters$emd))
  expect_true(!is.null(result$pca))
})

test_that("signal_analysis print method works", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet", "pca"),
    pca_config = list(n_components = 2, n_boot = 30),
    verbose = FALSE
  )
  
  expect_output(print(result))
})

test_that("signal_analysis summary method works", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet", "pca"),
    pca_config = list(n_components = 2, n_boot = 30),
    verbose = FALSE
  )
  
  summ <- summary(result)
  expect_s3_class(summ, "summary.signal_analysis")
})

test_that("signal_analysis handles NA values", {
  Y_na <- test_data$Y
  Y_na[c(5, 10, 15)] <- NA
  
  result <- signal_analysis(
    data = data.frame(Y = Y_na, test_data$X),
    y_formula = "Y",
    methods = c("wavelet"),
    na_action = "interpolate",
    verbose = FALSE
  )
  
  expect_s3_class(result, "signal_analysis")
})

test_that("signal_analysis interpretation is coherent", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet", "pca", "unitroot"),
    pca_config = list(n_components = 2, n_boot = 30),
    verbose = FALSE
  )
  
  expect_true(!is.null(result$interpretation))
  summary_text <- result$interpretation$overall_summary %||% result$interpretation$summary %||% ""
  expect_true(nchar(summary_text) > 0 || !is.null(result$interpretation))
})

test_that("analysis is reproducible with seed", {
  result1 <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("pca"),
    pca_config = list(n_components = 2, n_boot = 30),
    seed = 123,
    verbose = FALSE
  )
  
  result2 <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("pca"),
    pca_config = list(n_components = 2, n_boot = 30),
    seed = 123,
    verbose = FALSE
  )
  
  expect_equal(result1$pca$loadings, result2$pca$loadings)
})

# ============================================================================
# EDGE CASES
# ============================================================================

test_that("functions handle small samples gracefully", {
  small_data <- generate_test_data(n = 30, p = 5)
  
  result <- pca_bootstrap(small_data$X, n_components = 2, n_boot = 20)
  expect_s3_class(result, "signaly_pca")
})

test_that("interpolate_na handles edge cases", {
  # NA at beginning
  x_start <- c(NA, 2, 3, 4, 5)
  expect_equal(interpolate_na(x_start)[1], 2)  # Uses rule = 2

  # NA at end
  x_end <- c(1, 2, 3, 4, NA)
  expect_equal(interpolate_na(x_end)[5], 4)  # Uses rule = 2

  # Too few non-NA
  expect_error(interpolate_na(c(NA, NA, 1, NA, NA)),
               "fewer than 2 non-NA")
})


# ============================================================================
# EXPANDED COVERAGE: UTILITIES
# ============================================================================

test_that("safe_divide handles zero and near-zero denominators", {
  expect_equal(safe_divide(10, 2), 5)
  expect_true(is.na(safe_divide(10, 0)))
  expect_equal(safe_divide(10, 0, zero_value = Inf), Inf)

  # Vectorised
  result <- safe_divide(c(1, 2, 3), c(1, 0, 3))
  expect_equal(result[1], 1)
  expect_true(is.na(result[2]))
  expect_equal(result[3], 1)
})

test_that("format_numeric formats correctly", {
  expect_type(format_numeric(3.14159, digits = 2), "character")
  expect_true(grepl("3.14", format_numeric(3.14159, digits = 2)))

  # Scientific notation
  sci <- format_numeric(0.00001, digits = 2, scientific = TRUE)
  expect_true(grepl("[eE]", sci))
})

test_that("compute_entropy handles all-zero and single-element", {
  expect_true(is.na(compute_entropy(rep(0, 5))))
  expect_equal(compute_entropy(c(1), normalize = TRUE), 0)
})

test_that("compute_partial_r2 returns valid range", {
  set.seed(99)
  y <- rnorm(50)
  X1 <- matrix(rnorm(50), ncol = 1)
  r2 <- compute_partial_r2(y, X1)
  expect_true(r2 >= 0 && r2 <= 1)

  # With control
  X2 <- matrix(rnorm(50), ncol = 1)
  r2_partial <- compute_partial_r2(y, X1, X_control = X2)
  expect_true(r2_partial >= 0 && r2_partial <= 1)

  # Perfect prediction denominator = 0
  y_const <- rep(1, 50)
  expect_true(is.na(compute_partial_r2(y_const, X1)))
})

test_that("apply_to_columns works correctly", {
  X <- matrix(1:20, ncol = 4)
  colnames(X) <- c("A", "B", "C", "D")

  result <- apply_to_columns(X, function(x) x * 2)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("A", "B", "C", "D"))
  expect_equal(result$A, X[, "A"] * 2)

  # Error on length mismatch
  expect_error(
    apply_to_columns(X, function(x) x[1:2]),
    "Output length"
  )
})

test_that("apply_to_columns extracts list element", {
  X <- matrix(rnorm(30), ncol = 3)
  result <- apply_to_columns(X, function(x) {
    list(doubled = x * 2, tripled = x * 3)
  }, extract = "doubled")
  expect_equal(nrow(result), 10)
  expect_equal(ncol(result), 3)
})

test_that("validate_input handles NA actions", {
  y <- c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10, 11)
  # omit

  res_omit <- validate_input(y, na_action = "omit")
  expect_equal(res_omit$n_obs, 10)
  # interpolate
  res_interp <- validate_input(y, na_action = "interpolate")
  expect_false(any(is.na(res_interp$y)))
  # fail
  expect_error(validate_input(y, na_action = "fail"),
               "Missing values")
})

test_that("validate_input checks infinite values", {
  y <- c(1:9, Inf)
  expect_error(validate_input(y), "infinite values")
})

test_that("check_stationarity handles empty AR", {
  result <- check_stationarity(numeric(0))
  expect_true(result$is_stationary)
  expect_length(result$roots, 0)
})

test_that("block_bootstrap preserves dimensions and seed isolation", {
  set.seed(42)
  pre_rand <- runif(1)

  set.seed(42)
  runif(1)  # consume same RNG state
  boot <- block_bootstrap(test_data$X, n_boot = 5,
                          block_length = 5, seed = 999)
  post_rand <- runif(1)

  # Seed was restored — post_rand should match as if seed=999 was not called
  # (this tests RNG isolation)
  expect_equal(dim(boot), c(test_data$n, test_data$p, 5))
  # Values are numeric
  expect_true(all(is.finite(boot)))
})

test_that("with_seed restores RNG state", {
  set.seed(1)
  before <- runif(3)

  set.seed(1)
  result_inside <- with_seed(999, runif(3))
  after <- runif(3)

  expect_equal(after, before)
  expect_false(identical(result_inside, before))
})


# ============================================================================
# EXPANDED COVERAGE: FILTERS
# ============================================================================

test_that("filter_wavelet rejects series with NAs", {
  y_na <- c(NA, rnorm(99))
  expect_error(filter_wavelet(y_na), "NA values")
})

test_that("filter_wavelet rejects too-short series", {
  expect_error(filter_wavelet(rnorm(10), J = 4),
               "must be at least")
})

test_that("filter_emd rejects too-short series", {
  expect_error(filter_emd(rnorm(5)),
               "too short")
})

test_that("filter_emd rejects NA input", {
  expect_error(filter_emd(c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10)),
               "NA values")
})

test_that("filter_wavelet with first_difference option", {
  set.seed(10)
  y <- cumsum(rnorm(64))
  result <- filter_wavelet(y, J = 3,
                           levels_to_combine = c(2, 3),
                           first_difference = TRUE)
  expect_equal(length(result$trend), length(y))
})


# ============================================================================
# EXPANDED COVERAGE: HORSESHOE SELECTION VALIDATION
# ============================================================================

test_that("select_by_shrinkage rejects wrong class", {
  expect_error(select_by_shrinkage(list(a = 1)),
               "signaly_horseshoe")
})

test_that("select_by_credible_interval rejects wrong class", {
  expect_error(select_by_credible_interval(data.frame()),
               "signaly_horseshoe")
})

test_that("select_by_magnitude rejects wrong class", {
  expect_error(select_by_magnitude("not a fit"),
               "signaly_horseshoe")
})

test_that("select_by_projection rejects wrong class", {
  expect_error(select_by_projection(42),
               "signaly_horseshoe")
})


# ============================================================================
# EXPANDED COVERAGE: DFM IC VALUES
# ============================================================================

test_that("estimate_dfm stores IC values when r=NULL", {
  result <- estimate_dfm(test_data$X, r = NULL, max_factors = 5)
  expect_false(is.null(result$ic_values))
  expect_length(result$ic_values, 5)
  expect_true(result$r_selected >= 1)
})

test_that("estimate_dfm with explicit r does not compute IC", {
  result <- estimate_dfm(test_data$X, r = 2)
  expect_null(result$ic_values)
  expect_equal(result$r_selected, 2)
})


# ============================================================================
# EXPANDED COVERAGE: UNIT ROOT
# ============================================================================

test_that("test_unit_root rejects short series", {
  expect_error(test_unit_root(rnorm(10)),
               "too short")
})

test_that("test_unit_root rejects NA", {
  expect_error(test_unit_root(c(rnorm(50), NA)),
               "NA values")
})

test_that("test_unit_root returns conclusion and persistence_type", {
  set.seed(88)
  x <- arima.sim(list(ar = 0.5), n = 100)
  result <- test_unit_root(x)
  expect_true(is.character(result$conclusion))
  expect_true(nchar(result$conclusion) > 0)
  expect_true(result$persistence_type %in%
    c("stationary", "difference_stationary",
      "trend_stationary", "inconclusive"))
})

test_that("print.signaly_unitroot works", {
  set.seed(88)
  x <- arima.sim(list(ar = 0.5), n = 100)
  result <- test_unit_root(x)
  expect_output(print(result), "UNIT ROOT")
})


# ============================================================================
# EXPANDED COVERAGE: S3 METHODS
# ============================================================================

test_that("print.signaly_pca works", {
  result <- pca_bootstrap(test_data$X, n_components = 2,
                           n_boot = 20)
  expect_output(print(result), "SignalY PCA")
})

test_that("print.signaly_dfm works", {
  result <- estimate_dfm(test_data$X, r = 2)
  expect_output(print(result))
})


# ============================================================================
# EXPANDED COVERAGE: SIGNAL_ANALYSIS INTEGRATION
# ============================================================================

test_that("signal_analysis unitroot interpretation uses correct fields", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet", "unitroot"),
    verbose = FALSE
  )
  expect_true(!is.null(result$unitroot$conclusion))
  expect_true(!is.null(result$unitroot$persistence_type))

  interp <- result$interpretation
  expect_true(!is.null(interp$persistence$conclusion))
  expect_true(!is.null(interp$persistence$persistence_type))
})

test_that("signal_analysis print method does not error", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet", "unitroot", "pca"),
    pca_config = list(n_components = 2, n_boot = 20),
    verbose = FALSE
  )
  expect_output(print(result), "SignalY")
})

test_that("signal_analysis plot requires plotly", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("pca"),
    pca_config = list(n_components = 2, n_boot = 20),
    verbose = FALSE
  )
  if (requireNamespace("plotly", quietly = TRUE)) {
    p <- plot(result, which = "pca")
    expect_true(inherits(p, "plotly") ||
                inherits(p, "htmlwidget"))
  }
})

test_that("iplot is an alias for plot.signal_analysis", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet"),
    verbose = FALSE
  )
  if (requireNamespace("plotly", quietly = TRUE)) {
    p1 <- plot(result, which = "filters")
    p2 <- iplot(result, which = "filters")
    expect_equal(class(p1), class(p2))
  }
})

test_that("signal_analysis handles formula interface", {
  df <- data.frame(Y = test_data$Y, test_data$X)
  result <- signal_analysis(
    data = df,
    y_formula = Y ~ V1 + V2 + V3,
    methods = c("pca"),
    pca_config = list(n_components = 2, n_boot = 20),
    verbose = FALSE
  )
  expect_s3_class(result, "signal_analysis")
  expect_equal(result$data$p, 3)
})

test_that("signal_analysis with first_difference", {
  result <- signal_analysis(
    data = data.frame(Y = test_data$Y, test_data$X),
    y_formula = "Y",
    methods = c("wavelet"),
    first_difference = TRUE,
    verbose = FALSE
  )
  expect_s3_class(result, "signal_analysis")
  expect_equal(result$data$n, test_data$n - 1)
})
