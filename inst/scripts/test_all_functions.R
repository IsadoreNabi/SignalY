#!/usr/bin/env Rscript
# ===========================================================================
# SignalY — Comprehensive Function Test Script
# ===========================================================================
# Tests EVERY exported and key internal function with real computations.
# Run from the package root:
#   Rscript inst/scripts/test_all_functions.R
#
# Requirements:
#   - waveslim, EMD, urca (Imports — mandatory)
#   - plotly (Suggests — optional, tested if available)
#   - cmdstanr + CmdStan (Suggests — optional, tested if available)
# ===========================================================================

cat("\n")
cat("###########################################################################\n")
cat("#        SignalY: COMPREHENSIVE FUNCTION TEST SUITE                       #\n")
cat("###########################################################################\n\n")

# Definimos la ruta de la carpeta del paquete fuente
ruta_paquete <- "C:/Users/ROG/Documents/SignalY"

# Load package from source (works without formal installation)
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(path = ruta_paquete, quiet = TRUE)
  
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(path = ruta_paquete, quiet = TRUE)
  
} else {
  # IMPORTANTE: library() solo funciona con paquetes ya instalados.
  # Si el paquete está "instalado" en esa carpeta específica (y no es solo código fuente),
  # puedes indicarle a R dónde buscarlo usando 'lib.loc':
  library(SignalY, lib.loc = "C:/Users/ROG/Documents")
}

passed  <- 0L
failed  <- 0L
skipped <- 0L
errors  <- character(0)

# ---------------------------------------------------------------------------
# Helper: run a single test block, capture errors
# ---------------------------------------------------------------------------
run_test <- function(name, expr) {
  cat(sprintf("  [TEST] %-58s", name))
  tryCatch({
    result <- withCallingHandlers(
      eval.parent(substitute(expr)),
      warning = function(w) invokeRestart("muffleWarning")
    )
    cat(" PASS\n")
    passed <<- passed + 1L
    invisible(result)
  }, error = function(e) {
    cat(" FAIL\n")
    msg <- sprintf("    >> %s: %s", name, conditionMessage(e))
    cat(msg, "\n")
    errors <<- c(errors, msg)
    failed <<- failed + 1L
    invisible(NULL)
  })
}

skip_test <- function(name, reason) {
  cat(sprintf("  [SKIP] %-58s %s\n", name, reason))
  skipped <<- skipped + 1L
}

# ===========================================================================
# 0. GENERATE TEST DATA
# ===========================================================================
cat("--- Generating test data -------------------------------------------------\n")

set.seed(42)
N   <- 120
P   <- 15
K_true <- 4

factors  <- matrix(rnorm(N * 2), N, 2)
loadings <- matrix(runif(P * 2, -1, 1), P, 2)
X <- factors %*% t(loadings) + matrix(rnorm(N * P, 0, 0.3), N, P)
colnames(X) <- paste0("V", seq_len(P))

beta_true <- c(rep(1.5, K_true), rep(0, P - K_true))
Y <- X %*% beta_true + rnorm(N, 0, 0.5)

panel <- data.frame(Y = as.numeric(Y), X)

y_ts <- cumsum(rnorm(N)) + sin(seq(0, 6 * pi, length.out = N))

y_stationary <- as.numeric(arima.sim(list(ar = 0.5), n = N))
y_nonstat    <- cumsum(rnorm(N))

cat("  N =", N, " P =", P, " true sparsity =", K_true, "/", P, "\n\n")

# ===========================================================================
# 1. UTILITY FUNCTIONS (internal)
# ===========================================================================
cat("=== 1. UTILITY FUNCTIONS =================================================\n")

run_test("safe_divide: normal division", {
  res <- SignalY:::safe_divide(10, 2)
  stopifnot(res == 5)
})

run_test("safe_divide: zero denominator -> NA", {
  res <- SignalY:::safe_divide(10, 0)
  stopifnot(is.na(res))
})

run_test("safe_divide: zero denominator -> custom value", {
  res <- SignalY:::safe_divide(10, 0, zero_value = 0)
  stopifnot(res == 0)
})

run_test("safe_divide: vectorized", {
  res <- SignalY:::safe_divide(c(1, 2, 3), c(1, 0, 3))
  stopifnot(res[1] == 1, is.na(res[2]), res[3] == 1)
})

run_test("format_numeric: basic", {
  res <- SignalY:::format_numeric(3.14159, digits = 2)
  stopifnot(is.character(res), grepl("3.14", res))
})

run_test("format_numeric: scientific notation", {
  res <- SignalY:::format_numeric(0.000123, digits = 2, scientific = TRUE)
  stopifnot(grepl("e|E", res))
})

run_test("compute_entropy: uniform distribution -> max entropy", {
  e <- compute_entropy(rep(1, 10), normalize = TRUE)
  stopifnot(abs(e - 1) < 1e-10)
})

run_test("compute_entropy: concentrated -> low entropy", {
  e <- compute_entropy(c(10, rep(0.01, 9)), normalize = TRUE)
  stopifnot(e < 0.3)
})

run_test("compute_entropy: all zeros -> NA with warning", {
  e <- suppressWarnings(compute_entropy(rep(0, 5)))
  stopifnot(is.na(e))
})

run_test("compute_entropy: single element", {
  e <- compute_entropy(5, normalize = TRUE)
  stopifnot(is.numeric(e), !is.na(e))
})

run_test("compute_partial_r2: simple R2", {
  r2 <- SignalY:::compute_partial_r2(Y, X[, 1:K_true])
  stopifnot(is.numeric(r2), r2 > 0, r2 <= 1)
})

run_test("compute_partial_r2: partial R2 with controls", {
  r2 <- SignalY:::compute_partial_r2(
    Y,
    X_interest = X[, 1:2],
    X_control  = X[, 3:5]
  )
  stopifnot(is.numeric(r2), r2 >= 0, r2 <= 1)
})

run_test("interpolate_na: no NAs -> identity", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  res <- SignalY:::interpolate_na(x)
  stopifnot(all(abs(res - x) < 1e-12))
})

run_test("interpolate_na: interior NAs", {
  x <- c(1, NA, 3, NA, 5)
  res <- SignalY:::interpolate_na(x)
  stopifnot(!any(is.na(res)), abs(res[2] - 2) < 1e-10)
})

run_test("interpolate_na: leading/trailing NAs", {
  x <- c(NA, NA, 3, 4, NA)
  res <- SignalY:::interpolate_na(x)
  stopifnot(!any(is.na(res)), res[1] == 3, res[5] == 4)
})

run_test("check_stationarity: stationary AR(1)", {
  res <- SignalY:::check_stationarity(0.5)
  stopifnot(res$is_stationary)
})

run_test("check_stationarity: unit root AR(1)", {
  res <- SignalY:::check_stationarity(1.0)
  stopifnot(!res$is_stationary)
})

run_test("check_stationarity: explosive AR(1)", {
  res <- SignalY:::check_stationarity(1.2)
  stopifnot(!res$is_stationary)
})

run_test("check_stationarity: empty coefficients", {
  res <- SignalY:::check_stationarity(numeric(0))
  stopifnot(res$is_stationary, length(res$roots) == 0)
})

run_test("check_stationarity: AR(2) stationary", {
  res <- SignalY:::check_stationarity(c(0.5, -0.2))
  stopifnot(res$is_stationary)
})

run_test("validate_input: good data", {
  res <- SignalY:::validate_input(Y, X)
  stopifnot(res$n_obs == N, res$n_vars == P)
})

run_test("validate_input: too few observations", {
  ok <- tryCatch({
    SignalY:::validate_input(1:3, min_obs = 10)
    FALSE
  }, error = function(e) grepl("Insufficient", e$message))
  stopifnot(ok)
})

run_test("validate_input: infinite values rejected", {
  y_inf <- c(1:19, Inf)
  ok <- tryCatch({
    SignalY:::validate_input(y_inf, min_obs = 10)
    FALSE
  }, error = function(e) grepl("infinite", e$message))
  stopifnot(ok)
})

run_test("validate_input: NA action = omit", {
  y_na <- Y
  y_na[c(1, 5)] <- NA
  res <- SignalY:::validate_input(y_na, na_action = "omit", verbose = FALSE)
  stopifnot(res$n_obs == N - 2)
})

run_test("validate_input: NA action = interpolate", {
  y_na <- Y
  y_na[c(1, 5)] <- NA
  res <- SignalY:::validate_input(y_na, na_action = "interpolate", verbose = FALSE)
  stopifnot(res$n_obs == N, !any(is.na(res$y)))
})

run_test("validate_input: NA action = fail", {
  y_na <- Y; y_na[1] <- NA
  ok <- tryCatch({
    SignalY:::validate_input(y_na, na_action = "fail")
    FALSE
  }, error = function(e) grepl("Missing values", e$message))
  stopifnot(ok)
})

run_test("validate_input: dimension mismatch X vs y", {
  ok <- tryCatch({
    SignalY:::validate_input(Y[1:50], X)
    FALSE
  }, error = function(e) grepl("mismatch", e$message, ignore.case = TRUE))
  stopifnot(ok)
})

run_test("validate_horseshoe_fit: rejects wrong class", {
  ok <- tryCatch({
    SignalY:::validate_horseshoe_fit(list(a = 1))
    FALSE
  }, error = function(e) grepl("signaly_horseshoe", e$message))
  stopifnot(ok)
})

run_test("safe_cores: returns positive integer", {
  cores <- SignalY:::safe_cores(reserve = 1L)
  stopifnot(is.numeric(cores), cores >= 1L)
})

run_test("with_seed: preserves RNG state", {
  set.seed(99)
  before <- rnorm(1)
  set.seed(99)
  SignalY:::with_seed(12345, rnorm(10))
  after <- rnorm(1)
  stopifnot(abs(before - after) < 1e-15)
})

run_test("with_seed: produces deterministic output", {
  a <- SignalY:::with_seed(777, sample(1:100, 5))
  b <- SignalY:::with_seed(777, sample(1:100, 5))
  stopifnot(identical(a, b))
})

run_test("block_bootstrap: dimensions correct", {
  bb <- SignalY:::block_bootstrap(X[, 1:3], n_boot = 5, block_length = 10)
  stopifnot(
    is.array(bb),
    dim(bb)[1] == N,
    dim(bb)[2] == 3,
    dim(bb)[3] == 5
  )
})

run_test("block_bootstrap: seed reproducibility", {
  b1 <- SignalY:::block_bootstrap(X[, 1:2], n_boot = 3, seed = 42)
  b2 <- SignalY:::block_bootstrap(X[, 1:2], n_boot = 3, seed = 42)
  stopifnot(identical(b1, b2))
})

run_test("block_bootstrap: RNG isolation (seed does not leak)", {
  set.seed(99)
  before <- rnorm(1)
  set.seed(99)
  SignalY:::block_bootstrap(X[, 1:2], n_boot = 3, seed = 42)
  after <- rnorm(1)
  stopifnot(abs(before - after) < 1e-15)
})

# ===========================================================================
# 2. apply_to_columns (exported)
# ===========================================================================
cat("\n=== 2. apply_to_columns ==================================================\n")

run_test("apply_to_columns: simple function", {
  res <- apply_to_columns(X[, 1:4], cumsum)
  stopifnot(is.data.frame(res), nrow(res) == N, ncol(res) == 4)
})

run_test("apply_to_columns: with extract parameter", {
  res <- apply_to_columns(
    X[, 1:3],
    function(x) list(trend = cumsum(x), other = "foo"),
    extract = "trend"
  )
  stopifnot(is.data.frame(res), nrow(res) == N)
})

run_test("apply_to_columns: mismatched output length -> error", {
  ok <- tryCatch({
    apply_to_columns(X[, 1:2], function(x) x[1:5])
    FALSE
  }, error = function(e) grepl("length", e$message))
  stopifnot(ok)
})

run_test("apply_to_columns: preserves column names", {
  M <- X[, 1:3]
  colnames(M) <- c("alpha", "beta", "gamma")
  res <- apply_to_columns(M, cumsum)
  stopifnot(identical(names(res), c("alpha", "beta", "gamma")))
})

# ===========================================================================
# 3. FILTER: WAVELET (filter_wavelet)
# ===========================================================================
cat("\n=== 3. WAVELET FILTER ====================================================\n")

run_test("filter_wavelet: basic decomposition", {
  res <- filter_wavelet(y_ts, wf = "la8", J = 4, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_wavelet"),
    length(res$trend) == length(y_ts),
    !is.null(res$mra),
    !is.null(res$diagnostics)
  )
})

run_test("filter_wavelet: trend + details = original", {
  res <- filter_wavelet(y_ts, wf = "la8", J = 4)
  reconstructed <- rowSums(res$detail_levels) + res$smooth_level
  stopifnot(max(abs(y_ts - reconstructed)) < 1e-10)
})

run_test("filter_wavelet: different filters", {
  for (wf in c("haar", "d4", "la8")) {
    res <- filter_wavelet(y_ts, wf = wf, J = 3, verbose = FALSE)
    stopifnot(length(res$trend) == length(y_ts))
  }
})

run_test("filter_wavelet: custom levels_to_combine", {
  res <- filter_wavelet(y_ts, J = 4, levels_to_combine = c(2, 3, 4))
  stopifnot(grepl("D2", res$combined_levels))
})

run_test("filter_wavelet: first_difference = TRUE", {
  res <- filter_wavelet(y_ts, J = 3, first_difference = TRUE, verbose = FALSE)
  stopifnot(length(res$trend) == length(y_ts))
})

run_test("filter_wavelet: rejects NA input", {
  y_bad <- y_ts; y_bad[5] <- NA
  ok <- tryCatch({
    filter_wavelet(y_bad)
    FALSE
  }, error = function(e) grepl("NA", e$message))
  stopifnot(ok)
})

run_test("filter_wavelet: rejects too-short series", {
  ok <- tryCatch({
    filter_wavelet(1:10, J = 4)
    FALSE
  }, error = function(e) grepl("length", e$message, ignore.case = TRUE))
  stopifnot(ok)
})

run_test("filter_wavelet: variance proportions sum to ~1", {
  res <- filter_wavelet(y_ts, J = 4)
  diag <- res$diagnostics
  stopifnot(abs(sum(diag$variance_proportion) - 1) < 0.01)
})

# ===========================================================================
# 4. FILTER: EMD (filter_emd)
# ===========================================================================
cat("\n=== 4. EMD FILTER ========================================================\n")

run_test("filter_emd: basic decomposition", {
  res <- filter_emd(y_ts, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_emd"),
    length(res$trend) == length(y_ts),
    res$n_imfs >= 1,
    is.matrix(res$imfs)
  )
})

run_test("filter_emd: IMFs + residue = original", {
  res <- filter_emd(y_ts)
  reconstructed <- rowSums(res$imfs) + res$residue
  stopifnot(max(abs(y_ts - reconstructed)) < 1e-8)
})

run_test("filter_emd: rejects NA input", {
  y_bad <- y_ts; y_bad[3] <- NA
  ok <- tryCatch({
    filter_emd(y_bad)
    FALSE
  }, error = function(e) grepl("NA", e$message))
  stopifnot(ok)
})

run_test("filter_emd: rejects too-short series", {
  ok <- tryCatch({
    filter_emd(1:5)
    FALSE
  }, error = function(e) grepl("short", e$message, ignore.case = TRUE))
  stopifnot(ok)
})

run_test("filter_emd: max_imf limits extraction", {
  res <- filter_emd(y_ts, max_imf = 3)
  stopifnot(res$n_imfs <= 3)
})

run_test("filter_emd: different boundary conditions", {
  for (bd in c("periodic", "symmetric")) {
    res <- filter_emd(y_ts, boundary = bd, verbose = FALSE)
    stopifnot(length(res$trend) == N)
  }
})

# ===========================================================================
# 5. FILTER: HP-GC (filter_hpgc)
#    NOTE: Each MCMC run is memory-heavy. We fit ONCE and test all
#    properties on that single result, then free memory immediately.
# ===========================================================================
cat("\n=== 5. HP-GC FILTER ======================================================\n")

run_test("filter_hpgc: rejects NA", {
  y_bad <- y_ts; y_bad[10] <- NA
  ok <- tryCatch({
    filter_hpgc(y_bad)
    FALSE
  }, error = function(e) grepl("NA", e$message))
  stopifnot(ok)
})

run_test("filter_hpgc: rejects too-short series", {
  ok <- tryCatch({
    filter_hpgc(1:10)
    FALSE
  }, error = function(e) grepl("short", e$message, ignore.case = TRUE))
  stopifnot(ok)
})

run_test("filter_hpgc: weak prior (all properties)", {
  res <- filter_hpgc(y_ts, prior_config = "weak",
                     n_chains = 2, iterations = 2000,
                     burnin = 500, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_hpgc"),
    length(res$trend) == N,
    length(res$cycle) == N,
    length(res$trend_lower) == N,
    length(res$trend_upper) == N,
    all(res$trend_lower <= res$trend_upper),
    !is.null(res$diagnostics),
    !is.null(res$dic),
    cor(y_ts, res$trend + res$cycle) > 0.9
  )
  rm(res); gc(verbose = FALSE)
})

run_test("filter_hpgc: informative prior", {
  res <- filter_hpgc(y_ts, prior_config = "informative",
                     n_chains = 2, iterations = 2000,
                     burnin = 500, verbose = FALSE)
  stopifnot(length(res$trend) == N)
  rm(res); gc(verbose = FALSE)
})

run_test("filter_hpgc: empirical prior", {
  res <- filter_hpgc(y_ts, prior_config = "empirical",
                     n_chains = 2, iterations = 2000,
                     burnin = 500, verbose = FALSE)
  stopifnot(length(res$trend) == N)
  rm(res); gc(verbose = FALSE)
})

# ===========================================================================
# 6. FILTER: filter_all (combined)
# ===========================================================================
cat("\n=== 6. COMBINED FILTER (filter_all) ======================================\n")

run_test("filter_all: runs all 3 filters + correlation", {
  res <- filter_all(y_ts, hpgc_chains = 2, hpgc_iterations = 2000,
                    hpgc_burnin = 500, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_multifilter"),
    !is.null(res$wavelet),
    !is.null(res$emd),
    !is.null(res$hpgc),
    is.data.frame(res$comparison),
    nrow(res$comparison) == N,
    is.matrix(res$correlations),
    nrow(res$correlations) == 3,
    ncol(res$correlations) == 3
  )
  rm(res); gc(verbose = FALSE)
})

# ===========================================================================
# 7. PCA WITH BOOTSTRAP (pca_bootstrap)
# ===========================================================================
cat("\n=== 7. PCA WITH BOOTSTRAP ================================================\n")

run_test("pca_bootstrap: basic (auto components)", {
  res <- pca_bootstrap(X, n_boot = 50, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_pca"),
    res$n_components >= 1,
    is.matrix(res$loadings),
    is.matrix(res$scores),
    length(res$eigenvalues) == P,
    length(res$variance_explained) == P,
    is.matrix(res$significant_loadings),
    is.matrix(res$p_values),
    !is.null(res$entropy)
  )
})

run_test("pca_bootstrap: fixed n_components", {
  res <- pca_bootstrap(X, n_components = 3, n_boot = 50, verbose = FALSE)
  stopifnot(res$n_components == 3, ncol(res$loadings) == 3)
})

run_test("pca_bootstrap: varimax rotation", {
  res <- pca_bootstrap(X, n_components = 3, rotation = "varimax",
                       n_boot = 50, verbose = FALSE)
  stopifnot(nrow(res$loadings) == P, ncol(res$loadings) == 3)
})

run_test("pca_bootstrap: no rotation", {
  res <- pca_bootstrap(X, n_components = 2, rotation = "none",
                       n_boot = 50, verbose = FALSE)
  stopifnot(ncol(res$loadings) == 2)
})

run_test("pca_bootstrap: variance explained sums <= 1", {
  res <- pca_bootstrap(X, n_boot = 50, verbose = FALSE)
  stopifnot(sum(res$variance_explained) <= 1 + 1e-10)
  stopifnot(all(res$cumulative_variance <= 1 + 1e-10))
})

run_test("pca_bootstrap: entropy within [0, max]", {
  res <- pca_bootstrap(X, n_components = 2, n_boot = 50, verbose = FALSE)
  max_ent <- log(P, base = 2)
  stopifnot(all(res$entropy >= 0), all(res$entropy <= max_ent + 1))
})

run_test("pca_bootstrap: p_values in [0, 1]", {
  res <- pca_bootstrap(X, n_components = 2, n_boot = 50, verbose = FALSE)
  stopifnot(all(res$p_values >= 0), all(res$p_values <= 1))
})

run_test("pca_bootstrap: assignments data frame", {
  res <- pca_bootstrap(X, n_components = 2, n_boot = 50, verbose = FALSE)
  stopifnot(is.data.frame(res$assignments), nrow(res$assignments) == P)
})

run_test("pca_bootstrap: summary_by_component", {
  res <- pca_bootstrap(X, n_components = 2, n_boot = 50, verbose = FALSE)
  stopifnot(is.data.frame(res$summary_by_component))
})

run_test("pca_bootstrap: custom block_length", {
  res <- pca_bootstrap(X, n_boot = 30, block_length = 5, verbose = FALSE)
  stopifnot(res$n_components >= 1)
})

run_test("pca_bootstrap: rejects NA in X", {
  X_bad <- X; X_bad[1, 1] <- NA
  ok <- tryCatch({
    pca_bootstrap(X_bad, n_boot = 10)
    FALSE
  }, error = function(e) grepl("NA", e$message))
  stopifnot(ok)
})

run_test("print.signaly_pca: runs without error", {
  res <- pca_bootstrap(X, n_components = 2, n_boot = 30, verbose = FALSE)
  invisible(capture.output(print(res)))
})

# ===========================================================================
# 8. DYNAMIC FACTOR MODEL (estimate_dfm)
# ===========================================================================
cat("\n=== 8. DYNAMIC FACTOR MODEL ==============================================\n")

run_test("estimate_dfm: auto-select r via IC2", {
  res <- estimate_dfm(X, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_dfm"),
    res$r_selected >= 1,
    is.matrix(res$factors),
    is.matrix(res$loadings),
    nrow(res$factors) == N,
    ncol(res$factors) == res$r_selected,
    !is.null(res$ic_values)
  )
})

run_test("estimate_dfm: fixed r", {
  res <- estimate_dfm(X, r = 3, verbose = FALSE)
  stopifnot(
    res$r_selected == 3,
    ncol(res$factors) == 3,
    nrow(res$loadings) == P,
    ncol(res$loadings) == 3
  )
})

run_test("estimate_dfm: IC1 criterion", {
  res <- estimate_dfm(X, ic = "IC1", verbose = FALSE)
  stopifnot(res$r_selected >= 1)
})

run_test("estimate_dfm: IC3 criterion", {
  res <- estimate_dfm(X, ic = "IC3", verbose = FALSE)
  stopifnot(res$r_selected >= 1)
})

run_test("estimate_dfm: max_factors limits search", {
  res <- estimate_dfm(X, max_factors = 3, verbose = FALSE)
  stopifnot(res$r_selected <= 3)
})

run_test("estimate_dfm: VAR dynamics (p = 2)", {
  res <- estimate_dfm(X, r = 2, p = 2, verbose = FALSE)
  stopifnot(!is.null(res$var_coefficients))
})

run_test("estimate_dfm: ic_values stored when auto-selecting", {
  res <- estimate_dfm(X, verbose = FALSE)
  stopifnot(is.numeric(res$ic_values), length(res$ic_values) >= 1)
})

run_test("estimate_dfm: ic_values NULL when r fixed", {
  res <- estimate_dfm(X, r = 2, verbose = FALSE)
  stopifnot(is.null(res$ic_values))
})

run_test("estimate_dfm: fitted values and residuals", {
  res <- estimate_dfm(X, r = 2, verbose = FALSE)
  stopifnot(
    is.matrix(res$fitted_values),
    is.matrix(res$residuals),
    all(dim(res$fitted_values) == dim(X)),
    all(dim(res$residuals) == dim(X))
  )
})

run_test("estimate_dfm: no standardize", {
  res <- estimate_dfm(X, r = 2, standardize = FALSE, verbose = FALSE)
  stopifnot(ncol(res$factors) == 2)
})

run_test("print.signaly_dfm: runs without error", {
  res <- estimate_dfm(X, r = 2, verbose = FALSE)
  invisible(capture.output(print(res)))
})

# ===========================================================================
# 9. UNIT ROOT TESTS (test_unit_root)
# ===========================================================================
cat("\n=== 9. UNIT ROOT TESTS ===================================================\n")

run_test("test_unit_root: stationary series", {
  res <- test_unit_root(y_stationary, verbose = FALSE)
  stopifnot(
    inherits(res, "signaly_unitroot"),
    !is.null(res$conclusion),
    !is.null(res$persistence_type),
    is.data.frame(res$summary)
  )
  stopifnot(grepl("STATIONARY|stationary", res$conclusion, ignore.case = TRUE))
})

run_test("test_unit_root: nonstationary series", {
  res <- test_unit_root(y_nonstat, verbose = FALSE)
  stopifnot(
    !is.null(res$conclusion),
    !is.null(res$persistence_type)
  )
})

run_test("test_unit_root: ADF results present", {
  res <- test_unit_root(y_stationary, verbose = FALSE)
  stopifnot(!is.null(res$adf), length(res$adf) >= 1)
})

run_test("test_unit_root: ERS results present", {
  res <- test_unit_root(y_stationary, verbose = FALSE)
  stopifnot(!is.null(res$ers), length(res$ers) >= 1)
})

run_test("test_unit_root: KPSS results present", {
  res <- test_unit_root(y_stationary, verbose = FALSE)
  stopifnot(!is.null(res$kpss), length(res$kpss) >= 1)
})

run_test("test_unit_root: PP results present", {
  res <- test_unit_root(y_stationary, verbose = FALSE)
  stopifnot(!is.null(res$pp), length(res$pp) >= 1)
})

run_test("test_unit_root: custom max_lags", {
  res <- test_unit_root(y_stationary, max_lags = 3, verbose = FALSE)
  stopifnot(!is.null(res$conclusion))
})

run_test("test_unit_root: verbose prints", {
  out <- capture.output(test_unit_root(y_stationary, verbose = TRUE))
  stopifnot(any(grepl("UNIT ROOT|CONCLUSION|ADF", out)))
})

run_test("print.signaly_unitroot: runs without error", {
  res <- test_unit_root(y_stationary, verbose = FALSE)
  invisible(capture.output(print(res)))
})

run_test("test_unit_root: rejects short series", {
  ok <- tryCatch({
    test_unit_root(1:10, verbose = FALSE)
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
})

# ===========================================================================
# 10. HORSESHOE REGRESSION (fit_horseshoe + select_by_*)
# ===========================================================================
cat("\n=== 10. HORSESHOE REGRESSION =============================================\n")

has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
if (has_cmdstanr) {
  has_cmdstan <- tryCatch({
    p <- cmdstanr::cmdstan_path()
    nchar(p) > 0
  }, error = function(e) FALSE)
} else {
  has_cmdstan <- FALSE
}

if (has_cmdstanr && has_cmdstan) {

  run_test("fit_horseshoe: basic fit", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    stopifnot(
      inherits(hs, "signaly_horseshoe"),
      is.data.frame(hs$coefficients),
      nrow(hs$coefficients) == P,
      !is.null(hs$sparsity),
      !is.null(hs$diagnostics),
      !is.null(hs$fit),
      !is.null(hs$stan_data)
    )
  })

  run_test("fit_horseshoe: custom p0", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X, p0 = 3,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    stopifnot(inherits(hs, "signaly_horseshoe"))
  })

  run_test("fit_horseshoe: coefficients structure", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    stopifnot(
      "variable" %in% names(hs$coefficients),
      "beta_mean" %in% names(hs$coefficients),
      "kappa_mean" %in% names(hs$coefficients)
    )
  })

  run_test("print.signaly_horseshoe: runs without error", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    invisible(capture.output(print(hs)))
  })

  run_test("summary.signaly_horseshoe: runs without error", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    invisible(capture.output(summary(hs)))
  })

  # -- select_by_shrinkage --
  run_test("select_by_shrinkage: basic", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    sel <- select_by_shrinkage(hs, threshold = 0.5, verbose = FALSE)
    stopifnot(
      is.list(sel),
      is.character(sel$selected),
      sel$method == "shrinkage",
      is.data.frame(sel$details)
    )
  })

  run_test("select_by_shrinkage: custom threshold", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    sel_loose <- select_by_shrinkage(hs, threshold = 0.9, verbose = FALSE)
    sel_tight <- select_by_shrinkage(hs, threshold = 0.1, verbose = FALSE)
    stopifnot(length(sel_loose$selected) >= length(sel_tight$selected))
  })

  # -- select_by_credible_interval --
  run_test("select_by_credible_interval: basic", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    sel <- select_by_credible_interval(hs, prob = 0.95, verbose = FALSE)
    stopifnot(
      is.list(sel),
      is.character(sel$selected),
      sel$method == "credible_interval"
    )
  })

  run_test("select_by_credible_interval: different prob levels", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    s90 <- select_by_credible_interval(hs, prob = 0.90, verbose = FALSE)
    s99 <- select_by_credible_interval(hs, prob = 0.99, verbose = FALSE)
    stopifnot(length(s90$selected) >= length(s99$selected))
  })

  # -- select_by_magnitude --
  run_test("select_by_magnitude: basic", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    sel <- select_by_magnitude(hs, threshold = 0.1, verbose = FALSE)
    stopifnot(
      is.list(sel),
      is.character(sel$selected),
      sel$method == "magnitude"
    )
  })

  run_test("select_by_magnitude: threshold ordering", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    s_low  <- select_by_magnitude(hs, threshold = 0.01, verbose = FALSE)
    s_high <- select_by_magnitude(hs, threshold = 1.0,  verbose = FALSE)
    stopifnot(length(s_low$selected) >= length(s_high$selected))
  })

  # -- select_by_projection --
  run_test("select_by_projection: basic (no CV)", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    proj <- select_by_projection(hs, cv_k = 0, n_draws_max = 50,
                                  verbose = FALSE)
    stopifnot(
      is.list(proj),
      is.character(proj$selected),
      proj$method == "projection",
      is.data.frame(proj$elpd),
      is.character(proj$ranking),
      length(proj$ranking) == P
    )
  })

  run_test("select_by_projection: with CV", {
    hs <- fit_horseshoe(
      y = as.numeric(Y), X = X,
      chains = 2, iter_warmup = 200, iter_sampling = 200,
      seed = 42, verbose = FALSE
    )
    proj <- select_by_projection(hs, cv_k = 3, n_draws_max = 50,
                                  verbose = FALSE)
    stopifnot(
      proj$suggested_size >= 1,
      proj$suggested_size <= P,
      is.numeric(proj$elpd_full)
    )
  })

  # -- class validation on select_by_* --
  run_test("select_by_shrinkage: rejects wrong class", {
    ok <- tryCatch({
      select_by_shrinkage(list(a = 1))
      FALSE
    }, error = function(e) grepl("signaly_horseshoe", e$message))
    stopifnot(ok)
  })

  run_test("select_by_credible_interval: rejects wrong class", {
    ok <- tryCatch({
      select_by_credible_interval(42)
      FALSE
    }, error = function(e) grepl("signaly_horseshoe", e$message))
    stopifnot(ok)
  })

  run_test("select_by_magnitude: rejects wrong class", {
    ok <- tryCatch({
      select_by_magnitude("not a fit")
      FALSE
    }, error = function(e) grepl("signaly_horseshoe", e$message))
    stopifnot(ok)
  })

  run_test("select_by_projection: rejects wrong class", {
    ok <- tryCatch({
      select_by_projection(data.frame(x = 1))
      FALSE
    }, error = function(e) grepl("signaly_horseshoe", e$message))
    stopifnot(ok)
  })

} else {
  for (nm in c(
    "fit_horseshoe: basic fit",
    "fit_horseshoe: custom p0",
    "fit_horseshoe: coefficients structure",
    "print.signaly_horseshoe",
    "summary.signaly_horseshoe",
    "select_by_shrinkage: basic",
    "select_by_shrinkage: custom threshold",
    "select_by_credible_interval: basic",
    "select_by_credible_interval: different prob levels",
    "select_by_magnitude: basic",
    "select_by_magnitude: threshold ordering",
    "select_by_projection: basic (no CV)",
    "select_by_projection: with CV",
    "select_by_shrinkage: rejects wrong class",
    "select_by_credible_interval: rejects wrong class",
    "select_by_magnitude: rejects wrong class",
    "select_by_projection: rejects wrong class"
  )) {
    skip_test(nm, "(cmdstanr/CmdStan not available)")
  }
}

# ===========================================================================
# 11. signal_analysis (MASTER PIPELINE)
#     Strategy: run ONE full pipeline with hpgc, reuse the result for
#     print/summary. All other signal_analysis tests use lightweight
#     methods (wavelet, pca, dfm, unitroot) to avoid MCMC memory pressure.
# ===========================================================================
cat("\n=== 11. SIGNAL ANALYSIS (MASTER PIPELINE) ================================\n")
gc(verbose = FALSE)

# -- Single heavy run: ALL non-horseshoe methods including hpgc --
run_test("signal_analysis: full pipeline (no horseshoe)", {
  sa_full <<- signal_analysis(
    data = panel,
    y_formula = "Y",
    methods = c("wavelet", "emd", "hpgc", "pca", "dfm", "unitroot"),
    filter_config = list(hpgc_chains = 2, hpgc_iterations = 2000,
                         hpgc_burnin = 500),
    pca_config = list(n_boot = 30),
    verbose = FALSE
  )
  stopifnot(
    inherits(sa_full, "signal_analysis"),
    !is.null(sa_full$filters$wavelet),
    !is.null(sa_full$filters$emd),
    !is.null(sa_full$filters$hpgc),
    !is.null(sa_full$pca),
    !is.null(sa_full$dfm),
    !is.null(sa_full$unitroot),
    !is.null(sa_full$interpretation)
  )
})

run_test("print.signal_analysis: runs without error", {
  out <- capture.output(print(sa_full))
  stopifnot(any(grepl("SignalY", out)))
})

run_test("summary.signal_analysis: runs without error", {
  out <- capture.output(summary(sa_full))
  stopifnot(any(grepl("Detailed|Summary|OVERVIEW", out, ignore.case = TRUE)))
})

run_test("signal_analysis: interpretation structure", {
  interp <- sa_full$interpretation
  stopifnot(
    is.list(interp),
    !is.null(interp$signal_characteristics),
    !is.null(interp$factor_structure),
    !is.null(interp$persistence),
    is.character(interp$overall_summary),
    nchar(interp$overall_summary) > 0
  )
})

rm(sa_full); gc(verbose = FALSE)

# -- Lightweight tests (no hpgc MCMC) --
run_test("signal_analysis: formula interface", {
  res <- signal_analysis(
    data = panel,
    y_formula = Y ~ V1 + V2 + V3 + V4 + V5,
    methods = c("pca", "dfm", "unitroot"),
    pca_config = list(n_boot = 30),
    verbose = FALSE
  )
  stopifnot(inherits(res, "signal_analysis"), res$data$p == 5)
})

run_test("signal_analysis: subset of methods", {
  res <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("wavelet", "unitroot"),
    verbose = FALSE
  )
  stopifnot(
    !is.null(res$filters$wavelet),
    is.null(res$filters$emd),
    is.null(res$horseshoe),
    is.null(res$pca),
    !is.null(res$unitroot)
  )
})

run_test("signal_analysis: first_difference = TRUE", {
  res <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("wavelet", "pca"),
    first_difference = TRUE,
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )
  stopifnot(res$data$n == N - 1)
})

run_test("signal_analysis: standardize = FALSE", {
  res <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("pca", "unitroot"),
    standardize = FALSE,
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )
  stopifnot(inherits(res, "signal_analysis"))
})

run_test("signal_analysis: seed reproducibility", {
  r1 <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("pca"),
    pca_config = list(n_boot = 30),
    seed = 123, verbose = FALSE
  )
  r2 <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("pca"),
    pca_config = list(n_boot = 30),
    seed = 123, verbose = FALSE
  )
  stopifnot(
    identical(r1$pca$loadings, r2$pca$loadings),
    identical(r1$pca$eigenvalues, r2$pca$eigenvalues)
  )
})

run_test("signal_analysis: NA interpolation", {
  panel_na <- panel
  panel_na$V1[c(3, 7)] <- NA
  panel_na$Y[5] <- NA
  res <- signal_analysis(
    data = panel_na, y_formula = "Y",
    methods = c("pca", "unitroot"),
    na_action = "interpolate",
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )
  stopifnot(inherits(res, "signal_analysis"))
})

run_test("signal_analysis: NA omit", {
  panel_na <- panel
  panel_na$V1[c(3, 7)] <- NA
  res <- signal_analysis(
    data = panel_na, y_formula = "Y",
    methods = c("pca", "unitroot"),
    na_action = "omit",
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )
  stopifnot(res$data$n == N - 2)
})

# ===========================================================================
# 12. signal_analysis WITH HORSESHOE (if available)
# ===========================================================================
cat("\n=== 12. SIGNAL ANALYSIS + HORSESHOE ======================================\n")

gc(verbose = FALSE)

if (has_cmdstanr && has_cmdstan) {

  run_test("signal_analysis: horseshoe via pipeline", {
    res <- signal_analysis(
      data = panel, y_formula = "Y",
      methods = c("horseshoe", "pca", "unitroot"),
      horseshoe_config = list(chains = 2, iter_warmup = 200,
                              iter_sampling = 200),
      pca_config = list(n_boot = 20),
      verbose = FALSE
    )
    stopifnot(
      !is.null(res$horseshoe),
      !is.null(res$horseshoe$selection),
      length(res$horseshoe$selection$selected) >= 0
    )
    rm(res); gc(verbose = FALSE)
  })

  run_test("signal_analysis: horseshoe selection methods", {
    for (method in c("shrinkage", "credible_interval", "magnitude")) {
      res <- signal_analysis(
        data = panel, y_formula = "Y",
        methods = "horseshoe",
        horseshoe_config = list(
          chains = 2, iter_warmup = 200, iter_sampling = 200,
          selection_method = method
        ),
        verbose = FALSE
      )
      stopifnot(!is.null(res$horseshoe$selection))
      rm(res); gc(verbose = FALSE)
    }
  })

} else {
  skip_test("signal_analysis: horseshoe via pipeline",
            "(cmdstanr not available)")
  skip_test("signal_analysis: horseshoe selection methods",
            "(cmdstanr not available)")
}

# ===========================================================================
# 13. PLOT + IPLOT (plotly)
# ===========================================================================
cat("\n=== 13. PLOT / IPLOT (plotly) =============================================\n")

has_plotly <- requireNamespace("plotly", quietly = TRUE)

if (has_plotly) {

  # Build ONE result with lightweight methods for all plot tests
  sa_plot <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("wavelet", "emd", "pca", "dfm", "unitroot"),
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )

  run_test("plot.signal_analysis: returns plotly htmlwidget", {
    p <- plot(sa_plot)
    stopifnot(inherits(p, "plotly"), inherits(p, "htmlwidget"))
  })

  run_test("iplot: alias returns same class", {
    p <- iplot(sa_plot)
    stopifnot(inherits(p, "plotly"), inherits(p, "htmlwidget"))
  })

  run_test("plot.signal_analysis: filters panel present", {
    p <- plot(sa_plot)
    stopifnot(inherits(p, "plotly"))
  })

  rm(sa_plot); gc(verbose = FALSE)

  run_test("plot.signal_analysis: PCA-only result", {
    res <- signal_analysis(
      data = panel, y_formula = "Y",
      methods = c("pca"),
      pca_config = list(n_boot = 20),
      verbose = FALSE
    )
    p <- plot(res)
    stopifnot(inherits(p, "plotly"))
  })

  run_test("plot.signal_analysis: DFM-only result", {
    res <- signal_analysis(
      data = panel, y_formula = "Y",
      methods = c("dfm"),
      verbose = FALSE
    )
    p <- plot(res)
    stopifnot(inherits(p, "plotly"))
  })

  if (has_cmdstanr && has_cmdstan) {
    run_test("plot.signal_analysis: with horseshoe", {
      res <- signal_analysis(
        data = panel, y_formula = "Y",
        methods = c("horseshoe"),
        horseshoe_config = list(chains = 2, iter_warmup = 200,
                                iter_sampling = 200),
        verbose = FALSE
      )
      p <- plot(res)
      stopifnot(inherits(p, "plotly"))
      rm(res); gc(verbose = FALSE)
    })
  } else {
    skip_test("plot.signal_analysis: with horseshoe",
              "(cmdstanr not available)")
  }

} else {
  for (nm in c(
    "plot.signal_analysis: returns plotly htmlwidget",
    "iplot: alias returns same class",
    "plot.signal_analysis: with filters only",
    "plot.signal_analysis: with PCA only",
    "plot.signal_analysis: with DFM only",
    "plot.signal_analysis: with horseshoe"
  )) {
    skip_test(nm, "(plotly not installed)")
  }
}

# ===========================================================================
# 14. PLOT WITHOUT PLOTLY (graceful error)
# ===========================================================================
cat("\n=== 14. EDGE CASES =======================================================\n")

run_test("plot.signal_analysis: error msg when plotly missing", {
  res <- signal_analysis(
    data = panel, y_formula = "Y",
    methods = c("pca"),
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )
  if (!has_plotly) {
    ok <- tryCatch({
      plot(res)
      FALSE
    }, error = function(e) grepl("plotly", e$message, ignore.case = TRUE))
    stopifnot(ok)
  } else {
    TRUE
  }
})

run_test("signal_analysis: invalid y_formula -> error", {
  ok <- tryCatch({
    signal_analysis(data = panel, y_formula = "NONEXISTENT", verbose = FALSE)
    FALSE
  }, error = function(e) grepl("not found", e$message))
  stopifnot(ok)
})

run_test("signal_analysis: invalid method name -> error", {
  ok <- tryCatch({
    signal_analysis(data = panel, y_formula = "Y",
                    methods = "magic", verbose = FALSE)
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
})

run_test("signal_analysis: matrix input (not data.frame)", {
  mat <- as.matrix(panel)
  res <- signal_analysis(
    data = mat, y_formula = "Y",
    methods = c("pca", "unitroot"),
    pca_config = list(n_boot = 20),
    verbose = FALSE
  )
  stopifnot(inherits(res, "signal_analysis"))
})

run_test("filter_wavelet: via apply_to_columns", {
  X_small <- X[, 1:3]
  res <- apply_to_columns(X_small, filter_wavelet, extract = "trend",
                          J = 3, verbose = FALSE)
  stopifnot(is.data.frame(res), nrow(res) == N, ncol(res) == 3)
})

run_test("filter_emd: via apply_to_columns", {
  X_small <- X[, 1:3]
  res <- apply_to_columns(X_small, filter_emd, extract = "trend",
                          verbose = FALSE)
  stopifnot(is.data.frame(res), nrow(res) == N, ncol(res) == 3)
})

# ===========================================================================
# FINAL REPORT
# ===========================================================================
cat("\n")
cat("###########################################################################\n")
cat("#                          RESULTS SUMMARY                               #\n")
cat("###########################################################################\n\n")

total <- passed + failed + skipped
cat(sprintf("  Total tests:   %d\n", total))
cat(sprintf("  Passed:        %d\n", passed))
cat(sprintf("  Failed:        %d\n", failed))
cat(sprintf("  Skipped:       %d\n", skipped))
cat("\n")

if (failed > 0) {
  cat("  FAILURES:\n")
  for (e in errors) cat(e, "\n")
  cat("\n")
}

if (failed == 0) {
  cat("  *** ALL TESTS PASSED ***\n")
} else {
  cat(sprintf("  *** %d TEST(S) FAILED ***\n", failed))
}

cat("\n###########################################################################\n\n")

quit(status = failed, save = "no")
