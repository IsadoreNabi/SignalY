#' @title Unit Root and Stationarity Tests
#' @description
#' Comprehensive suite of unit root and stationarity tests for characterizing
#' the persistence properties of time series and extracted signals.
#'
#' @name unit_root
NULL


#' Comprehensive Unit Root Test Suite
#'
#' @description
#' Applies multiple unit root and stationarity tests to a time series,
#' providing an integrated assessment of persistence properties. Implements
#' Augmented Dickey-Fuller (ADF), Elliott-Rothenberg-Stock (ERS),
#' Kwiatkowski-Phillips-Schmidt-Shin (KPSS), and Phillips-Perron tests.
#'
#' @param y Numeric vector of the time series to test.
#' @param max_lags Maximum number of lags for ADF-type tests. If NULL,
#'   defaults to \code{floor(12 * (length(y)/100)^0.25)}.
#' @param significance_level Significance level for hypothesis testing.
#'   Default is 0.05.
#' @param verbose Logical indicating whether to print detailed results.
#'
#' @return A list of class "signaly_unitroot" containing:
#' \describe{
#'   \item{adf}{Results from ADF tests (none, drift, trend specifications)}
#'   \item{ers}{Results from ERS tests (DF-GLS and P-test)}
#'   \item{kpss}{Results from KPSS tests (level and trend)}
#'   \item{pp}{Results from Phillips-Perron tests}
#'   \item{summary}{Data frame summarizing all test results}
#'   \item{conclusion}{Integrated conclusion about stationarity}
#'   \item{persistence_type}{Classification: stationary, trend-stationary,
#'     difference-stationary, or inconclusive}
#' }
#'
#' @details
#' The battery of tests addresses different null hypotheses and specifications:
#'
#' **Augmented Dickey-Fuller (ADF)** tests the null of a unit root against
#' the alternative of stationarity. Three specifications are tested:
#' \itemize{
#'   \item **none**: No constant, no trend (random walk)
#'   \item **drift**: Constant included (random walk with drift)
#'   \item **trend**: Constant and linear trend
#' }
#'
#' **Elliott-Rothenberg-Stock (ERS)** tests provide more power than ADF by
#' using GLS detrending. Two variants:
#' \itemize{
#'   \item **DF-GLS**: GLS-detrended Dickey-Fuller test
#'   \item **P-test**: Point-optimal test
#' }
#'
#' **KPSS** reverses the hypotheses: null is stationarity, alternative is
#' unit root. This allows testing the stationarity hypothesis directly.
#'
#' **Phillips-Perron** uses non-parametric corrections for serial correlation,
#' avoiding lag selection issues.
#'
#' @section Interpretation Strategy:
#' The function synthesizes results using the following logic:
#'
#' \enumerate{
#'   \item If ADF/ERS reject unit root AND KPSS fails to reject stationarity:
#'     Series is likely **stationary**
#'   \item If ADF/ERS fail to reject AND KPSS rejects stationarity:
#'     Series likely has **unit root** (difference-stationary)
#'   \item If only trend-ADF rejects: Series is likely **trend-stationary**
#'   \item Conflicting results indicate **inconclusive** or structural breaks
#' }
#'
#' @references
#' Dickey, D. A., & Fuller, W. A. (1979). Distribution of the Estimators for
#' Autoregressive Time Series with a Unit Root. Journal of the American
#' Statistical Association, 74(366), 427-431.
#'
#' Elliott, G., Rothenberg, T. J., & Stock, J. H. (1996). Efficient Tests for
#' an Autoregressive Unit Root. Econometrica, 64(4), 813-836.
#'
#' Kwiatkowski, D., Phillips, P. C. B., Schmidt, P., & Shin, Y. (1992).
#' Testing the null hypothesis of stationarity against the alternative of a
#' unit root. Journal of Econometrics, 54(1-3), 159-178.
#'
#' Phillips, P. C. B., & Perron, P. (1988). Testing for a unit root in time
#' series regression. Biometrika, 75(2), 335-346.
#'
#' @examples
#' set.seed(123)
#' stationary <- arima.sim(list(ar = 0.5), n = 100)
#' result <- test_unit_root(stationary)
#' print(result$conclusion)
#'
#' nonstationary <- cumsum(rnorm(100))
#' result2 <- test_unit_root(nonstationary)
#' print(result2$conclusion)
#'
#' @seealso \code{\link[urca]{ur.df}}, \code{\link[urca]{ur.ers}},
#'   \code{\link[urca]{ur.kpss}}, \code{\link[urca]{ur.pp}}
#'
#' @export
test_unit_root <- function(y,
                           max_lags = NULL,
                           significance_level = 0.05,
                           verbose = FALSE) {

  if (!requireNamespace("urca", quietly = TRUE)) {
    stop("Package 'urca' is required. Please install it.", call. = FALSE)
  }

  y <- as.numeric(y)
  n <- length(y)

  if (any(is.na(y))) {
    stop("Series contains NA values. Please handle missing data first.",
         call. = FALSE)
  }

  if (n < 20) {
    stop("Series too short for reliable unit root testing (minimum 20 obs).",
         call. = FALSE)
  }

  if (is.null(max_lags)) {
    max_lags <- floor(12 * (n / 100)^0.25)
  }

  level <- switch(
    as.character(significance_level),
    "0.01" = "1pct",
    "0.05" = "5pct",
    "0.1" = "10pct",
    "5pct"
  )

  results <- list()

  adf_none <- urca::ur.df(y, type = "none", lags = max_lags, selectlags = "AIC")
  adf_drift <- urca::ur.df(y, type = "drift", lags = max_lags, selectlags = "AIC")
  adf_trend <- urca::ur.df(y, type = "trend", lags = max_lags, selectlags = "AIC")

  results$adf <- list(
    none = interpret_adf(adf_none, level),
    drift = interpret_adf(adf_drift, level),
    trend = interpret_adf(adf_trend, level)
  )

  ers_dfgls_const <- urca::ur.ers(y, type = "DF-GLS", model = "constant",
                                  lag.max = max_lags)
  ers_dfgls_trend <- urca::ur.ers(y, type = "DF-GLS", model = "trend",
                                  lag.max = max_lags)
  ers_p_const <- urca::ur.ers(y, type = "P-test", model = "constant",
                              lag.max = max_lags)
  ers_p_trend <- urca::ur.ers(y, type = "P-test", model = "trend",
                              lag.max = max_lags)

  results$ers <- list(
    dfgls_constant = interpret_ers(ers_dfgls_const, level),
    dfgls_trend = interpret_ers(ers_dfgls_trend, level),
    ptest_constant = interpret_ers_ptest(ers_p_const, level),
    ptest_trend = interpret_ers_ptest(ers_p_trend, level)
  )

  kpss_mu_short <- urca::ur.kpss(y, type = "mu", lags = "short")
  kpss_mu_long <- urca::ur.kpss(y, type = "mu", lags = "long")
  kpss_tau_short <- urca::ur.kpss(y, type = "tau", lags = "short")
  kpss_tau_long <- urca::ur.kpss(y, type = "tau", lags = "long")

  results$kpss <- list(
    level_short = interpret_kpss(kpss_mu_short, level),
    level_long = interpret_kpss(kpss_mu_long, level),
    trend_short = interpret_kpss(kpss_tau_short, level),
    trend_long = interpret_kpss(kpss_tau_long, level)
  )

  pp_alpha_short <- urca::ur.pp(y, type = "Z-alpha", model = "trend",
                                lags = "short")
  pp_tau_short <- urca::ur.pp(y, type = "Z-tau", model = "trend",
                              lags = "short")
  pp_alpha_long <- urca::ur.pp(y, type = "Z-alpha", model = "trend",
                               lags = "long")
  pp_tau_long <- urca::ur.pp(y, type = "Z-tau", model = "trend",
                             lags = "long")

  results$pp <- list(
    z_alpha_short = interpret_pp(pp_alpha_short, level),
    z_tau_short = interpret_pp(pp_tau_short, level),
    z_alpha_long = interpret_pp(pp_alpha_long, level),
    z_tau_long = interpret_pp(pp_tau_long, level)
  )

  summary_df <- create_unitroot_summary(results, level)
  results$summary <- summary_df

  conclusion <- synthesize_unitroot_results(results, significance_level)
  results$conclusion <- conclusion$text
  results$persistence_type <- conclusion$type

  class(results) <- c("signaly_unitroot", "list")

  if (verbose) {
    print_unitroot_results(results)
  }

  results
}


#' Interpret ADF Test Results
#'
#' @keywords internal
interpret_adf <- function(adf_obj, level) {
  model_type <- adf_obj@model

  if (model_type == "none") {
    stat_name <- "tau1"
    crit <- adf_obj@cval[stat_name, level]
    stat <- adf_obj@teststat["statistic", stat_name]
    reject <- stat < crit
    interpretation <- ifelse(reject, "Unit root rejected (stationary)",
                             "Unit root not rejected")
  } else if (model_type == "drift") {
    tau2_crit <- adf_obj@cval["tau2", level]
    tau2_stat <- adf_obj@teststat["statistic", "tau2"]
    phi1_crit <- adf_obj@cval["phi1", level]
    phi1_stat <- adf_obj@teststat["statistic", "phi1"]

    tau2_reject <- tau2_stat < tau2_crit
    phi1_reject <- phi1_stat > phi1_crit

    stat <- tau2_stat
    crit <- tau2_crit
    reject <- tau2_reject

    if (tau2_reject) {
      interpretation <- "Unit root rejected (stationary)"
    } else {
      if (phi1_reject) {
        interpretation <- "Unit root not rejected; drift present"
      } else {
        interpretation <- "Unit root not rejected; no drift"
      }
    }
  } else {
    tau3_crit <- adf_obj@cval["tau3", level]
    tau3_stat <- adf_obj@teststat["statistic", "tau3"]
    phi2_crit <- adf_obj@cval["phi2", level]
    phi2_stat <- adf_obj@teststat["statistic", "phi2"]
    phi3_crit <- adf_obj@cval["phi3", level]
    phi3_stat <- adf_obj@teststat["statistic", "phi3"]

    tau3_reject <- tau3_stat < tau3_crit
    phi2_reject <- phi2_stat > phi2_crit
    phi3_reject <- phi3_stat > phi3_crit

    stat <- tau3_stat
    crit <- tau3_crit
    reject <- tau3_reject

    if (tau3_reject) {
      interpretation <- "Unit root rejected (trend-stationary)"
    } else {
      if (phi3_reject) {
        interpretation <- "Unit root not rejected; trend may be present"
      } else {
        interpretation <- "Unit root not rejected; no trend"
      }
    }
  }

  list(
    test = paste("ADF", model_type),
    statistic = stat,
    critical_value = crit,
    reject_null = reject,
    interpretation = interpretation
  )
}


#' Interpret ERS DF-GLS Results
#'
#' @keywords internal
interpret_ers <- function(ers_obj, level) {
  stat <- ers_obj@teststat
  crit <- ers_obj@cval[, level]

  if (length(crit) > 1) crit <- crit[length(crit)]

  reject <- stat < crit

  list(
    test = paste("ERS DF-GLS", ers_obj@model),
    statistic = as.numeric(stat),
    critical_value = as.numeric(crit),
    reject_null = reject,
    interpretation = ifelse(reject, "Unit root rejected", "Unit root not rejected")
  )
}


#' Interpret ERS P-test Results
#'
#' @keywords internal
interpret_ers_ptest <- function(ers_obj, level) {
  stat <- ers_obj@teststat
  crit <- ers_obj@cval[, level]

  if (length(crit) > 1) crit <- crit[1]

  reject <- stat < crit

  list(
    test = paste("ERS P-test", ers_obj@model),
    statistic = as.numeric(stat),
    critical_value = as.numeric(crit),
    reject_null = reject,
    interpretation = ifelse(reject, "Unit root rejected", "Unit root not rejected")
  )
}


#' Interpret KPSS Results
#'
#' @keywords internal
interpret_kpss <- function(kpss_obj, level) {
  stat <- kpss_obj@teststat
  crit <- kpss_obj@cval[, level]

  if (length(crit) > 1) crit <- crit[1]

  reject <- stat > crit

  list(
    test = paste("KPSS", kpss_obj@type),
    statistic = as.numeric(stat),
    critical_value = as.numeric(crit),
    reject_null = reject,
    interpretation = ifelse(reject, "Stationarity rejected (unit root)",
                            "Stationarity not rejected")
  )
}


#' Interpret Phillips-Perron Results
#'
#' @keywords internal
interpret_pp <- function(pp_obj, level) {
  stat <- pp_obj@teststat[1]
  cval <- pp_obj@cval
  crit <- tryCatch({
    if (!is.null(dim(cval)) && !is.null(dimnames(cval))) {
      cval[1, level]
    } else if (!is.null(names(cval))) {
      cval[level]
    } else {
      cval[1]
    }
  }, error = function(e) NA)

  if (length(crit) > 1) crit <- crit[1]

  reject <- stat < crit

  list(
    test = paste("PP", pp_obj@type),
    statistic = as.numeric(stat),
    critical_value = as.numeric(crit),
    reject_null = reject,
    interpretation = ifelse(reject, "Unit root rejected", "Unit root not rejected")
  )
}


#' Create Summary Data Frame
#'
#' @keywords internal
create_unitroot_summary <- function(results, level) {
  rows <- list()

  for (test_type in c("adf", "ers", "kpss", "pp")) {
    for (test_name in names(results[[test_type]])) {
      test <- results[[test_type]][[test_name]]
      rows[[length(rows) + 1]] <- data.frame(
        category = toupper(test_type),
        test = test$test,
        statistic = test$statistic,
        critical_value = test$critical_value,
        reject_null = test$reject_null,
        interpretation = test$interpretation,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}


#' Synthesize Unit Root Results
#'
#' @keywords internal
synthesize_unitroot_results <- function(results, alpha) {
  adf_trend_reject <- results$adf$trend$reject_null
  adf_drift_reject <- results$adf$drift$reject_null
  ers_reject <- results$ers$dfgls_trend$reject_null || results$ers$ptest_trend$reject_null
  kpss_reject <- results$kpss$level_long$reject_null || results$kpss$trend_long$reject_null
  pp_reject <- results$pp$z_tau_long$reject_null

  unit_root_tests_reject <- sum(c(adf_trend_reject, adf_drift_reject, ers_reject, pp_reject))
  total_ur_tests <- 4

  if (unit_root_tests_reject >= 3 && !kpss_reject) {
    type <- "stationary"
    text <- paste(
      "STATIONARY: Strong evidence against unit root.",
      sprintf("%d/%d unit root tests reject the null.", unit_root_tests_reject, total_ur_tests),
      "KPSS fails to reject stationarity.",
      "The series appears to be stationary or trend-stationary."
    )
  } else if (unit_root_tests_reject <= 1 && kpss_reject) {
    type <- "difference_stationary"
    text <- paste(
      "UNIT ROOT PRESENT: Strong evidence for non-stationarity.",
      sprintf("Only %d/%d unit root tests reject.", unit_root_tests_reject, total_ur_tests),
      "KPSS rejects stationarity.",
      "The series likely has a unit root (difference-stationary)."
    )
  } else if (adf_trend_reject && !adf_drift_reject && !kpss_reject) {
    type <- "trend_stationary"
    text <- paste(
      "TREND-STATIONARY: Evidence for deterministic trend.",
      "ADF with trend rejects, but ADF with drift does not.",
      "The series may be stationary around a deterministic trend."
    )
  } else {
    type <- "inconclusive"
    text <- paste(
      "INCONCLUSIVE: Mixed evidence from different tests.",
      sprintf("Unit root tests rejecting: %d/%d", unit_root_tests_reject, total_ur_tests),
      sprintf("KPSS rejects stationarity: %s", ifelse(kpss_reject, "Yes", "No")),
      "Consider structural breaks or near-unit-root behavior."
    )
  }

  list(type = type, text = text)
}


#' Print Unit Root Results
#'
#' @keywords internal
print_unitroot_results <- function(x) {
  cat("\n")
  cat("================================================================\n")
  cat("UNIT ROOT AND STATIONARITY TEST RESULTS\n")
  cat("================================================================\n\n")

  cat("--- AUGMENTED DICKEY-FULLER TESTS ---\n")
  for (test in x$adf) {
    cat(sprintf("  %s: stat=%.3f, crit=%.3f, %s\n",
                test$test, test$statistic, test$critical_value,
                ifelse(test$reject_null, "REJECT", "FAIL TO REJECT")))
  }

  cat("\n--- ELLIOTT-ROTHENBERG-STOCK TESTS ---\n")
  for (test in x$ers) {
    cat(sprintf("  %s: stat=%.3f, crit=%.3f, %s\n",
                test$test, test$statistic, test$critical_value,
                ifelse(test$reject_null, "REJECT", "FAIL TO REJECT")))
  }

  cat("\n--- KPSS TESTS (null = stationarity) ---\n")
  for (test in x$kpss) {
    cat(sprintf("  %s: stat=%.3f, crit=%.3f, %s\n",
                test$test, test$statistic, test$critical_value,
                ifelse(test$reject_null, "REJECT", "FAIL TO REJECT")))
  }

  cat("\n--- PHILLIPS-PERRON TESTS ---\n")
  for (test in x$pp) {
    cat(sprintf("  %s: stat=%.3f, crit=%.3f, %s\n",
                test$test, test$statistic, test$critical_value,
                ifelse(test$reject_null, "REJECT", "FAIL TO REJECT")))
  }

  cat("\n================================================================\n")
  cat("CONCLUSION\n")
  cat("================================================================\n")
  cat(strwrap(x$conclusion, width = 65), sep = "\n")
  cat(sprintf("\nPersistence type: %s\n", toupper(x$persistence_type)))
  cat("================================================================\n\n")

  invisible(x)
}


#' @export
print.signaly_unitroot <- function(x, ...) {
  print_unitroot_results(x)
}

#' @export
print.unitroot_battery <- function(x, ...) {
  cat("Unit Root Test Battery\n")
  cat(sprintf("  Conclusion: %s\n", x$synthesis$conclusion))
  invisible(x)
}