#' @title SignalY: Signal Extraction from Panel Data via Bayesian Sparse
#'   Regression and Spectral Decomposition
#'
#' @description
#' SignalY provides a comprehensive methodological framework for extracting
#' latent signals from panel data through the integration of spectral
#' decomposition methods, Bayesian variable selection, and automated technical
#' interpretation. The package is designed for researchers working with
#' multivariate time series who seek to distinguish underlying structural
#' dynamics from phenomenological noise.
#'
#' @section Philosophical Foundation:
#' The package operationalizes a distinction between **latent structure** and

#' **phenomenological dynamics**. In complex systems, observed variables often
#' represent the superposition of: (1) underlying generative processes that
#' exhibit persistent, structured behavior; and (2) transient perturbations,
#' measurement noise, and stochastic fluctuations. SignalY provides tools to
#' decompose this mixture and identify which candidate variables contribute
#' meaningfully to the latent structure of a target signal.
#'
#' This framework recognizes that panel data exhibit **multivariate non-linear
#' interdependence**: the relationships between variables may be complex,
#' non-additive, and evolve over time. The methods implemented here are robust
#' to such complexities while remaining interpretable.
#'
#' @section Core Methodological Components:
#'
#' **1. Spectral Decomposition (Signal Filtering)**
#'
#' The package implements three complementary approaches to extract trend
#' components from time series:
#'
#' \itemize{
#'   \item **Wavelet Multiresolution Analysis**: Using the maximal overlap
#'     discrete wavelet transform (MODWT) with configurable Daubechies wavelets,
#'     the signal is decomposed into scale-specific components. Lower-frequency
#'     detail levels (e.g., D3, D4) capture structural dynamics while
#'     higher-frequency levels capture transient noise.
#'   \item **Empirical Mode Decomposition (EMD)**: A data-adaptive method that
#'     decomposes signals into intrinsic mode functions (IMFs) without
#'     requiring pre-specified basis functions. The residual component captures
#'     the underlying trend.
#'   \item **Grant-Chan Embedded Hodrick-Prescott Filter**: A Bayesian
#'     implementation embedding the HP filter within an unobserved components
#'     model, allowing for principled uncertainty quantification around the
#'     extracted trend via Markov Chain Monte Carlo sampling.
#' }
#'
#' **2. Bayesian Variable Selection (Horseshoe Regression)**
#'
#' When the target signal Y is constructed from or influenced by a set of
#' candidate variables X, identifying which candidates are structurally
#' relevant versus informationally redundant is crucial. The regularized
#' Horseshoe prior provides:
#'
#' \itemize{
#'   \item **Adaptive shrinkage**: Coefficients for irrelevant variables are
#'     strongly shrunk toward zero (high kappa), while relevant variables
#'     escape shrinkage (low kappa).
#'   \item **Uncertainty quantification**: Full posterior distributions over
#'     coefficients enable credible interval construction.
#'   \item **Automatic sparsity detection**: The effective number of non-zero
#'     coefficients (m_eff) is estimated as part of the model.
#' }
#'
#' **3. Dimensionality Reduction and Factor Analysis**
#'
#' For high-dimensional panels, the package provides:
#'
#' \itemize{
#'   \item **Principal Component Analysis (PCA)**: With bootstrap significance
#'     testing to identify which variables load significantly on each component.
#'   \item **Dynamic Factor Models (DFM)**: For extracting common factors that
#'     drive co-movement in the panel.
#'   \item **Entropy-based interpretation**: Shannon entropy of loadings
#'     distinguishes between diffuse systemic movement (high entropy) and
#'     concentrated structural signals (low entropy).
#' }
#'
#' **4. Unit Root and Stationarity Testing**
#'
#' Comprehensive suite of tests to characterize the persistence properties of
#' extracted signals:
#'
#' \itemize{
#'   \item Augmented Dickey-Fuller (ADF) tests with drift and trend options
#'   \item Elliott-Rothenberg-Stock (ERS) DF-GLS and P-tests
#'   \item Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests
#'   \item Phillips-Perron tests
#' }
#'
#' @section Interpretation Framework:
#'
#' SignalY generates automated technical interpretations based on:
#'
#' \itemize{
#'   \item **Signal smoothness**: Comparing variance of second differences
#'     between original and filtered series
#'   \item **Trend persistence**: Whether extracted trends are deterministic
#'     or stochastic based on unit root tests
#'   \item **Information topology**: Entropy and distributional fit of PCA
#'     loadings indicating structural concentration
#'   \item **Sparsity ratio**: Proportion of candidate variables shrunk to
#'     zero under Horseshoe regression
#'   \item **Regime detection**: Identification of structural breakpoints in
#'     mean or volatility
#' }
#'
#' @section Important Caveats:
#'
#' SignalY provides **methodology**, not **theory**. The statistical
#' identification of relevant variables does not establish causal or structural
#' relationships without supporting domain theory. Users must:
#'
#' \enumerate{
#'   \item Justify variable inclusion based on domain knowledge
#'   \item Interpret sparsity results in theoretical context
#'   \item Recognize that statistical significance is necessary but not
#'     sufficient for structural claims
#' }
#'
#' @references
#'
#' Daubechies, I. (1992). Ten Lectures on Wavelets. SIAM.
#'
#' Grant, A. L., & Chan, J. C. C. (2017). Reconciling output gaps: Unobserved
#' components model and Hodrick-Prescott filter. Journal of Economic Dynamics
#' and Control, 75, 114-121.
#'
#' Huang, N. E., et al. (1998). The empirical mode decomposition and the Hilbert
#' spectrum for nonlinear and non-stationary time series analysis. Proceedings
#' of the Royal Society A, 454(1971), 903-995.
#'
#' Percival, D. B., & Walden, A. T. (2000). Wavelet Methods for Time Series
#' Analysis. Cambridge University Press.
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of
#' Statistics, 11(2), 5018-5051.
#'
#' @author Jose Mauricio Gomez Julian \email{isadore.nabi@@pm.me}
#' @seealso
#' \itemize{
#'   \item \code{\link{signal_analysis}}: Master function for complete analysis
#'   \item \code{\link{filter_wavelet}}: Wavelet multiresolution analysis
#'   \item \code{\link{filter_emd}}: Empirical mode decomposition
#'   \item \code{\link{filter_hpgc}}: Grant-Chan HP filter
#'   \item \code{\link{fit_horseshoe}}: Regularized Horseshoe regression
#' }
#'
#' @docType package
#' @name SignalY-package
#' @aliases SignalY
#' @keywords package
#'
#' @importFrom stats complete.cases model.frame model.matrix model.response sd var
#' @importFrom utils head modifyList
"_PACKAGE"
