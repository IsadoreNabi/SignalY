// Regularized Horseshoe Prior for Sparse Regression
// Based on Piironen & Vehtari (2017): Sparsity information and regularization
// in the horseshoe and other shrinkage priors
//
// This model implements:
// - Non-centered parameterization for improved NUTS efficiency
// - Optional QR decomposition for multicollinearity handling
// - Slab regularization via Student-t distribution
// - Calibrated global shrinkage based on expected number of non-zeros

data {
  int<lower=1> N;                     // Number of observations
  int<lower=1> P;                     // Number of predictors
  vector[N] y;                        // Response variable
  matrix[N, P] X;                     // Predictor matrix (possibly QR-transformed)
  
  // Horseshoe hyperparameters
  real<lower=0> scale_global;         // Scale for global shrinkage (tau)
  real<lower=1> nu_global;            // Degrees of freedom for half-t prior on tau
  real<lower=1> nu_local;             // Degrees of freedom for half-t prior on lambdas
  real<lower=0> slab_scale;           // Slab scale for regularization (c)
  real<lower=1> slab_df;              // Slab degrees of freedom
}

transformed data {
  // Compute slab scale squared for efficiency
  real slab_scale2 = square(slab_scale);
}

parameters {
  // Intercept
  real alpha;
  
  // Non-centered parameterization for coefficients
  vector[P] z;                        // Standard normal draws
  
  // Global shrinkage parameter (half-t)
  real<lower=0> tau;                  // Global shrinkage scale
  
  // Local shrinkage parameters (half-t)
  vector<lower=0>[P] lambda;          // Local shrinkage scales
  
  // Slab regularization
  real<lower=0> caux;                 // Auxiliary for slab scale
  
  // Noise scale
  real<lower=0> sigma;
}

transformed parameters {
  // Regularized slab scale (regularized horseshoe)
  real<lower=0> c = slab_scale * sqrt(caux);
  
  // Regularized local scales (Eq. 2 in Piironen & Vehtari 2017)
  vector<lower=0>[P] lambda_tilde;
  
  // Actual regression coefficients
  vector[P] beta;
  
  // Linear predictor
  vector[N] mu;
  
  // Compute regularized local scales
  // lambda_tilde^2 = c^2 * lambda^2 / (c^2 + tau^2 * lambda^2)
  for (p in 1:P) {
    real lambda2 = square(lambda[p]);
    real tau2 = square(tau);
    lambda_tilde[p] = sqrt(c^2 * lambda2 / (c^2 + tau2 * lambda2));
  }
  
  // Non-centered parameterization: beta = tau * lambda_tilde .* z
  beta = tau * lambda_tilde .* z;
  
  // Linear predictor
  mu = alpha + X * beta;
}

model {
  // --- Priors ---
  
  // Intercept: weakly informative
  alpha ~ normal(0, 2);
  
  // Standard normal for non-centered parameterization
  z ~ std_normal();
  
  // Global shrinkage: half-t(nu_global, 0, scale_global)
  // Following the recommendation: scale_global = p0/(P-p0) * sigma/sqrt(N)
  tau ~ student_t(nu_global, 0, scale_global);
  
  // Local shrinkage: half-t(nu_local, 0, 1)
  // nu_local = 1 gives half-Cauchy (original horseshoe)
  // nu_local = 3-5 provides heavier regularization
  lambda ~ student_t(nu_local, 0, 1);
  
  // Slab scale auxiliary: inverse-gamma(slab_df/2, slab_df/2)
  // This gives c^2 ~ inverse-gamma(slab_df/2, slab_df * slab_scale^2 / 2)
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  
  // Noise scale: half-t(3, 0, 2)
  sigma ~ student_t(3, 0, 2);
  
  // --- Likelihood ---
  y ~ normal(mu, sigma);
}

generated quantities {
  // Shrinkage factors (kappa): how much each coefficient is shrunk toward zero
  // kappa = 1 / (1 + tau^2 * lambda^2 / c^2) when using regularized horseshoe
  // For standard horseshoe: kappa = 1 / (1 + tau^2 * lambda^2)
  vector<lower=0, upper=1>[P] kappa;
  
  // Effective number of non-zero coefficients
  real m_eff;
  
  // Log likelihood for LOO-CV
  vector[N] log_lik;
  
  // Posterior predictive samples
  vector[N] y_rep;
  
  // Compute shrinkage factors
  for (p in 1:P) {
    real lambda2 = square(lambda[p]);
    real tau2 = square(tau);
    // Regularized shrinkage factor
    kappa[p] = 1.0 / (1.0 + tau2 * lambda2 / (c^2));
  }
  
  // Effective number of non-zeros: sum of (1 - kappa)
  m_eff = P - sum(kappa);
  
  // Log likelihood contributions
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
  }
  
  // Posterior predictive
  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma);
  }
}
