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
