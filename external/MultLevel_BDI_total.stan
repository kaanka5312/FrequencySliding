data {
  int<lower=1> n_subj;              // number of subjects
  matrix[n_subj, 2] predictors;     // two predictors per subject (e.g., brain deviations)
  vector[n_subj] BDI_score;         // total BDI score per subject (standardized if needed)
  int<lower=0, upper=1> prior_only; // 1 = prior predictive only
}

parameters {
  // Fixed effects
  real intercept;
  vector[2] beta;

  // Random effects (subject-specific intercept and slopes)
  vector[3] z_random[n_subj];              // non-centered parameterization
  vector<lower=0>[3] sigma_random;         // SDs of random intercept and slopes
  cholesky_factor_corr[3] L_random;        // correlation structure

  // Residual
  real<lower=0> sigma;
  real<lower=2> nu;  // degrees of freedom for Student-t
}

transformed parameters {
  matrix[n_subj, 3] random_effects;
  for (i in 1:n_subj) {
    random_effects[i] = (diag_pre_multiply(sigma_random, L_random) * z_random[i])';
  }
}

model {
  vector[n_subj] mu;

  // Priors
  intercept ~ normal(0, 0.5);
  beta ~ normal(0, 0.25);
  sigma_random ~ exponential(2);
  L_random ~ lkj_corr_cholesky(2);
  sigma ~ exponential(2);
  nu ~ gamma(8, 1);

  for (i in 1:n_subj) {
    z_random[i] ~ normal(0, 1);
  }

  for (i in 1:n_subj) {
    mu[i] = intercept + random_effects[i, 1] +
      (beta[1] + random_effects[i, 2]) * predictors[i, 1] +
      (beta[2] + random_effects[i, 3]) * predictors[i, 2];
  }

  if (prior_only == 0) {
    BDI_score ~ student_t(nu, mu, sigma);
  }
}

generated quantities {
  vector[n_subj] y_rep;
  matrix[3, 3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_random);

  for (i in 1:n_subj) {
    real mu_i = intercept + random_effects[i, 1] +
      (beta[1] + random_effects[i, 2]) * predictors[i, 1] +
      (beta[2] + random_effects[i, 3]) * predictors[i, 2];

    y_rep[i] = student_t_rng(nu, mu_i, sigma);
  }
}