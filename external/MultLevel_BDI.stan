data {
  int<lower=1> n_obs;            // total observations (subject × subscale)
  int<lower=1> n_subj;           // number of subjects
  int<lower=1> n_subscales;      // number of subscales (probably 3)
  int<lower=1> subj_id[n_obs];   // subject ID for each observation
  int<lower=1> subscale_id[n_obs]; // subscale ID for each observation
  matrix[n_obs, 2] predictors;   // predictors: Two brain deviation variables
  vector[n_obs] BDI_score;       // observed BDI subscale score (long format)
}

parameters {
  // Fixed effects (shared across subscales)
  vector[n_subscales] intercept;
  matrix[n_subscales, 2] beta; // slopes per subscale

  // Random effects (subject-level deviations per subscale for intercept + slopes)
  vector[3] z_random[n_subj]; // 3 random effects per subject (1 intercept + 2 slopes)
  vector<lower=0>[3] sigma_random; // SDs of random effects
  cholesky_factor_corr[3] L_random; // Correlations among random effects. Mean rho is 3x3 matrix

  real<lower=0> sigma; // residual SD
  real<lower=2> nu;      // degrees of freedom for Student-t (must be >2 for finite variance)
}

transformed parameters {
  matrix[n_subj, 3] random_effects;
  for (i in 1:n_subj) {
    random_effects[i] = (diag_pre_multiply(sigma_random, L_random) * z_random[i])';
  }
}

model {
  vector[n_obs] mu;

  // Priors
  intercept ~ normal(0, 5);
  to_vector(beta) ~ normal(0, 2);
  sigma_random ~ exponential(1);
  L_random ~ lkj_corr_cholesky(2);
  sigma ~ exponential(1);
  nu ~ gamma(2, 0.1);  // favors moderate values, avoids heavy tails dominating

  // Random effect standard normal prior
  for (i in 1:n_subj) {
    z_random[i] ~ normal(0, 1);
  }

  // Observation model
  for (i in 1:n_obs) {
    int subj = subj_id[i];
    int subscale = subscale_id[i];

    mu[i] = intercept[subscale_id[i]] + random_effects[subj_id[i], 1] +
        (beta[subscale_id[i], 1] + random_effects[subj_id[i], 2]) * predictors[i, 1] +
        (beta[subscale_id[i], 2] + random_effects[subj_id[i], 3]) * predictors[i, 2];

  }

  BDI_score ~ student_t(nu, mu, sigma);
}

generated quantities{
  vector[n_obs] y_rep;
  matrix[3,3] Rho;
  Rho = multiply_lower_tri_self_transpose(L_random);
  
  for (i in 1:n_obs) {
    int subj = subj_id[i];
    int subscale = subscale_id[i];

    real mu_i = intercept[subscale] + random_effects[subj, 1] +
      (beta[subscale, 1] + random_effects[subj, 2]) * predictors[i, 1] +
      (beta[subscale, 2] + random_effects[subj, 3]) * predictors[i, 2];

    y_rep[i] = student_t_rng(nu, mu_i, sigma);
  }
}

