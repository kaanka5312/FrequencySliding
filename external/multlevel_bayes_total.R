set.seed(123)

# === PARAMETERS ===
n_subj <- 37
n_predictors <- 2

# === STEP 1: Generate predictors (centered) ===
predictors <- matrix(rnorm(n_subj * n_predictors), ncol = n_predictors)
predictors <- scale(predictors, center = TRUE, scale = FALSE)

# === STEP 2: True fixed effects ===
intercept <- 0.5              # from N(0, 0.5)
beta <- c(0.4, -0.3)         # from N(0, 0.25)

# === STEP 3: Random effects ===
sigma_random <- c(0.4, 0.25, 0.25)   # from Exp(2)
L_random <- chol(matrix(c(1, 0.3, 0.2,
                          0.3, 1, 0.4,
                          0.2, 0.4, 1), 
                        nrow = 3))   # LKJ(2)

z_random <- matrix(rnorm(n_subj * 3), nrow = n_subj, ncol = 3)
random_effects <- z_random %*% diag(sigma_random) %*% t(L_random)

# === STEP 4: Compute linear predictor ===
mu <- numeric(n_subj)
for (i in 1:n_subj) {
  mu[i] <- intercept + random_effects[i, 1] +
    (beta[1] + random_effects[i, 2]) * predictors[i, 1] +
    (beta[2] + random_effects[i, 3]) * predictors[i, 2]
}

# === STEP 5: Add noise ===
sigma <- 0.6   # from Exp(2)
nu <- 8        # from Gamma(8,1)
BDI_score <- mu + sigma * rt(n_subj, df = nu)

# === STEP 6: Return list for Stan model ===
synthetic_data <- list(
  n_subj = n_subj,
  predictors = predictors,
  BDI_score = BDI_score,
  prior_only = 0
)

library(rstan)

stan_model <- stan_model(file = "./external/MultLevel_BDI_total.stan")

# Fit model
fit_syn <- sampling(
  stan_model,
  data = synthetic_data,  # your list of data
  chains = 4,
  iter = 2000,
  seed = 123
)

posterior_syn <- as_draws_df(fit_syn)

library(bayesplot)
y_rep <- rstan::extract(fit_syn, pars = "y_rep")$y_rep


# Overlay observed and simulated densities (first 100 simulated sets)
ppc_dens_overlay(synthetic_data$BDI_score, y_rep[1:1000, ]) + xlim(c(-5,5))
print(fit_syn,pars=c("intercept","beta"))
mcmc_areas(
  posterior_syn,
  pars = c("intercept", "beta[1]", "beta[2]"),
  prob = 0.95  # 95% credible intervals
)

# DATA # 
bp_data <- subject_data_list[[5]] %>%
  filter(group == "BP") %>%
  mutate(
    subj_numeric = as.integer(factor(subj_id)),  # Numeric subject IDs for Stan
    subscale_numeric = as.integer(rep(1, times =37)),  # Numeric subscale IDs
    BDI_std = scale(BDI_total_score)[, 1]  # Z-score standardization
  )

# Prepare Stan Data List
stan_data_bp <- list(
  n_subj = length(unique(bp_data$subj_numeric)),
  predictors = as.matrix(bp_data[, c("zdev_uni_mean", "zdev_trans_mean")]),
  BDI_score = bp_data$BDI_std,
  prior_only = 0
)

# Inspect
str(stan_data_bp)

# Fit model
fit_data <- sampling(
  stan_model,
  data = stan_data_bp,  # your list of data
  chains = 4,
  iter = 2000,
  seed = 123
)

posterior_data <- as_draws_df(fit_data)
y_rep <- rstan::extract(fit_data, pars = "y_rep")$y_rep
ppc_dens_overlay(stan_data_bp$BDI_score, y_rep[1:100, ])
print(fit_data,pars=c("intercept","beta"))
mcmc_areas(
  posterior_data,
  pars = c("intercept", "beta[1]", "beta[2]"),
  prob = 0.95  # 95% credible intervals
)