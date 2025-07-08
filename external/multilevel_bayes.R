# windows
setwd("C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/")

# MacOS
setwd("~/projects/FrequencySliding/")

set.seed(123)

# Settings
n_subj <- 150
n_subscales <- 3
n_obs <- n_subj * n_subscales
n_predictors <- 2

# Simulate predictors
predictors <- matrix(rnorm(n_obs * n_predictors), ncol = n_predictors)
predictors <- scale(predictors, center = TRUE, scale = FALSE)  # Centering

# Simulate fixed intercepts & slopes (true values)
intercepts <- c(10, 15, 20)
beta <- matrix(c(2, -1, 1.5, 0.5, -2, 1), nrow = n_subscales, byrow = TRUE)

# Simulate random effects covariance structure
L_random <- chol(matrix(c(1, 0.7, 0.7,
                          0.7, 1, 0.7,
                          0.7, 0.7, 1), nrow = 3))  # Cholesky of random effects correlation

sigma_random <- c(2, 2, 2)  # SDs of random effects (intercept, slopes)

# Simulate random effects for subjects
z_random <- matrix(rnorm(n_subj * 3), nrow = n_subj, ncol = 3)
random_effects <- z_random %*% diag(sigma_random) %*% t(L_random)

# Subject & subscale indices
subj_id <- rep(1:n_subj, each = n_subscales)
subscale_id <- rep(1:n_subscales, times = n_subj)

# Generate outcome (BDI_score)
mu <- numeric(n_obs)
for (i in 1:n_obs) {
  subj <- subj_id[i]
  subscale <- subscale_id[i]
  pred <- predictors[i, ]
  
  mu[i] <- intercepts[subscale] + random_effects[subj, 1] + 
    (beta[subscale, 1] + random_effects[subj, 2]) * predictors[i, 1] + 
    (beta[subscale, 2] + random_effects[subj, 3]) * predictors[i, 2]
  
}

# Simulate with Student-t noise
sigma <- 2     # Scale
nu <- 5        # Degrees of freedom (heavier tails for small nu)
BDI_score <- mu + (sigma * rt(n_obs, df = nu))  # t-distributed noise

# Final data list for Stan
synthetic_data <- list(
  n_obs = n_obs,
  n_subj = n_subj,
  n_subscales = n_subscales,
  subj_id = subj_id,
  subscale_id = subscale_id,
  predictors = predictors,
  BDI_score = BDI_score
)

str(synthetic_data)  # Check structure

library(rstan)

stan_model <- stan_model(file = "./external/MultLevel_BDI.stan")

# Fit model
fit <- sampling(
  stan_model,
  data = synthetic_data,  # your list of data
  chains = 4,
  iter = 2000,
  seed = 123
)

library(bayesplot)

# Extract posterior samples from rstan fit
posterior <- as.array(fit)  # 3D array: iterations × chains × parameters

# Density plots (posterior distributions)
mcmc_areas(
  posterior,
  pars = c("Rho[2,1]", "Rho[3,1]"),
  prob = 0.95
)

