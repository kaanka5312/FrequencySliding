# windows
setwd("C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/")

# MacOS
setwd("~/projects/FrequencySliding/")

set.seed(123)

# Constants
n_subj <- 30
n_subscales <- 3
n_predictors <- 2
n_obs <- n_subj * n_subscales

# Predictors: centered
predictors <- matrix(rnorm(n_obs * n_predictors), ncol = n_predictors)
predictors <- scale(predictors, center = TRUE, scale = FALSE)

# === STEP 1: Fix parameters closer to center of priors ===
intercepts <- c(0.5, 0, -0.5)  # Well within normal(0, 0.5)
beta <- matrix(c(0.3, -0.3, 0.2, 0.2, -0.2, 0.2), nrow = n_subscales, byrow = TRUE)  # normal(0, 0.5)

# === STEP 2: Simulate random effects ===
sigma_random <- c(0.5, 0.3, 0.3)  # Median of exponential(1) is ~0.69
L_random <- chol(matrix(c(1, 0.3, 0.4,
                          0.3, 1, 0.4,
                          0.4, 0.4, 1), nrow = 3))  # Reduced correlation

z_random <- matrix(rnorm(n_subj * 3), nrow = n_subj, ncol = 3)
random_effects <- z_random %*% diag(sigma_random) %*% t(L_random)

# === STEP 3: Subject & subscale indices ===
subj_id <- rep(1:n_subj, each = n_subscales)
subscale_id <- rep(1:n_subscales, times = n_subj)

# === STEP 4: Generate mu ===
mu <- numeric(n_obs)
for (i in 1:n_obs) {
  subj <- subj_id[i]
  subscale <- subscale_id[i]
  pred <- predictors[i, ]
  
  mu[i] <- intercepts[subscale] + random_effects[subj, 1] +
    (beta[subscale, 1] + random_effects[subj, 2]) * pred[1] +
    (beta[subscale, 2] + random_effects[subj, 3]) * pred[2]
}

# === STEP 5: Add noise ===
sigma <- 0.7  # Expected from exponential(1)
nu <- 8       # Gamma(4,1) center; avoid heavy tails

BDI_score <- mu + sigma * rt(n_obs, df = nu)

# === Check output ===
summary(BDI_score)
plot(density(BDI_score), col = "darkblue", lwd = 2, main = "BDI_score Prior Predictive")
abline(v = 0, col = "red")

# ==== STEP 6: Final Stan data ====
synthetic_data <- list(
  n_obs = n_obs,
  n_subj = n_subj,
  n_subscales = n_subscales,
  subj_id = subj_id,
  subscale_id = subscale_id,
  predictors = predictors,
  BDI_score = BDI_score,
  prior_only = 0
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

posterior <- rstan::extract(fit)
posterior <- rstan::extract(fit_data)
library(bayesplot)
y_rep <- rstan::extract(fit, pars = "y_rep")$y_rep
y_rep <- rstan::extract(fit_data, pars = "y_rep")$y_rep

# Overlay observed and simulated densities (first 100 simulated sets)
ppc_dens_overlay(synthetic_data$BDI_score, y_rep[1:100, ]) + xlim(c(-5,5))
ppc_dens_overlay(stan_data_bp$BDI_score, y_rep[1:100, ])

library(dplyr)
library(ggplot2)

# Assumes:
# - y is the observed BDI scores (length 450)
# - y_rep is a matrix of shape (n_draws x 450)
# - synthetic_data$subj_id gives subject IDs

y <- synthetic_data$BDI_score
y <- stan_data_bp$BDI_score

posterior_mean <- matrixStats::colMeans2(y_rep)
posterior_lower <- apply(y_rep, 2, quantile, 0.025)
posterior_upper <- apply(y_rep, 2, quantile, 0.975)

# Combine into a data frame
pp_check_df <- data.frame(
  subj_id = stan_data_bp$subj_id,
  subscale_id = stan_data_bp$subscale_id,
  obs = y,
  pred = posterior_mean,
  lower = posterior_lower,
  upper = posterior_upper
)

ggplot(pp_check_df, aes(x = factor(subscale_id), y = obs)) +
  geom_point(color = "black", size = 1.5) +  # observed
  geom_point(aes(y = pred), color = "blue", shape = 1) +  # predicted mean
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "blue") +  # 95% CI
  facet_wrap(~ subj_id, scales = "free_y", ncol = 10) +  # 150 subjects = 15x10 grid
  labs(
    x = "Subscale ID",
    y = "BDI score",
    title = "Posterior Predictive Checks by Subject"
  ) +
  theme_minimal(base_size = 9)

## Parameters
# Intercepts
intercept_mean <- apply(posterior$intercept, 2, mean)
intercept_lower <- apply(posterior$intercept, 2, quantile, 0.025)
intercept_upper <- apply(posterior$intercept, 2, quantile, 0.975)

# Betas
beta_mean <- apply(posterior$beta, c(2,3), mean)
beta_lower <- apply(posterior$beta, c(2,3), quantile, 0.025)
beta_upper <- apply(posterior$beta, c(2,3), quantile, 0.975)

# True values
true_intercepts <- intercepts  # Since BDI_score was standardized
true_beta <- beta  # Adjust beta similarly

intercept_df <- data.frame(
  type = "intercept",
  index = paste0("intercept[", 1:3, "]"),
  mean = intercept_mean,
  lower = intercept_lower,
  upper = intercept_upper,
  true = true_intercepts
)

# Beta data frame
beta_df <- data.frame(
  type = "beta",
  subscale = rep(1:3, each = 2),
  predictor = rep(1:2, times = 3),
  mean = as.vector(beta_mean),
  lower = as.vector(beta_lower),
  upper = as.vector(beta_upper),
  true = as.vector(true_beta)
)

library(ggplot2)

# Combine all
plot_df <- dplyr::bind_rows(
  intercept_df %>% rename(index = index),
  beta_df %>% mutate(index = paste0("beta[", subscale, ",", predictor, "]"))
)

ggplot(plot_df, aes(x = index, y = mean)) +
  geom_point(color = "blue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, color = "blue") +
  geom_point(aes(y = true), color = "red", shape = 4, size = 3) +
  coord_flip() +
  labs(
    title = "Posterior Estimates vs. Ground Truth",
    x = "Parameter",
    y = "Estimate (Standardized)"
  ) +
  theme_minimal()

### For data ###
mcmc_areas(
  posterior,
  pars = c("intercept[1]"),
  prob = 0.95  # 95% credible intervals
) +
  ggtitle("Posterior distributions of loadings and beta")