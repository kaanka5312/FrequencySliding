library(brms)
PF_long <- PF_long_list[[5]] |>
  bind_rows() 

PF_subject <- PF_long %>%                                   # long table (all recordings)
  group_by(subj_id, group, sex, age, BDI, region_type) %>%  # keep covariates
  summarise(PF = mean(value), .groups = "drop") %>%         # average PF within subj × region-type
  pivot_wider(names_from = region_type, values_from = PF)   # → columns: unimodal, transmodal

## ──────────────────────────────────────────────────────────────────────────
## 2.  Scale continuous variables & code factors
## ──────────────────────────────────────────────────────────────────────────
dat <- PF_subject %>% 
  mutate(
    age_z        = as.numeric(scale(age)),
    uni_self_z   = as.numeric(scale(uni_self)),
    trans_self_z = as.numeric(scale(trans_self)),
    uni_nonself_z   = as.numeric(scale(uni_nonself)),
    trans_nonself_z = as.numeric(scale(trans_nonself)),
    BDI_z        = as.numeric(scale(BDI)),
    group        = factor(group),      # “BP”, “HC”
    sex          = factor(sex)         # e.g. “1”, “2”
  )

## ---- 1. formula -------------------------------------------------------
f_BDI <- bf(
  BDI_z ~                           # mean model
    uni_self_z + uni_nonself_z +
    trans_self_z + trans_nonself_z +
    group + sex + age_z +
    (uni_self_z + uni_nonself_z +
       trans_self_z + trans_nonself_z):group,   # interactions w/ group
  
  sigma ~ 0 + group,                # σBP , σHC  (log-link)
  nu    ~ 0 + group                 # νBP , νHC  (log-link, Student-t)
)

## ---- 2. one compact priors object ------------------------------------
priors <- c(
  # intercept for the mean equation
  prior(normal(0, 0.5),                class = "Intercept"),
  
  # ALL mean-model slopes (main & interactions)
  prior(normal(0, 0.30),               class = "b"),
  
  # BOTH log-σ coefficients (b_sigma_groupBP & b_sigma_groupHC)
  prior(normal(0, 0.30),               class = "b", dpar = "sigma"),
  
  # BOTH log-ν coefficients (b_nu_groupBP & b_nu_groupHC)
  prior(normal(log(10), 0.50),         class = "b", dpar = "nu")
)

## ---- 3. fit -----------------------------------------------------------
fit_grp <- brm(
  formula = f_BDI,
  data    = dat,          # your PF_subject data frame
  family  = student(),    # enables σ & ν sub-models
  prior   = priors,
  chains  = 4, iter = 4000, cores = 4,
  seed    = 123
)

pp_check(fit_grp, ndraws=1e2)
pp_check(fit_grp, type = "error_scatter_avg")

fit_grp_prior <- update(fit_grp, sample_prior = "only")  # runs the full model
pp_check(fit_grp_prior, ndraws=1e2)
#################
###############################################################################
##  STEP-BY-STEP PIPELINE:  collapse to UNIMODAL / TRANSMODAL  ▸  fit brms   ##
###############################################################################
library(dplyr)      # data wrangling
library(tidyr)      # pivot_wider()
library(forcats)    # fct_collapse()
library(brms)       # Bayesian model

## ──────────────────────────────────────────────────────────────────────────
## 1.  Collapse the four region labels → two levels, then average per subject
## ──────────────────────────────────────────────────────────────────────────
PF_long <- PF_long_list[[5]] |>
  bind_rows() 

PF_subject <- PF_long %>%                                   # long table (all recordings)
  mutate(
    region_type = fct_collapse(                             # uni_* → “unimodal”, trans_* → “transmodal”
      region_type,
      unimodal   = c("uni_self",  "uni_nonself"),
      transmodal = c("trans_self","trans_nonself")
    )
  ) %>% 
  group_by(subj_id, group, sex, age, BDI, region_type) %>%  # keep covariates
  summarise(PF = mean(value), .groups = "drop") %>%         # average PF within subj × region-type
  pivot_wider(names_from = region_type, values_from = PF)   # → columns: unimodal, transmodal

## ──────────────────────────────────────────────────────────────────────────
## 2.  Scale continuous variables & code factors
## ──────────────────────────────────────────────────────────────────────────
dat <- PF_subject %>% 
  mutate(
    age_z        = as.numeric(scale(age)),
    unimodal_z   = as.numeric(scale(unimodal)),
    transmodal_z = as.numeric(scale(transmodal)),
    BDI_z        = as.numeric(scale(BDI)),
    group        = factor(group),      # “BP”, “HC”
    sex          = factor(sex)         # e.g. “1”, “2”
  )

## ──────────────────────────────────────────────────────────────────────────
## 3.  Model: group-specific σ and ν (Student-t residuals)
## ──────────────────────────────────────────────────────────────────────────
f_BDI <- bf(
  BDI_z ~ unimodal_z + transmodal_z + group + sex + age_z +        # mean slopes
    unimodal_z:group + transmodal_z:group,                   # allow slopes to differ by group
  sigma ~ 0 + group,                                         # σ_BP , σ_HC
  nu    ~ 0 + group                                                # ν_BP , ν_HC
)

priors <- c(
  prior(normal(0, 0.50),            class = "Intercept"),          # mean intercept
  prior(normal(0, 0.30),            class = "b"),                 # all mean slopes & interactions
  prior(normal(0, 0.30),            class = "b", dpar = "sigma"), # log-σ coefficients (one per group)
  prior(normal(log(10), 0.50),      class = "b", dpar = "nu")     # log-ν coefficients (one per group)
)

fit_grp <- brm(
  formula = f_BDI,
  data    = dat,
  family  = student(),          # enables σ & ν sub-models
  prior   = priors,
  chains  = 4, iter = 4000, cores = 4,
  seed    = 123
)

## ──────────────────────────────────────────────────────────────────────────
## 4.  Posterior-predictive check by group (should capture the two peaks)
## ──────────────────────────────────────────────────────────────────────────

fit_grp_prior <- update(fit_grp, sample_prior = "only")  # runs the only prior
pp_check(fit_grp_prior, ndraws = 1e2) + ggtitle("PRIOR")

fit_grp
pp_check(fit_grp, type = "dens_overlay_grouped",
         group = "group", ndraws = 1e2) + ggtitle("POSTERIOR")


pp_check(fit_grp, type = "error_scatter_avg")

posterior <- as_draws_df(fit_grp)
mcmc_areas(
  posterior, # Extract posterior draws
  pars =  names(posterior)[1:13],
  prob = 0.95
) 

##### MIXTURE MODEL TO SOLVE THAT error_scatter_ave

## ──────────────────────────────────────────────────────────────────────────
## 1.  Collapse the four region labels → two levels, then average per subject
## ──────────────────────────────────────────────────────────────────────────
PF_long <- PF_long_list[[5]] |>
  bind_rows() 

PF_subject <- PF_long %>%                                   # long table (all recordings)
  mutate(
    region_type = fct_collapse(                             # uni_* → “unimodal”, trans_* → “transmodal”
      region_type,
      unimodal   = c("uni_self",  "uni_nonself"),
      transmodal = c("trans_self","trans_nonself")
    )
  ) %>% 
  group_by(subj_id, group, sex, age, BDI, region_type) %>%  # keep covariates
  summarise(PF = mean(value), .groups = "drop") %>%         # average PF within subj × region-type
  pivot_wider(names_from = region_type, values_from = PF) %>%   # → columns: unimodal, transmodal
  filter(group == "BP")
## ──────────────────────────────────────────────────────────────────────────
## 2.  Scale continuous variables & code factors
## ──────────────────────────────────────────────────────────────────────────
dat <- PF_subject %>% 
  mutate(
    age_z        = as.numeric(scale(age)),
    unimodal_z   = as.numeric(scale(unimodal)),
    transmodal_z = as.numeric(scale(transmodal)),
    BDI_z        = as.numeric(scale(BDI)),
    group        = factor(group),      # “BP”, “HC”
    sex          = factor(sex)         # e.g. “1”, “2”
  )

## ──────────────────────────────────────────────────────────────────────────
## 3.  Model: group-specific σ and ν (Student-t residuals)
## ──────────────────────────────────────────────────────────────────────────
f_BDI <- bf(
  BDI_z ~ unimodal_z + transmodal_z  + sex + age_z,
  sigma ~ unimodal_z + transmodal_z
)

priors <- c(
  prior(normal(0, 0.5),            class = "Intercept"),          # mean intercept
  prior(normal(0, 1.0),            class = "b")
)

fit_grp <- brm(
  formula = f_BDI,
  data    = dat,
  family  = student(),          # enables σ & ν sub-models
  prior   = priors,
  chains  = 4, iter = 4000, cores = 4,
  seed    = 123
)

## ──────────────────────────────────────────────────────────────────────────
## 4.  Posterior-predictive check by group (should capture the two peaks)
## ──────────────────────────────────────────────────────────────────────────

fit_grp_prior <- update(fit_grp, sample_prior = "only")  # runs the only prior
pp_check(fit_grp_prior, ndraws = 1e2) + ggtitle("PRIOR")

fit_grp
pp_check(fit_grp, type = "dens_overlay_grouped",
         group = "group", ndraws = 1e2) + ggtitle("POSTERIOR")


pp_check(fit_grp, type = "error_scatter_avg")

posterior <- as_draws_df(fit_grp)
mcmc_areas(
  posterior, # Extract posterior draws
  pars =  names(posterior)[1:13],
  prob = 0.95
) 