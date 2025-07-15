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

library(dplyr)

dat_cls <- PF_subject %>%                       # ← your collapsed table
  mutate(
    subj_id       = as.factor(subj_id),
    group_bin     = if_else(group == "BP", 1L, 0L),  # 1 = BP, 0 = HC
    unimodal_z    = as.numeric(scale(unimodal)),
    transmodal_z  = as.numeric(scale(transmodal))
  ) %>%
  select(subj_id, group_bin, unimodal_z, transmodal_z, sex, age)

library(brms)

f_cls <- bf(
  group_bin ~ unimodal_z * transmodal_z + sex + age,  # no intercept → brms adds it
  family = bernoulli(link = "logit")
)

priors <- c(
  prior(normal(0, 1.0),  class = "Intercept"),  # log-odds baseline
  prior(normal(0, 1.0),  class = "b")           # slopes
)

fit_cls <- brm(
  f_cls, data = dat_cls,
  prior  = priors, sample_prior = "only",
  chains = 4, iter = 4000, cores = 4,
  seed   = 2025
)

pp_check(fit_cls, ndraws = 1e2)
full_fit <- update(fit_cls, sample_prior = "yes")

p_draws <- posterior_epred(full_fit, newdata = dat_cls)   # S × N matrix
prob_BP <- colMeans(p_draws)                             # mean over draws
dat_cls$prob_BP <- prob_BP

dat_cls$pred_label <- if_else(prob_BP > 0.5, 1L, 0L)
table(Predicted = dat_cls$pred_label, True = dat_cls$group_bin)


library(pROC)

auc_obj <- roc(dat_cls$group_bin, prob_BP)
auc_obj$auc          # posterior mean AUC, point estimate
ci.auc(auc_obj, conf.level = 0.95)        # 95 % CrI by DeLong method
plot(auc_obj)

# Leave-one-out accuracy
set.seed(123)
# ── 1.  helper: fit-once-per-subject  ─────────────────────────────────
loso_probs <- map_dfr(unique(dat_cls$subj_id), function(id) {
  
  train_dat <- filter(dat_cls, subj_id != id)
  test_dat  <- filter(dat_cls, subj_id == id)
  
  fit_i <- brm(
    f_cls, data = train_dat,
    prior  = priors,
    chains = 4, iter = 4000, cores = 4,   # light settings: 71 fits!
    seed   = 123, refresh = 0, silent = 2
  )
  
  p_i <- posterior_epred(fit_i, newdata = test_dat) |> mean()
  
  tibble(subj_id = id,
         true    = test_dat$group_bin,
         prob_BP = p_i)
})

# ── 2.  evaluation  ───────────────────────────────────────────────────
loso_probs <- loso_probs %>%
  mutate(pred    = if_else(prob_BP > 0.5, 1L, 0L))

# confusion matrix & accuracy
print(table(Predicted = loso_probs$pred, True = loso_probs$true))
acc <- mean(loso_probs$pred == loso_probs$true)
cat("LOSO accuracy:", round(acc * 100, 1), "%\n")

# ROC / AUC with 95 % CI
auc_obj <- roc(loso_probs$true, loso_probs$prob_BP)
print(auc_obj$auc)
print(ci.auc(auc_obj, conf.level = 0.95))

library(readr)          # install.packages("readr") if needed

write_csv(loso_probs, "./data/output/loso_probs.csv")

#########
#### LONG FORMAT SUCKED ######
df_long <- PF_collapsed_list[[5]] %>%                    # already long
  mutate(group_bin = (group == "BP") * 1) # 1 = BP, 0 = HC


f_cls_long <- bf(
  group_bin ~ value * region_category +        # PF value & its interaction
    age + sex +                     # covariates
    (region_category | subj_id) +   # varying intercept & slope
    (1 | region),                   # repeated parcels
  family = bernoulli(link = "logit")
)

priors_long <- c(
  prior(normal(0, 1),  class = "Intercept"),
  prior(normal(0, 0.7), class = "b"),          # slopes
  prior(exponential(2), class = "sd")          # random‐effect SDs
)

fit_long <- brm(
  f_cls_long, data = df_long,
  prior  = priors_long, sample_prior = "only",
  chains = 4, iter = 4000, cores = 4, seed = 123
)

update(fit_long, sample_prior = "yes")

# ── 3.  posterior probabilities & AUC  ───────────────────────────────
# one mean probability per row
p_hat   <- posterior_epred(fit_long, ndraws = 200) |> colMeans()
df_long$prob_BP <- p_hat

# aggregate to subject-level (optional)
df_subj <- df_long %>%
  group_by(subj_id, group_bin) %>%
  summarise(prob_BP = mean(prob_BP), .groups = "drop")

auc_obj <- roc(df_subj$group_bin, df_subj$prob_BP)
print(auc_obj$auc)                        # point estimate
ci.auc(auc_obj, conf.level = 0.95)        # 95 % CrI by DeLong method

# ── 4.  posterior-predictive check (density overlay)  ────────────────
pp_check(fit_long, type = "bars_grouped",
         group = "group_bin", ndraws = 200) +
  ggtitle("Posterior predictive: BP vs HC class counts")