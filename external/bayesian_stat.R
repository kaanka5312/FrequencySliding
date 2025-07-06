# This script is for the data preparation for the R 
library(R.matlab)
library(NeuroMyelFC)
library(reshape2)
library(afex)

# MacOS
setwd("~/projects/FrequencySliding/")

# windows
setwd("C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/")


###### Matrix preperation #######
library(R.matlab)

# Read files
BP <- readMat("./data/output/PF_BP.mat")
HC <- readMat("./data/output/PF_HC.mat")

n_bp <- 38
n_hc <- 33

# Define all component names
components <- c("PF.peak", "PF.trough", "PF.rise", "PF.fall", "PF.whole")

# Initialize result list
PF_all_scaled <- list()

for (comp in components) {
  # Extract component from each subject (1st dimension is subject, 2nd is region)
  PF_bp_mat <- t(sapply(BP$PF.all.slow4[1:n_bp], function(x) {
    x[[1]][comp, 1, 1][[1]]
  }))
  
  PF_hc_mat <- t(sapply(HC$PF.all.slow4[1:n_hc], function(x) {
    x[[1]][comp, 1, 1][[1]]
  }))
  
  # Combine and scale
  PF_all_scaled[[comp]] <- scale(rbind(PF_bp_mat, PF_hc_mat))
}


library(tidyverse)

# Add region type info
# Your region labels
region_names <- paste0("V", 1:360)
category_vector <- integer(360)  # initialize with zeros

category_vector[NeuroMyelFC::uni_nonself] <- 1
category_vector[NeuroMyelFC::uni_self] <- 2
category_vector[NeuroMyelFC::trans_nonself] <- 3
category_vector[NeuroMyelFC::trans_self] <- 4

region_types <- data.frame(
  region = region_names,
  region_type = factor(category_vector,
                       labels = c("uni_nonself", "uni_self", "trans_nonself", "trans_self"))  # adjust if needed
)

# Initialize the result list
PF_long_list <- list()

# Loop over each PF component
for (comp in names(PF_all_scaled)) {
  df <- as.data.frame(PF_all_scaled[[comp]])
  colnames(df) <- region_names
  df$subj_id <- 1:71
  df$group <- rep(c("BP", "HC"), times = c(38, 33))
  
  # Convert to long format
  df_long <- pivot_longer(df, cols = all_of(region_names),
                          names_to = "region", values_to = "value")
  
  # Add region_type
  df_long <- left_join(df_long, region_types, by = "region")
  
  # Save in the list
  PF_long_list[[comp]] <- df_long
}

PF_collapsed_list <- lapply(PF_long_list, function(df) {
  df %>%
    mutate(region_category = case_when(
      str_starts(region_type, "uni") ~ "unimodal",
      str_starts(region_type, "trans") ~ "transmodal",
      TRUE ~ NA_character_
    ))
})
#########
library(brms)

# Initialize empty list to hold model fits
fit_list <- list()

# Loop through each component's long data frame
for (comp in names(PF_long_list)) {
  message("Fitting model for: ", comp)
  
  fit_list[[comp]] <- brm(
    value ~ group * region_type + (1 | subj_id) + (1 | region),
    data = PF_long_list[[comp]],
    family = gaussian(),
    chains = 4,
    iter = 2000,
    cores = 6,
    seed = 123  # Set seed for reproducibility
  )
}

# Loop through each component's long data frame
for (comp in names(PF_collapsed_list)) {
  message("Fitting model for: ", comp)
  
  fit_list[[comp]] <- brm(
    value ~ group * region_category + (1 | subj_id) + (1 | region),
    data = PF_collapsed_list[[comp]],
    family = gaussian(),
    chains = 4,
    iter = 2000,
    cores = 6,
    seed = 123  # Set seed for reproducibility
  )
}

saveRDS(fit_list, file = "./data/output/unitrans_slow4_PF_brms_models.rds")

##
ce_plot <- plot(conditional_effects(fit_list[[1]], effects = "region_type:group"), plot=FALSE)[[1]]
ce_plot + ggtitle("test")

plot_list <- lapply(seq_along(fit_list), function (x) {
  plot(conditional_effects(fit_list[[x]], effects = "region_type:group"), plot=FALSE)[[1]] + 
    ggtitle(names(PF_all_scaled)[x]) +
    scale_x_discrete(labels = c(
      "1" = "Uni-Nonself",
      "2" = "Uni-Self",
      "3" = "Trans-Nonself",
      "4" = "Trans-Self"
    )) +
    ylab("Standardized Peak Frequency")
})

plot_list <- lapply(seq_along(fit_list), function (x) {
  plot(conditional_effects(fit_list[[x]], effects = "region_category:group"), plot=FALSE)[[1]] + 
    ggtitle(names(PF_all_scaled)[x]) +
    ylab("Standardized Peak Frequency")
})


library(ggpubr)
combined_plot <- ggarrange(plotlist = plot_list,common.legend = TRUE)
final_plot <- annotate_figure(combined_plot,
                top = text_grob("Slow 5", face = "bold", size = 16))
# Save to file
ggsave("./figures/slow5_combined_unitrans_plot.png", final_plot, width = 10, height = 8, dpi = 300,bg="white")

###### RELATION TO BDI #####

PF_summary_list <- lapply(seq_along(PF_long_list), function(x) {PF_long_list[[x]] %>%
    group_by(subj_id, region_type) %>%
    summarise(mean_PF = mean(value), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = region_type, values_from = mean_PF)} )

subjects <- readLines("./data/raw/bpb_strings.txt")
subjects_split <- c(strsplit(subjects, "\t")[[1]],paste0("sub-",43:75,".results"))

# Extract numbers from the strings
numbers <- as.numeric(gsub("sub-(\\d+)\\.results", "\\1", subjects_split))

# Sort based on the numbers
subjects_sorted <- subjects_split[order(numbers)]

print(subjects_sorted)

PF_summary_list <- lapply(seq_along(PF_summary_list), function(x) {
  PF_summary_list[[x]]$subj_id <- subjects_sorted
  PF_summary_list[[x]]  # Return the modified element
  })

PF_summary_list_processed <- lapply(PF_summary_list, process_pf_summary)


###############

subject_data_list <- lapply(seq_along(PF_summary_list_processed), function(x) {
  PF_summary_list_processed[[x]] %>% select(subj_id, 
                                            group, 
                                            BDI_total_score, 
                                            NA.Score, 
                                            PD.Score, 
                                            SM.Score,  
                                            starts_with("zdev_"))
})

subject_data_list <- lapply(seq_along(subject_data_list), function(x) {
  subject_data_list[[x]] %>%
    mutate(
      zdev_trans_mean = (zdev_trans_self + zdev_trans_nonself) / 2,
      zdev_uni_mean = (zdev_uni_self + zdev_uni_nonself) / 2
    )
})
####### PRIOR PREDICTIVE CHECK ############
# For student distribution univariate total BDI
priors <- c(
  prior(normal(30, 5), class = "Intercept"),   # More tightly around typical BDI
  prior(normal(0, 2), class = "b"),            # Smaller plausible effects per SD
  prior(exponential(1), class = "sigma")    # Less noise allowed
)

fit_prior <- brm(
  formula = BDI_total_score  ~ zdev_uni_nonself + zdev_uni_self + zdev_trans_nonself + zdev_trans_self,
  data = subject_data %>% filter(group == "BP"),
  family = student(),
  prior = priors,
  sample_prior = "only",
  chains = 4,
  iter = 2000,
  cores = 4
)

prior_preds <- posterior_predict(fit_prior)

library(bayesplot)

ppc_dens_overlay(
  y = subject_data$BDI_total_score[subject_data$group == "BP"],
  yrep = prior_preds[1:1000, ]  # Use first 100 draws to avoid overplotting
)

######## DUE TO HIGH COLLINEARITY, GATHERED UNDER UNI-TRANS #####
priors <- c(
  prior(normal(30, 5), class = "Intercept"),   # More tightly around typical BDI
  prior(normal(0, 2), class = "b"),            # Smaller plausible effects per SD
  prior(exponential(1), class = "sigma")    # Less noise allowed
)

fit_prior <- brm(
  formula = BDI_total_score  ~ zdev_uni_mean + zdev_trans_mean,
  data = subject_data_list[[5]] %>% filter(group == "BP"),
  family = student(),
  prior = priors,
  sample_prior = "only",
  chains = 4,
  iter = 2000,
  cores = 4
)

prior_preds <- posterior_predict(fit_prior)

library(bayesplot)

ppc_dens_overlay(
  y = subject_data_list[[5]]$BDI_total_score[subject_data_list[[5]]$group == "BP"],
  yrep = prior_preds[1:1000, ]  # Use first 100 draws to avoid overplotting
)


fit_unv_list <- lapply(seq_along(subject_data_list), function(x){
  brm(
    BDI_total_score  ~ zdev_uni_mean + zdev_trans_mean,
    data = subject_data_list[[x]] %>% filter(group=="BP"),
    family = student(),
    prior = priors,
    chains = 4,
    iter = 2000,
    cores = 4
  )
})

summary(fit_unv_list[[5]])

# Extract posterior draws
posterior <- as_draws_df(fit_unv_list[[4]])

# Forest plot for fixed effects
p1 <- mcmc_areas(
  posterior,
  pars = c("b_Intercept", "b_zdev_uni_nonself", 
           "b_zdev_uni_self", "b_zdev_trans_nonself", "b_zdev_trans_self"),
  prob = 0.95
) +
  ggtitle("Peak Frequency - Fall") 

ggsave("./figures/PF_n22_Fall.png", p1, width = 6, height = 4, dpi = 300,bg="white")

###########
# Multivariate Model
bp <- subject_data_list[[1]][subject_data_list[[1]]$group=="BP",]
cor(bp[, c("zdev_uni_nonself", "zdev_uni_self", "zdev_trans_nonself", "zdev_trans_self")])

bp <- subject_data_list[[1]][subject_data_list[[1]]$group=="BP",]
cor(bp[, c("zdev_uni_mean", "zdev_trans_mean")])

fit_multv_list <- lapply(seq_along(subject_data_list), function(x){
  brm(
  mvbind(NA.Score, PD.Score, SM.Score) ~ 1 + zdev_uni_nonself + zdev_uni_self + zdev_trans_nonself + zdev_trans_self,
  data = subject_data_list[[x]] %>% filter(group=="BP"),
  family = student(),
  prior =  priors <- c(
    prior(normal(0, 5), class = "Intercept", resp = "NAScore"),
    prior(normal(0, 5), class = "Intercept", resp = "PDScore"),
    prior(normal(0, 5), class = "Intercept", resp = "SMScore"),
    
    prior(normal(0, 2), class = "b", resp = "NAScore"),
    prior(normal(0, 2), class = "b", resp = "PDScore"),
    prior(normal(0, 2), class = "b", resp = "SMScore"),
    
    prior(exponential(1), class = "sigma", resp = "NAScore"),
    prior(exponential(1), class = "sigma", resp = "PDScore"),
    prior(exponential(1), class = "sigma", resp = "SMScore"),
    
    prior(lkj(2), class = "rescor")  # Residual correlation prior
  ),
  chains = 4,
  iter = 2000,
  cores = 4
)
})

fit_multv_list <- lapply(seq_along(subject_data_list), function(x){
  brm(
    mvbind(NA.Score, PD.Score, SM.Score) ~ 1 + zdev_uni_mean + zdev_trans_mean,
    data = subject_data_list[[x]] %>% filter(group=="BP"),
    family = student(),
    prior =  priors <- c(
      prior(normal(0, 5), class = "Intercept", resp = "NAScore"),
      prior(normal(0, 5), class = "Intercept", resp = "PDScore"),
      prior(normal(0, 5), class = "Intercept", resp = "SMScore"),
      
      prior(normal(0, 2), class = "b", resp = "NAScore"),
      prior(normal(0, 2), class = "b", resp = "PDScore"),
      prior(normal(0, 2), class = "b", resp = "SMScore"),
      
      prior(exponential(1), class = "sigma", resp = "NAScore"),
      prior(exponential(1), class = "sigma", resp = "PDScore"),
      prior(exponential(1), class = "sigma", resp = "SMScore"),
      
      prior(lkj(2), class = "rescor")  # Residual correlation prior
    ),
    chains = 4,
    iter = 2000,
    cores = 4
  )
})

summary(fit_subject)

fit_multv_list[[5]]

# Hierarhical Modelling 
library(dplyr)
library(tidyr)
long_data_list <- lapply(seq_along(subject_data_list), function(x) {
  subject_data_list[[x]] %>%
    pivot_longer(
      cols = c(NA.Score, PD.Score, SM.Score),
      names_to = "subscale",
      values_to = "BDI_score"
    )
  })

fit_hierarchical_list <- lapply(seq_along(subject_data_list), function(x){
  brm(
    BDI_score ~ 1 + zdev_uni_nonself + zdev_uni_self + zdev_trans_nonself + zdev_trans_self + 
      (1 + zdev_uni_nonself + zdev_uni_self + zdev_trans_nonself + zdev_trans_self | subscale),
    data = long_data_list[[x]] %>% filter(group=="BP"),
    family = student(),
    prior = c(
      prior(normal(0, 5), class = "Intercept"),
      prior(normal(0, 2), class = "b"),
      prior(exponential(1), class = "sd"),
      prior(lkj(2), class = "cor")  # Allows correlations between subscale-level slopes
    ),
    chains = 4,
    iter = 2000,
    cores = 4
  )
})

fit_hierarchical_list <- lapply(seq_along(subject_data_list), function(x){
  brm(
    BDI_score ~ 1 + zdev_uni_mean + zdev_trans_mean + 
      (1 + zdev_uni_mean + zdev_trans_mean | subscale),
    data = long_data_list[[x]] %>% filter(group=="BP"),
    family = student(),
    prior = c(
      prior(normal(0, 5), class = "Intercept"),
      prior(normal(0, 2), class = "b"),
      prior(exponential(1), class = "sd"),
      prior(lkj(2), class = "cor")  # Allows correlations between subscale-level slopes
    ),
    chains = 4,
    iter = 2000,
    cores = 4
  )
})


summary(fit_hierarchical_list[[5]])

###  P L O T ####
library(bayesplot)
hier_plot_list <- lapply(seq_along(fit_hierarchical_list), function(x){
  posterior <- as_draws_df(fit_hierarchical_list[[x]])
  # Select regression coefficients (exclude intercepts if desired)
  params_to_plot <- c("b_zdev_uni_mean", "b_zdev_trans_mean")
  mcmc_areas(posterior, pars = params_to_plot, prob = 0.95) +
    ggtitle(components[x]) +
    theme_minimal() + 
    geom_vline(xintercept = 0, color = "red", linetype = "dashed")
})

library(ggpubr)
combined_plot <- ggarrange(plotlist = hier_plot_list,common.legend = TRUE,nrow = 1,ncol = 5)
final_plot <- annotate_figure(combined_plot,
                              top = text_grob("Slow 5", face = "bold", size = 16))
# Save to file
ggsave("./figures/slow5_Freq_n38.png", final_plot, width = 20, height = 8, dpi = 300,bg="white")

####### CHECK WITH LOO ######
library(loo)
loo_mv <- loo(fit_subject)
loo_hier <- loo(fit_hierarchical)


