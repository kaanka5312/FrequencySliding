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

saveRDS(fit_list, file = "./data/output/slow4_PF_brms_models.rds")

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

library(ggpubr)
combined_plot <- ggarrange(plotlist = plot_list,common.legend = TRUE)
final_plot <- annotate_figure(combined_plot,
                top = text_grob("Slow 4", face = "bold", size = 16))
# Save to file
ggsave("./figures/slow4_combined_plot.png", final_plot, width = 10, height = 8, dpi = 300,bg="white")
