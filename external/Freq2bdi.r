# This script is for the data preparation for the R 
library(R.matlab)
library(NeuroMyelFC)
library(reshape2)
library(afex)
library(tidyverse)

# MacOS
#setwd("~/projects/FrequencySliding/")

# windows
setwd("C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/")


###### Matrix preperation #######
bp <- readMat("./data/output/PF_BP.mat")
n_bp <- 38
bp_rows <- scan("./bpb_strings.txt", what = "character", sep = "\t")

# Use sub() to replace ".results" with an empty string
bp_rows_new <- sub(".results", "", bp_rows)
pf_whole_bp <- list()
for (i in seq_along(names(bp))) {
  pf_whole_bp[[i]] <- t(sapply(bp[[i]][1:n_bp], function(x) {
    x[[1]]["PF.whole", 1, 1][[1]]
  })) %>% as.data.frame()
}
names(pf_whole_bp) <- names(bp)

pf_whole_bp <- lapply(pf_whole_bp, function(x) {
  # Assign the new row names to each data frame (df) in the list
  rownames(x) <- bp_rows_new
  # Return the modified data frame
  return(x)
})

pf_whole_bp_filtered <- lapply(pf_whole_bp, function(x) {
  x <- x %>%
    filter(!rownames(x) %in% c("sub-20", "sub-30", "sub-35", "sub-36", "sub-37", "sub-41", "sub-42"))
  return(x)
})

# Parsing through
trans <- NeuroMyelFC::trans_nonself | NeuroMyelFC::trans_self
uni <- !trans

mask_list <- list(trans = trans, uni = uni)

# Apply each mask to each dataset
result_nested <- lapply(pf_whole_bp_filtered, function(data) {
  lapply(mask_list, function(mask) {
    data.frame(
      mean_score = rowMeans(data[, mask])
    ) %>%
      rownames_to_column(var = "SubjectID")
  })
})

###### BDI anhedonia correlation ##########
library(readxl)
bdi_data <- read_excel(
  "C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/data/raw/BDI_ThreeFactor_MultiSubject.xlsx"
) %>%
  mutate(cog = as.numeric(`Item 3`) + # rumination
           as.numeric(`Item 5`) + # rumination
           as.numeric(`Item 6`) + # rumination
           as.numeric(`Item 7`) + # rumination
           as.numeric(`Item 8`) + # rumination
           as.numeric(`Item 11`) + # anxiety
           as.numeric(`Item 17`)) %>% # anxiety
  mutate(dep = as.numeric(`Item 4`) + # anhedonia
      as.numeric(`Item 12`) + # anhedonia
      as.numeric(`Item 21`) + # anhedonia
      as.numeric(`Item 1`) # depression
  )

dep_df <- data.frame(
  SubjectID = bdi_data$SubjectID,
  dep = bdi_data$dep,
  cog = bdi_data$cog
)

## Transmodal Correlation
library(dplyr)
library(tidyr)

# Assuming your nested list is named result_nested
# and the external dataframe is named dep_df

result_joined <- lapply(result_nested, function(sub_list) {
  lapply(sub_list, function(df) {
    df %>%
      left_join(dep_df, by = "SubjectID") %>%
      drop_na()
  })
})

library(dplyr)

# Assuming your list is named result_joined

correlation_results <- lapply(result_joined, function(sub_list) {
  lapply(sub_list, function(df) {
    # This inner loop iterates over 'dep' and 'cog'
    lapply(c("dep", "cog"), function(var_name) {
      # Perform the correlation test
      test_result <- cor.test(df$mean_score, df[[var_name]])
      
      # Extract the r and p-value and store them in a tibble
      tibble(
        variable = var_name,
        r = as.numeric(test_result$estimate),
        p_value = test_result$p.value
      )
    }) %>% bind_rows() # Combine the results for this dataframe
  })
})

# View the results
print(correlation_results)

