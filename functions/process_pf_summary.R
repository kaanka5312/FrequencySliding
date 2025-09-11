process_pf_summary <- function(PF_summary, bdi_file = "./data/raw/BDI_ThreeFactor_MultiSubject.xlsx") {
  library(readxl)
  library(dplyr)
  
  # Load and preprocess BDI data
  bdi <- readxl::read_xlsx(bdi_file, sheet=1) %>% 
    select(c("SubjectID", "NA Score", "PD Score", "SM Score")) %>% 
    mutate(across(c(`NA Score`, `PD Score`, `SM Score`), as.numeric)) %>%
    drop_na()
  
  bdi$SubjectID <- paste0(bdi$SubjectID, ".results")
  colnames(bdi) <- c("SubjectID", "NA.Score", "PD.Score", "SM.Score")
  
  # Add HC subjects with zero scores
  hc_frame <- data.frame(
    SubjectID = paste0("sub-", 43:75, ".results"),
    `NA.Score` = 0,
    `PD.Score` = 0,
    `SM.Score` = 0
  )
  
  bdi_all <- rbind(bdi, hc_frame)
  
  # Merge with PF_summary
  merged_data <- PF_summary %>%
    left_join(bdi_all, by = c("subj_id" = "SubjectID")) %>%
    mutate(group = c(rep("BP", 38), rep("HC", 33))) %>%
    drop_na()
  
  # Total BDI score and Z-score
  merged_data <- merged_data %>%
    mutate(BDI_total_score = NA.Score + PD.Score + SM.Score,
           BDI_total_z = scale(BDI_total_score))
  
  # HC stats
  hc_stats <- merged_data %>%
    filter(group == "HC") %>%
    summarise(across(c(uni_nonself, uni_self, trans_nonself, trans_self),
                     list(mean = mean, sd = sd), na.rm = TRUE),
              mean_BDI = mean(BDI_total_score, na.rm = TRUE),
              sd_BDI = sd(BDI_total_score, na.rm = TRUE))
  
  # Z-normalize to HC
  merged_data <- merged_data %>%
    mutate(across(c(uni_nonself, uni_self, trans_nonself, trans_self),
                  ~ (.x - hc_stats[[paste0(cur_column(), "_mean")]]) / hc_stats[[paste0(cur_column(), "_sd")]],
                  .names = "zdev_{.col}"))
  
  return(merged_data)
}
