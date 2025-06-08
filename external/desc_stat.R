# This script is for the data preparation for the R 
library(R.matlab)
library(NeuroMyelFC)
library(reshape2)
library(afex)
setwd("~/projects/FrequencySliding/")

###### Matrix preperation #######
BP <- readMat("/home/kaanka5312/projects/FrequencySliding/data/output/PF_BP.mat")
n_bp <- 20
PF_whole_bp <- t(sapply(BP$PF.all.slow4[1:n_bp], function(x) {
  x[[1]]["PF.whole", 1, 1][[1]]
}))

HC <- readMat("/home/kaanka5312/projects/FrequencySliding/data/output/PF_HC.mat")
n_hc <- 33
PF_whole_hc <- t(sapply(HC$PF.all.slow4[1:n_hc], function(x) {
  x[[1]]["PF.whole", 1, 1][[1]]
}))

PF_whole_all <- rbind(PF_whole_bp,PF_whole_hc)

###### Masking #######
mask_list <- list(trans_self = NeuroMyelFC::trans_self,
                  trans_nonself = NeuroMyelFC::trans_nonself,
                  uni_self = NeuroMyelFC::uni_self,
                  uni_nonself = NeuroMyelFC::uni_nonself )

PF_all_topo <- sapply(mask_list, function(mask) {
         rowMeans(PF_whole_all[,mask])
})

Peak_long <- reshape2::melt(PF_all_topo, varnames=c("Subject","Region"))
Peak_long$group <- rep(c(c(rep("BP",n_bp),rep("HC",n_hc))),4)

aov_result <- aov_car(
  value ~ group * Region + Error(Subject/Region),
  data = Peak_long,
  factorize = FALSE  # optional if factors are already set
)

summary(aov_result)

library(emmeans)
# Store the emmeans object
em <- emmeans(aov_result, ~ Region | group)

# Get pairwise comparisons WITHIN each group (e.g., BP: trans_self vs uni_self, etc.)
pairs(em)


library(ggplot2)
ggplot(Peak_long, aes(x = Region, y = value, color = group, group = group)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               position = position_dodge(width = 0.2)) +
  labs(y="Peak Frequency", title = "Slow4")+
  theme_minimal()

