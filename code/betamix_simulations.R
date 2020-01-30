library(doParallel)
registerDoParallel(cores=16)

library(tidyverse)
library(adaptMT)
source("multiple_testing_methods.R")
source("betamix_functions.R")


betamix_sim_combs <- rbind(expand.grid(mu_slope = seq(1, 3, length=7),
                              seed = 1:400),
                           data.frame(mu_slope =1, seed=401:1000)) 
                     # more simulations for betamix in first case to stabilize Monte Carlo FDP estimate


betamix_r <- foreach(i=1:nrow(betamix_sim_combs)) %dopar% {
  print(i);
  mapply(eval_beta_unif_sim,
         betamix_sim_combs$mu_slope[i],
         betamix_sim_combs$seed[i],
         lfdr_only=FALSE,
         SIMPLIFY=FALSE )}

saveRDS(betamix_r, file="betamix_simulations.Rds")

betamix_res <-  bind_rows(unlist(betamix_r, recursive = FALSE)) %>% 
  group_by(method, m, prob_one_sided) %>%
  summarize(FDR = mean(FDP), Power=mean(pow)) %>%
  arrange(prob_one_sided, desc(Power)) %>% 
  ungroup() 
