library(doParallel)
registerDoParallel(cores=16)

library(tidyverse)
library(adaptMT)
source(file="multiple_testing_methods.R")
source("betamix_functions.R")


betamix_sim_combs <- expand.grid(mu_slope = 2,
                                 seed = 1:400,
                                 prob_onesided = seq(0,0.1, length=6))



betamix_onesided_r <- foreach(i=1:nrow(betamix_sim_combs)) %dopar% {
  print(i);
  mapply(eval_beta_unif_sim,
         betamix_sim_combs$mu_slope[i],
         betamix_sim_combs$seed[i],
         TRUE,
         betamix_sim_combs$prob_onesided[i],
         lfdr_only=FALSE,
         SIMPLIFY=FALSE )}

saveRDS(betamix_onesided_r, file="betamix_onesided_simulations.Rds")


