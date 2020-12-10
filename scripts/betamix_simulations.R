library(doParallel)
registerDoParallel(cores=14)
library(doRNG)

library(tidyverse)
library(IHWStatsPaper)


betamix_sim_combs <- expand.grid(mu_slope = seq(1, 3, length=7),
                                 seed = 1:400)


set.seed(1)
betamix_r <- foreach(i=1:nrow(betamix_sim_combs)) %dorng% {
  print(paste0("simulation run:",i));
  eval_beta_unif_sim(
         betamix_sim_combs$mu_slope[i],
         betamix_sim_combs$seed[i],
         lfdr_only=FALSE)}

saveRDS(betamix_r, file="betamix_simulations.Rds")
