library(doParallel)
library(tidyverse)
library(IHWStatsPaper)

betamix_sim_combs <- rbind(expand.grid(mu_slope = seq(2, 3, length=2),
                                       seed = 1:2))


betamix_r <- foreach(i=1:nrow(betamix_sim_combs)) %dopar% {
  print(paste0("simulation run:",i));
  mapply(eval_beta_unif_sim,
         betamix_sim_combs$mu_slope[i],
         betamix_sim_combs$seed[i],
         lfdr_only=FALSE,
         SIMPLIFY=FALSE )}

betamix_r <- bind_rows(betamix_r)
#tmp <- group_by(betamix_r,  mu_slope, method) %>% summarize(n=n(), FDR = mean(FDP), se_FDR = sd(FDP)/sqrt(n), Power = mean(pow))
#saveRDS(betamix_r, file="betamix_simulations.Rds")

betamix_sim_combs <- expand.grid(mu_slope = 2,
                                 seed = 1:2,
                                 prob_onesided = seq(0,0.1, length=2))



betamix_onesided_r <- foreach(i=1:nrow(betamix_sim_combs)) %dopar% {
  print(paste0("simulation run:",i));
  mapply(eval_beta_unif_sim,
         betamix_sim_combs$mu_slope[i],
         betamix_sim_combs$seed[i],
         TRUE,
         betamix_sim_combs$prob_onesided[i],
         lfdr_only=FALSE,
         SIMPLIFY=FALSE )}

betamix_onesided_r <- bind_rows(betamix_onesided_r)

