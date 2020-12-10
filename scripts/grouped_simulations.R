library(doParallel)
registerDoParallel(cores=14)

library(doRNG)
library(tidyverse)
library(IHWStatsPaper)


group_sim_combs <- expand.grid(m = 20000,
                               K_coarse= c(2, 5, 10, 20, 40),
                               pi0_global = 0.9,
                               seed = 1:200)

set.seed(1)
group_sim_res <-foreach(i=1:nrow(group_sim_combs)) %dorng% {
                                  print(i);
                                  IHWStatsPaper:::eval_grouped_sim(
                                    group_sim_combs$m[i],
                                    group_sim_combs$K[i],
                                    group_sim_combs$pi0_global[i],
                                    4,
                                    group_sim_combs$seed[i])}

saveRDS(group_sim_res, file= "grouped_simulations_signal.Rds")
