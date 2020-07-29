library(IHWStatsPaper)
library(dplyr)

group_sim_combs <- expand.grid(m = 20000,
                               K_coarse= c(2, 5, 10, 20, 40),
                               pi0_global = 0.9,
                               seed = 1:2)

group_sim_res <- bind_rows(mapply(IHWStatsPaper:::eval_grouped_sim,
                                  group_sim_combs$m,
                                  group_sim_combs$K,
                                  group_sim_combs$pi0_global,
                                  4,
                                  group_sim_combs$seed,
                                  SIMPLIFY=FALSE ))
