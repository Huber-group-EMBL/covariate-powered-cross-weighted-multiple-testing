library(tidyverse)


source("multiple_testing_methods.R")
source("grouped_simulations_functions.R")


#------------------------------------------------
#          Main  simulations
#------------------------------------------------
group_sim_combs <- expand.grid(m = 20000,
                                 K_coarse= c(2, 5, 10, 20, 40),
                                 pi0_global = 0.9,
                                 seed = 1:200)

group_sim_res <- bind_rows(mapply(eval_grouped_sim,
                                    group_sim_combs$m,
                                    group_sim_combs$K,
                                    group_sim_combs$pi0_global,
                                    4,
                                    group_sim_combs$seed,
                                    SIMPLIFY=FALSE ))

saveRDS(group_sim_res, file= "grouped_simulations_signal.Rds")



