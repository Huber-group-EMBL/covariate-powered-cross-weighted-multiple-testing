library(doParallel)
registerDoParallel(cores=4)

library(tidyverse)
source("multiple_testing_methods.R")
source("grouped_simulations_functions.R")





#-------------------------------------------
#--------  Null simulation from Section 2
#-------------------------------------------

Ks <- c(5, 10, 25, 50, 100, 200, 500,1000)
m <- 10000
seeds <- 1:8000
group_null_combs <- expand.grid(m = m, K=Ks, seed = seeds)

set.seed(1)
group_null_res <- bind_rows(mapply(eval_grouped_global_null_sim, group_null_combs$m,
                                   group_null_combs$K, group_null_combs$seed, alpha=0.2, gbh_only=TRUE,
                                   SIMPLIFY=FALSE ))

saveRDS(group_null_res, file= "group_null_res_intro.Rds")




#-------------------------------------------
#--------  Null simulation from Section 5
#-------------------------------------------

group_sim_combs <- expand.grid(m = 20000,
                               K_coarse= c(2, 5, 10, 20, 40),
                               pi0_global = 1.0,
                               seed = 1:10000)



tmp_r <- foreach(i=1:nrow(group_sim_combs)) %dopar% {
          print(i);
          mapply(eval_grouped_sim,
                group_sim_combs$m[i],
                group_sim_combs$K[i],
                group_sim_combs$pi0_global[i],
                4,
                group_sim_combs$seed[i],
                SIMPLIFY=FALSE )}


saveRDS(tmp_r, file="grouped_null_simulations_all_methods.Rds", )  


