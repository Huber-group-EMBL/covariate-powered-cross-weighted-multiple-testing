library(doParallel)
registerDoParallel(cores=14)

library(doRNG)
library(tidyverse)
library(IHWStatsPaper)

#-------------------------------------------
#--------  Null simulation from Section 2
#-------------------------------------------

Ks <- c(5, 10, 25, 50, 100, 200, 500,1000)
m <- 10000
seeds <- 1:12000 # #not formal random seeds, just for printing
group_null_combs <- expand.grid(m = m, K=Ks, seed = seeds)

set.seed(1)
group_null_res <- foreach(i=1:nrow(group_null_combs)) %dorng% {
                      print(paste0("simulation run:",i));
                      IHWStatsPaper:::eval_grouped_global_null_sim(
                            group_null_combs$m[i],
                            group_null_combs$K[i],
                            group_null_combs$seed[i],
                            alpha=0.2,
                            gbh_only=TRUE) }

saveRDS(group_null_res, file= "group_null_res_intro.Rds")
