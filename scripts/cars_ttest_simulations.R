library(doParallel)
registerDoParallel(cores=14)


library(doRNG)
library(tidyverse)
library(IHWStatsPaper)


cars_sim_combs <- expand.grid(m = 10000,
                              k = c(30, 70, 161, 373, 864, 2000),
                              eff_size = 0.5,
                              mu_x_2 = 0.25,
                              seed = 1:400)


set.seed(1)

r <- foreach(i=1:nrow(cars_sim_combs)) %dorng% {
  print(i);
  IHWStatsPaper:::eval_cars_sim(
         cars_sim_combs$m[i],
         cars_sim_combs$k[i],
         cars_sim_combs$eff_size[i],
         cars_sim_combs$mu_x_2[i],
         cars_sim_combs$seed[i])}

saveRDS(r, file="cars_simulations.Rds")
