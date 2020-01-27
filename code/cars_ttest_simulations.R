library(doParallel)
registerDoParallel(cores=16)

library(tidyverse)
source("multiple_testing_methods.R")
source("cars_ttest_functions.R")


#ks_ballpark <- exp(seq(log(10), to=log(2000), length.out=5))

cars_sim_combs <- expand.grid(m = 10000,
                              k = c(30, 70, 161, 373, 864, 2000),
                              eff_size = 0.5,
                              mu_x_2 = 0.25,
                              seed = 1:400)



source("multiple_testing_methods.R")
source("cars_ttest_functions.R")


r <- foreach(i=1:nrow(cars_sim_combs)) %dopar% {
  print(i);
  mapply(eval_cars_sim,
         cars_sim_combs$m[i],
         cars_sim_combs$k[i],
         cars_sim_combs$eff_size[i],
         cars_sim_combs$mu_x_2[i],
         cars_sim_combs$seed[i],
         SIMPLIFY=FALSE )}

saveRDS(r, file="cars_simulations.Rds")

res_r <-  bind_rows(unlist(r, recursive = FALSE)) %>% 
           group_by(method, m, k, mu_x_1, mu_x_2) %>%
           summarize(FDR = mean(FDP), Power=mean(pow)) %>%
           arrange(k, mu_x_1, mu_x_2, desc(Power)) %>% 
           ungroup() 


