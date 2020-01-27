library(doParallel)
#registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
registerDoParallel(cores=4)

library(tidyverse)
source("multiple_testing_methods.R")
source("grouped_simulations_functions.R")

group_sim_combs <- expand.grid(m = 20000,
                               K_coarse= c(2, 5, 10, 20, 40),
                               pi0_global = 1.0,
                               seed = 1:12000)



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


summarize_group_null_sim_res <- group_by(tmp_r, method, K_coarse) %>%
  summarize(FDR = mean(FDP), Power=mean(pow), n_monte_carlo = n()) %>%
  arrange(K_coarse, desc(Power)) %>% 
  ungroup() %>%
  filter(method %in% method_names) %>%
  mutate(method = factor(method, levels=method_names)) 

library(cowplot)
fdr_grouped_null_plot <- ggplot(summarize_group_null_sim_res, aes(x=K_coarse, y=FDR,shape=method, col=method)) + 
  geom_line() + geom_point() + 
  scale_x_log10(breaks=c(2,5, 10, 20,40)) + 
  scale_color_manual(values=method_colors)+
  scale_shape_manual(values=method_shapes)+
  geom_hline(yintercept=0.1, linetype=2)+
  xlab("Number of groups (G)") +
  theme_cowplot() +
  theme(legend.position = "none")
fdr_grouped_null_plot
ggsave(fdr_grouped_null_plot, filename="fdr_grouped_null_plot.pdf", width=3.4,height=2.8)


