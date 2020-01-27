### colors used throughout
###
beyonce_colors <- c("#b72da0", "#7c5bd2", "#0097ed","#00c6c3",
                    "#9cd78a", "#f7f7a7", "#ebab5f", "#e24344",
                    "#04738d")#,"#d8cdc9")
beyonce_colors[6] <- c("#dbcb09") # thicker yellow
pretty_colors <- beyonce_colors[c(2,1,3:5)]


method_names <- c("BH", "SBH", "Clfdr", "GBH-Storey", "IHW-GBH-Storey", "IHW-Grenander-Storey", "SABHA")
method_colors <- c("Black", beyonce_colors[c(3:5, 2, 1, 6)])
names(method_colors) <- method_names

method_shapes <- c(20, 0,1,2, 12, 13,2 )
names(method_shapes) <- method_names


###############33
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

source("multiple_testing_methods.R")
source("grouped_simulations_functions.R")

# Settings to run through

tmp <- grouped_sim(20000, 2, 1.00)

tst <- groupwise_sabha(tmp$Ps, tmp$Xs, 0.1)
sum(tst)
tmp_ihw <- ihw_gbh(tmp$Ps, tmp$Xs, 0.1, Storey=TRUE, return_weights=TRUE)

#-------------------------------------------
#--------  Null simulation from Intro
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

summarize_group_null_sim_res <- group_by(group_null_res, method, K) %>%
  summarize(FDR = mean(FDP))

saveRDS(fdr_tbl, file= "global_null_bh_vs_ihwgbh_fdr_control.Rds")


fdr_null_plot <- ggplot(summarize_group_null_sim_res, 
                        aes(x=K, y=FDR,col=method, shape=method)) + geom_line() +
  geom_point() +
  xlab("Number of groups (G)") + 
  ylim(0, 0.4) + 
  geom_hline(yintercept=0.2, linetype=2) + 
  scale_x_log10() + 
  scale_color_manual(values=pretty_colors[1:2]) + 
  theme_cowplot()

fdr_null_plot

ggsave(fdr_null_plot, filename="intro_fdr_null_plot.pdf", width=6,height=2.5)

#---------------------------------------------
# More exhaustive null simulations
#---------------------------------------------

Ks <- c(5, 10, 25, 50, 100, 200)
m <- 10000
seeds <- 1:5000
group_null_combs <- expand.grid(m = m, K=Ks, seed = seeds)

group_null_res <- bind_rows(mapply(eval_grouped_global_null_sim, group_null_combs$m,
                                   group_null_combs$K, group_null_combs$seed, alpha=0.2, gbh_only=TRUE,
                                   SIMPLIFY=FALSE ))

summarize_group_null_sim_res <- group_by(group_null_res, method, K) %>%
  summarize(FDR = mean(FDP))

#------------------------------------------------
# Main  simulations
#------------------------------------------------
group_sim_combs <- expand.grid(m = 20000,
                                 K_coarse= c(2, 5, 10, 20, 40),
                                 pi0_global = 0.9,
                                 seed = 1:200)

set.seed(1)
group_sim_res <- bind_rows(mapply(eval_grouped_sim,
                                    group_sim_combs$m,
                                    group_sim_combs$K,
                                    group_sim_combs$pi0_global,
                                    4,
                                    group_sim_combs$seed,
                                    SIMPLIFY=FALSE ))
saveRDS(group_sim_res, file= "grouped_simulations_signal.Rds")
group_sim_res <- readRDS(file= "grouped_simulations_signal.Rds")

  
summarize_group_sim_res <- group_by(group_sim_res, method, K_coarse) %>%
                             summarize(FDR = mean(FDP), Power=mean(pow)) %>%
                             arrange(K_coarse, desc(Power)) %>% 
                             ungroup() %>%
                             filter(method %in% method_names) %>%
                             mutate(method = factor(method, levels=method_names))


fdr_grouped_signal_plot <- ggplot(summarize_group_sim_res, aes(x=K_coarse, y=FDR,shape=method, col=method)) + 
  geom_line() + geom_point() + 
  scale_color_manual(values=method_colors)+
  scale_shape_manual(values=method_shapes)+
  geom_hline(yintercept=0.1, linetype=2)+
  xlab("Number of groups (G)")+ 
  scale_x_log10(breaks=c(2,5, 10, 20,40)) +   
  theme_cowplot() +
  theme(legend.position = "none")
  
fdr_grouped_signal_plot
ggsave(fdr_grouped_signal_plot, filename="fdr_grouped_signal_plot.pdf", width=3.4,height=2.8)




power_grouped_signal_plot <- ggplot(summarize_group_sim_res, aes(x=K_coarse, y=Power, shape=method, col=method)) + 
                                geom_line() + geom_point() + 
                                scale_color_manual(values=method_colors) + 
                                scale_shape_manual(values=method_shapes) +
                                xlab("Number of groups (G)")  +
                                scale_x_log10(breaks=c(2,5, 10, 20,40)) +   
                                theme_cowplot() 
ggsave(power_grouped_signal_plot, filename="power_grouped_signal_plot.pdf", width=5.6,height=2.8)



tmp <- grouped_sim(800, 10, 0.75)


alpha_ADMM <- 10^2; 
beta <- 10^3;
eta <- 5;
max_iters <- 5000;
converge_thr <- 1e-4 # parameters for ADMM
ADMM_params <- c(alpha_ADMM,beta,eta,max_iters,converge_thr)

tmp2 <- groupwise_sabha(tmp$Ps, tmp$Xs, max_iters=100, beta=1000)

sum( (tmp$Ps > 0.5)/tmp2$q_all)
tmp_pi0s <- Solve_q_block(tmp$Ps, 0.5, 0.1, tmp$Xs, ADMM_params)



pi0_df <- mutate(tmp) %>% 
                 group_by(Xs) %>% 
                 summarize(pi0 = 1 - mean(Hs)) %>% mutate(pi0_sabha = tmp2)
source("groupwise_sabha_pi0.R")
source("multiple_testing_methods.R")

pi0_df2 <- mutate(tmp, sabha_pi0s = tmp_pi0s) %>% 
               group_by(Xs) %>% 
               summarize(pi0 = 1 - mean(Hs), 
                         sabha_pi0 = mean(sabha_pi0s), 
                         var_sabha = sd(sabha_pi0s)) %>%
              mutate(pi02 = tmp2)

