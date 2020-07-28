library(dplyr)
library(testthat)
library(IHWStatsPaper)
#group_sim_combs <- expand.grid(m = 20000,
#                               K_coarse= c(2, 5, 10, 20, 40),
#                               pi0_global = 0.9,
#                               seed = 1:200)

#group_sim_res <- bind_rows(mapply(eval_grouped_sim,
#                                  group_sim_combs$m,
#                                  group_sim_combs$K_coarse,
#                                  group_sim_combs$pi0_global,
#                                  4,
#                                  group_sim_combs$seed,
#                                  SIMPLIFY=FALSE ))

m <- 10000
K <- 40
null_sim_run <- grouped_global_null_sim(m, K)
expect_equal(length(unique(null_sim_run$Xs)), 40)
expect_equal(as.numeric(table(null_sim_run$Xs)), rep(m/K, K))

grouped_sim_run <- grouped_sim(20000, 10, 0.9, 4)
expect_equal(length(unique(grouped_sim_run$Xs)), 10)
expect_equal(length(unique(grouped_sim_run$Xs_tilde)), 40)

pi0s_df <- group_by(grouped_sim_run, Xs_tilde) %>% summarize(pi0 = min(pi0s), pi0max = max(pi0s),
                                                             mu = min(eff_sizes), mumax = max(eff_sizes))
expect_equal(pi0s_df$pi0, pi0s_df$pi0max)
expect_equal(pi0s_df$mu, pi0s_df$mumax)

try_pi0s <- (0.2  + 0.8*(0:39)/36)
expect_true(  all(pi0s_df$pi0 == 1 | ( abs( pi0s_df$pi0 - try_pi0s) < 1e-10 )) )

try_effsizes <- 2.5 - 2*(0:39)/36
expect_true(  all(pi0s_df$mu == 0 | ( abs( pi0s_df$mu - try_effsizes) < 1e-10 )) )




grouped_sim_run_5 <- grouped_sim(100000, 2, 0.9, 4)
expect_equal(length(unique(grouped_sim_run_5$Xs)), 2)
expect_equal(length(unique(grouped_sim_run_5$Xs_tilde)), 40)

sabha_res <- groupwise_sabha(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs, 0.1, return_fit=TRUE, max_iters=400 )
qs <- sabha_res$q
expect_true(sum( 2*(grouped_sim_run_5$Ps >= 0.5)/sabha_res$q_all) <= 100000.001)
expect_true(all(qs >= 0.1))
expect_true(all(qs <= 1.0))

true_pi0s <- group_by(grouped_sim_run_5, Xs) %>% summarize(pi0 = mean(1-Hs)) %>% mutate(sabha_qs = sabha_res$q)
expect_equal( as.numeric(rank(true_pi0s$pi0)), as.numeric(rank(true_pi0s$sabha_qs)))

sabha_res2 <- groupwise_sabha(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs_tilde, 0.1, return_fit=TRUE, max_iters=400 )
expect_true(sum( 2*(grouped_sim_run_5$Ps >= 0.5)/sabha_res2$q_all) <= 100000.001)
expect_true(all(sabha_res2$q >= 0.1))
expect_true(all(sabha_res2$q <= 1.0))
true_pi0s2 <- group_by(grouped_sim_run_5, Xs_tilde) %>% summarize(pi0 =min(pi0s)) %>% mutate(sabha_qs = sabha_res2$q)

grouped_sim_run_5 <- mutate(grouped_sim_run_5, sabha_q = sabha_res2$q_all)
grouped_sim_run_5_summary <- group_by(grouped_sim_run_5, Xs_tilde) %>% summarize(sabha_qmin = min(sabha_q) , sabha_qmax = max(sabha_q))
expect_equal(grouped_sim_run_5_summary$sabha_qmax, grouped_sim_run_5_summary$sabha_qmin)



# Test stratified BHq

sbhq_res <- stratified_bhq(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs, 0.05)
sbhq_manual <- group_by(grouped_sim_run_5, Xs) %>% summarize(rjs = sum(p.adjust(Ps, method="BH") <= 0.05))
expect_equal(sum(sbhq_res), sum(sbhq_manual$rjs))
expect_false(sum(sbhq_res) == sum(p.adjust(grouped_sim_run_5$Ps, method="BH") <= 0.05))
