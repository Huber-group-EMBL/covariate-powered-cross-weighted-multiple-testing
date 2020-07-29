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


# Test stratified CLfdr and some of the local FDR code

clfdr_res <- stratified_clfdr(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs, 0.05)
manual_local_fdrs <- group_by(grouped_sim_run_5, Xs) %>% mutate(local_fdrs = fdrtool::fdrtool(Ps, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr)

expect_lt(mean(manual_local_fdrs$local_fdrs[clfdr_res]), 0.05)
filtered_table <- filter(manual_local_fdrs, Xs == 1)
expect_equal( order(filtered_table$local_fdrs, filtered_table$Ps), order(filtered_table$Ps))

expect_false( all( order(manual_local_fdrs$local_fdrs, manual_local_fdrs$Ps) == order(manual_local_fdrs$Ps)))
expect_equal( oracle_local_fdr_test(manual_local_fdrs$Ps, manual_local_fdrs$local_fdrs, 0.05), clfdr_res)

# Test tau-weighted code
Ws <- rep(1, 100000)
tbh_res <- tau_weighted_bh(grouped_sim_run_5$Ps, Ws, 0.1)
expect_true( all( (tbh_res <= 0.1) ==  (p.adjust(grouped_sim_run_5$Ps, method="BH") <= 0.1)))

expect_equal(weighted_storey_pi0(grouped_sim_run_5$Ps, Ws), (1+sum(grouped_sim_run_5$Ps >= 0.5))*2/100000)
expect_equal(weighted_storey_pi0(grouped_sim_run_5$Ps, Ws, tau=0.9), (1+sum(grouped_sim_run_5$Ps >= 0.9))*10/100000)


Ws <- 2*(grouped_sim_run_5$Xs == 1)
tbh_res2 <- tau_weighted_bh(grouped_sim_run_5$Ps, Ws, 0.1)

expect_equal( tbh_res2[grouped_sim_run_5$Xs == 1] <= 0.1, p.adjust(grouped_sim_run_5$Ps[grouped_sim_run_5$Xs == 1], method="BH") <= 0.1)
expect_equal( tbh_res2[grouped_sim_run_5$Xs == 2], rep(1,50000) )
expect_equal(weighted_storey_pi0(grouped_sim_run_5$Ps, Ws), (1+sum(grouped_sim_run_5$Ps[grouped_sim_run_5$Xs == 1] >= 0.5))*2/50000)


tbh_res3 <- tau_weighted_bh(grouped_sim_run_5$Ps, Ws, 1.0)
expect_equal( tbh_res3[grouped_sim_run_5$Xs == 1], p.adjust(grouped_sim_run_5$Ps[grouped_sim_run_5$Xs == 1], method="BH"))

# Test GBH

gbh_res <- gbh_simple(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs, 0.05)
gbh_res_adaptive <- gbh_simple(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs, 0.05, Storey=TRUE)

expect_true(sum(gbh_res_adaptive) > sum(gbh_res))

# try to see if it is consistent with an approach based on grouped Storey tester.

ws_gbh <- grouped_storey_weighter(grouped_sim_run_5$Ps, grouped_sim_run_5$Xs, grouped_sim_run_5$Xs, 0.5, 0.2)
ws_gbh <- ws_gbh*length(ws_gbh)/sum(ws_gbh)
expect_equal(gbh_res, tau_weighted_bh(grouped_sim_run_5$Ps, ws_gbh, 0.5) <= 0.05)

pi0_st_bh <- weighted_storey_pi0(grouped_sim_run_5$Ps, ws_gbh)
expect_equal(gbh_res_adaptive, tau_weighted_bh(grouped_sim_run_5$Ps*pi0_st_bh, ws_gbh, 0.5) <= 0.05)
expect_equal(gbh_res_adaptive, tau_weighted_bh(grouped_sim_run_5$Ps, ws_gbh/pi0_st_bh, 0.5) <= 0.05)


