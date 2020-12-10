library(IHWStatsPaper)
library(dplyr)
library(testthat)


m <- 20000
K_coarse <- 40
pi0_global <- 0.9

set.seed(1)
sim_tmp <- IHWStatsPaper::grouped_sim(m, 2, pi0_global, 4)
expect_equal(length(unique(sim_tmp$Xs)),2)

set.seed(1)
sim <- IHWStatsPaper::grouped_sim(m, K_coarse, pi0_global, 4)

ihw_gbh_res <- ihw_gbh(sim$Ps, sim$Xs, 0.1, Storey=FALSE)
sum(ihw_gbh_res)
ihw_gbh_storey_res <- ihw_gbh(sim$Ps, sim$Xs, 0.1, Storey=TRUE)
sum(ihw_gbh_storey_res)
fdp_eval(sim$Hs, ihw_gbh_storey_res)

ihw_res <- IHWStatsPaper::ihw_nmeth_wrapper(sim$Ps, sim$Xs, 0.1, Storey=TRUE)

ihw_res_direct <- IHW::ihw(sim$Ps, as.factor(sim$Xs), 0.1, null_proportion=TRUE, null_proportion_level=0.5)
IHW::rejections(ihw_res_direct) == sum(ihw_res)

sum((1-sim$Hs)*ihw_res)  == fdp_eval(sim$Hs,ihw_res)$FDP * sum(ihw_res)
sum(sim$Hs*ihw_res)  == fdp_eval(sim$Hs,ihw_res)$pow * sum(sim$Hs)
fdp_eval(sim$Hs,ihw_res)
bh_rjs <- p.adjust(sim$Ps, method="BH") <= 0.1
sum(bh_rjs)
fdp_eval(sim$Hs, bh_rjs)

oracle_sabha <- p.adjust(sim$Ps*sim$pi0s, method="BH") <= 0.1
sum(oracle_sabha)

dd_sabha <- groupwise_sabha(sim$Ps, as.factor(sim$Xs), 0.1)
sum(dd_sabha)
fdp_eval(sim$Hs, dd_sabha)


gbh_storey_res <- gbh_simple(sim$Ps, sim$Xs, 0.1, Storey=TRUE)
sum(gbh_storey_res)
fdp_eval(sim$Hs, gbh_storey_res)


strat_bhq <- stratified_bhq(sim$Ps, sim$Xs, 0.1)
sum(strat_bhq)
fdp_eval(sim$Hs, strat_bhq)

clfdr_res <- stratified_clfdr(sim$Ps, sim$Xs, 0.1)
fdp_eval(sim$Hs, clfdr_res)

set.seed(1)
collect_res <- IHWStatsPaper:::eval_grouped_sim(20000,40,0.9,4,seed=1,alpha=0.1)

expect_equal(filter(collect_res, method=="BH")$rjs, sum(bh_rjs))
expect_equal(filter(collect_res, method=="Clfdr")$rjs, sum(clfdr_res))
expect_equal(filter(collect_res, method=="SABHA")$rjs, sum(dd_sabha))
expect_equal(filter(collect_res, method=="IHW-Grenander-Storey")$rjs, sum(ihw_res))
expect_equal(filter(collect_res, method=="GBH-Storey")$rjs, sum(gbh_storey_res))
expect_equal(filter(collect_res, method=="IHW-GBH-Storey")$rjs, sum(ihw_gbh_storey_res))
expect_equal(filter(collect_res, method=="IHW-GBH")$rjs, sum(ihw_gbh_res))
expect_equal(filter(collect_res, method=="SBH")$rjs, sum(strat_bhq))

set.seed(100)
sim_null <- grouped_global_null_sim(10000, 1000)
ihw_gbh_null <- ihw_gbh(sim_null$Ps, sim_null$Xs, 0.2, Storey=FALSE)
gbh_null <- gbh_simple(sim_null$Ps, sim_null$Xs, 0.2, Storey=FALSE)

set.seed(100)
eval_null <- IHWStatsPaper:::eval_grouped_global_null_sim(10000, 1000, 1, alpha=0.2, gbh_only=TRUE)
expect_equal(eval_null$FDP, c(1,0))
expect_equal(eval_null$rjs, c(sum(gbh_null), sum(ihw_gbh_null)))
