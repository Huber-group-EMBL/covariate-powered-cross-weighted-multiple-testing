library(tidyverse)
library(IHWStatsPaper)
library(testthat)

set.seed(10)
cars_data <- cars_ttest_data(m=10000, k=200, n_x=10000, n_y=10000, mu_x_1 = 0.5/10,
                             mu_x_2 = 0.25, mu_y_2 = 0.25)

cars_data_preprocess <- with(cars_data, CARS_preprocess(x, y, var_mat))

expect_true( max(abs( 2*(1-pnorm(abs(cars_data_preprocess$t_1))) - cars_data_preprocess$pvalue)) < 1e-10)

pv2 <- 2*(1-pnorm(abs(cars_data_preprocess$t_2)))
grps <- IHW::groups_by_filter(pv2, 3)
cars_data_preprocess$grps <- grps
cars_data_preprocess$secondary_pvalue <- pv2
#ggplot(cars_data_preprocess, aes(x=pvalue)) + geom_histogram(origin=0,binwidth=0.05) + facet_grid(.~grps)

set.seed(10)
main_cars <- CARS::CARS(cars_data$x, cars_data$y, 0.2, tau=0.9, variance=cars_data$var_mat, option="regular")
set.seed(10)
custom_cars <- IHWStatsPaper:::CARS_extended(cars_data_preprocess$t_1, cars_data_preprocess$t_2, 0.2, tau=0.9, option = "regular")

expect_equal(custom_cars$de, main_cars$de)
expect_true(abs(custom_cars$th - main_cars$th) <= 0.001)

expect_equal(as.logical(main_cars$de), oracle_local_fdr_test(cars_data_preprocess$pvalue, main_cars$cars, 0.2))
cars_fun_tmp <- custom_cars$cars_fun
cars_fun_tmp_eval <- cars_fun_tmp(cars_data_preprocess$t_1, cars_data_preprocess$t_2)
expect_equal(cars_fun_tmp_eval, custom_cars$cars)

t1_seq <- seq(1,5,length=100)
cars_fun_at_2 <- cars_fun_tmp(t1_seq, 2)
cars_fun_at_1 <- cars_fun_tmp(t1_seq, 1)
cars_fun_at_3 <- cars_fun_tmp(t1_seq, 3)
#plot(t1_seq, cars_fun_at_2)
#lines(t1_seq, cars_fun_at_1)
#lines(t1_seq, cars_fun_at_3, col='red')

set.seed(10)
ihw_cars_res <- ihw_cars(cars_example_preprocess$t_1, cars_example_preprocess$t_2, 0.2, return_weights=TRUE, kfolds=2, Storey=FALSE)
set.seed(10)
ihw_cars_rjs <- ihw_cars(cars_example_preprocess$t_1, cars_example_preprocess$t_2, 0.2, return_weights=FALSE, kfolds=2, Storey=FALSE)


expect_equal(length(unique(ihw_cars_res$folds)),2)
fold1_idx = ihw_cars_res$folds==1
fold2_idx = ihw_cars_res$folds==2


set.seed(10)
ihw_cars_refit <- ihw_cars(cars_example_preprocess$t_1, cars_example_preprocess$t_2, 0.2, return_weights=TRUE, Storey=FALSE, folds=ihw_cars_res$folds)
set.seed(10)
cars_wts_fold1 <- cars_weighter(cars_example_preprocess$t_1[fold2_idx], cars_example_preprocess$t_2[fold2_idx], cars_example_preprocess$t_2[fold1_idx],
                                tau=0.5, alpha=0.2)
cars_wts_fold2 <- cars_weighter(cars_example_preprocess$t_1[fold1_idx], cars_example_preprocess$t_2[fold1_idx], cars_example_preprocess$t_2[fold2_idx],
                                             tau=0.5, alpha=0.2)


expect_true(abs(sum(cars_wts_fold1) - sum(fold1_idx)) <= 0.0001)
expect_true(abs(sum(cars_wts_fold2) - sum(fold2_idx)) <= 0.0001)

all_wts <- rep(NA, 10000)
all_wts[fold1_idx] <- cars_wts_fold1
all_wts[fold2_idx] <- cars_wts_fold2

expect_true(max(abs(all_wts-ihw_cars_refit$ws)) <= 0.0000001)

expect_equal(tau_weighted_bh(cars_example_preprocess$pvalue, all_wts, 0.5) <= 0.2, ihw_cars_refit$rjs)

set.seed(1)
ihw_cars_storey_res <-ihw_cars(cars_example_preprocess$t_1, cars_example_preprocess$t_2, 0.2, return_weights=TRUE, Storey=TRUE, folds=ihw_cars_res$folds)

expect_true(sum(ihw_cars_storey_res$rjs) > sum(ihw_cars_refit$rjs))

expect_true(abs(sum(ihw_cars_storey_res$ws[fold1_idx]) - sum(fold1_idx)/ihw_cars_storey_res$storey_pi0s[1]) <= 0.0001)
expect_true(abs(sum(ihw_cars_storey_res$ws[fold2_idx]) - sum(fold2_idx)/ihw_cars_storey_res$storey_pi0s[2]) <= 0.0001)

ihw_df <- data.frame(Ps = cars_example_preprocess$pvalue,
                     Ws = ihw_cars_refit$ws,
                     Ws_storey = ihw_cars_storey_res$ws,
                     folds = ihw_cars_refit$folds,
                     adj_p = ihw_cars_refit$adj_p,
                     rjs = ihw_cars_refit$rjs)

ihw_df_summary <- group_by(ihw_df, folds) %>% summarize(sum_Ws = sum(Ws), n=n())
expect_equal(ihw_df_summary$sum_Ws, ihw_df_summary$n)
