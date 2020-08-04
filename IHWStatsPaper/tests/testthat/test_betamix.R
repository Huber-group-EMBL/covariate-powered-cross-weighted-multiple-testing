library(dplyr)
library(testthat)
library(IHWStatsPaper)

set.seed(1)
m <- 5000
alpha_dbeta <- 0.2
Ps_alt <- rbeta(m, alpha_dbeta, 1)

# part 1, sanity check
Ys_alt <- -log(Ps_alt)
glm_beta <- glm(Ys_alt ~ 1, family=Gamma(link="inverse"))
summary(glm_beta)

# part 2, check predictions
Ps_all <- sample( c(Ps_alt,  runif(m) ) )

# part 3, check EM algorithm without censoring
Xs_all <- data.frame(X1 = rep(1,2*m), X2=runif(2*m))
fit1<- gamma_glm_basic_em(Ps_all, Xs_all, formula_rhs="~X1")
expect_true( abs(fit1$alphas[1] - alpha_dbeta) <= 0.1)
expect_true( abs(fit1$pi1s[1] - 0.5) <= 0.05)
expect_equal(sd(fit1$alphas), 0)
expect_equal(sd(fit1$pi1s), 0)

fit2 <-  gamma_glm_basic_em(Ps_all, Xs_all)
expect_true( abs(fit2$alphas[1] - alpha_dbeta) <= 0.1)
expect_lt(sd(fit2$alphas), 0.01)
expect_true( abs(fit2$pi1s[1] - 0.5) <= 0.05)
expect_lt(sd(fit2$pi1s), 0.01)

# Check AdaPT EM algorithm
fit_adapt <- adapt_mtp(Ps_all, Xs_all, 0.3, formula_rhs="~X1", return_fit=TRUE)
expect_true(all(fit_adapt$fit$rejs == which(fit_adapt$rjs)))

# Check IHW \tau-censoring EM algorithm
tau_censor1 <- 1e-17
Ps_all_censor1 <- ifelse(Ps_all >= tau_censor1, Ps_all, 0)
fit1_censor1 <- gamma_glm_censored_em(Ps_all_censor1, Xs_all, tau_censor1, formula_rhs="~X1")

expect_true( abs(fit1$alphas[1] - fit1_censor1$alphas[1]) <= 0.02)
expect_true( abs(fit1$pi1s[1] -fit1_censor1$pi1s[1]) <= 0.01)

fit1_censor_nocovariates <- gamma_glm_censored_em(Ps_all_censor1, Xs_all, tau_censor1, formula_rhs="~1")
fit1_nocovariates <- gamma_glm_basic_em(Ps_all, Xs_all, formula_rhs="~1")

expect_equal( var(fit1_censor_nocovariates$pi1s), 0 )
expect_equal( var(fit1_censor_nocovariates$alphas), 0 )

expect_equal( var(fit1_nocovariates$pi1s), 0 )
expect_equal( var(fit1_nocovariates$alphas), 0 )

# Check local FDR code
lfdrs_nocovariates_model <- get_localfdrs_betamix(Ps_all, fit1_censor_nocovariates$pi1s, fit1_censor_nocovariates$alphas)
expect_equal(order(Ps_all), order(lfdrs_nocovariates_model))
rjs_nocovariates <- oracle_local_fdr_test(Ps_all, lfdrs_nocovariates_model, 0.2)
expect_true( mean(lfdrs_nocovariates_model[rjs_nocovariates]) <= 0.2)
max_t <- max(Ps_all[rjs_nocovariates])
tail_fdr <- get_tailfdr_betamix(max_t, fit1_censor_nocovariates$pi1s, fit1_censor_nocovariates$alphas, numerator_bh = FALSE)
tail_fdr2 <- get_tailfdr_betamix(max_t, fit1_censor_nocovariates$pi1s[1], fit1_censor_nocovariates$alphas[1], numerator_bh = FALSE)
expect_equal(tail_fdr, tail_fdr2)
expect_true(abs(tail_fdr - 0.2) <= 0.01)

tail_fdr_bh  <-  get_tailfdr_betamix(max_t, fit1_censor_nocovariates$pi1s[1], fit1_censor_nocovariates$alphas[1], numerator_bh = TRUE)
tail_fdr_bh2 <- get_tailfdr_betamix(rep(max_t, 10000), fit1_censor_nocovariates$pi1s, fit1_censor_nocovariates$alphas, numerator_bh = TRUE)
expect_equal(tail_fdr_bh, tail_fdr_bh2)
expect_true(tail_fdr_bh >= tail_fdr)


# betamix datadriven lfdr
lfdrs_nocovariates_nocensor_model <- get_localfdrs_betamix(Ps_all, fit1_nocovariates$pi1s, fit1_nocovariates$alphas)
rjs_nocovariates_nocensor <- oracle_local_fdr_test(Ps_all, lfdrs_nocovariates_nocensor_model, 0.2)

lfdr_betamix_em_rjs <- betamix_datadriven_lfdr(Ps_all, Xs_all,0.2, formula_rhs="~1", maxiter=50 )
expect_equal(lfdr_betamix_em_rjs, rjs_nocovariates_nocensor)

max_t_nocensor <- max(Ps_all[rjs_nocovariates_nocensor])
lfdr_nocensor <- get_localfdrs_betamix(max_t_nocensor, fit1_nocovariates$pi1s[1], fit1_nocovariates$alphas[1])

# test the threshold function and the weighting function
t_t <- get_thresholds_betamix(lfdr_nocensor, fit1_nocovariates$pi1s[1], fit1_nocovariates$alphas[1])
expect_equal(as.numeric(t_t), max_t_nocensor)



# Check betamix simulation
betamix_sim <- beta_unif_sim(20000)
manual_lfdrs <- get_localfdrs_betamix(betamix_sim$Ps, betamix_sim$pi1s, betamix_sim$alphas)
expect_equal(manual_lfdrs, betamix_sim$oracle_lfdrs)

oracle_rjs <- oracle_local_fdr_test(betamix_sim$Ps, betamix_sim$oracle_lfdrs, 0.2)
oracle_rjs2 <- betamix_oracle_lfdr(betamix_sim$Ps, betamix_sim$pi1s, betamix_sim$alphas, 0.2)
expect_equal(oracle_rjs, oracle_rjs)

expect_true( mean(manual_lfdrs[oracle_rjs]) <= 0.2)


all_ts <- get_thresholds_betamix(0.1, betamix_sim$pi1s, betamix_sim$alphas)
ts_to_lfdrs <- get_localfdrs_betamix(all_ts, betamix_sim$pi1s, betamix_sim$alphas)
expect_true(all( abs(ts_to_lfdrs - 0.1) <= 0.00001))
