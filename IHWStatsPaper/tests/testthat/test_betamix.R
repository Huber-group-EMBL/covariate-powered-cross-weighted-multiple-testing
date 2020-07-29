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


