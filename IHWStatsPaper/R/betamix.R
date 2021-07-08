#' Local fdr for beta-mixture model
#'
#' @param ts      Numeric vector of thresholds at which to compute the local fdr
#' @param pi1s   Numeric vector with probability i-th hypothesis is an alternative
#' @param alphas  Numeric vector of first parameter of alternative Beta distribution (Beta(a_i, 1)).
#'
#' @return Vector of local fdrs.
#' @export
get_localfdrs_betamix <- function(ts, pi1s, alphas){
  pi0s <- 1 - pi1s
  pi0s/(pi0s + pi1s*dbeta(ts, alphas, 1))
}

#' Tail fdr for beta-mixture model
#'
#' @param ts      Numeric vector of thresholds at which to compute the local fdr
#' @param pi1s   Numeric vector with probability i-th hypothesis is an alternative
#' @param alphas  Numeric vector of first parameter of alternative Beta distribution (Beta(a_i, 1)).
#' @param numerator_bh  If TRUE (default is FALSE) then replace pi0 in numerator by 1.
#'
#' @return Numeric with Tail fdr value.
#' @export
get_tailfdr_betamix <- function(ts, pi1s, alphas, numerator_bh=FALSE){
  m <- length(ts)
  pi0s <- 1 - pi1s
  if (numerator_bh){
    pi0s_num <- rep(1,m)
  } else {
    pi0s_num <- pi0s
  }
  num <- sum(pi0s_num * ts)
  denom <- sum(pi0s*ts + pi1s*pbeta(ts, alphas, 1))
  num/denom
}

#' Thresholds for beta-mixture model
#'
#' @param c       Numeric constant, the contour level ofthe local fdrs
#' @param pi1s   Numeric vector with probability i-th hypothesis is an alternative
#' @param alphas  Numeric vector of first parameter of alternative Beta distribution (Beta(a_i, 1)).
#'
#' @return Vector of thresholds t_i so that lfdr_i(t_i) = c
#' @export
get_thresholds_betamix <- function(c, pi1s, alphas){
  pi0s <- 1 - pi1s
  ((pi0s/c - pi0s)/(alphas*pi1s))^(1/(alphas-1))
}

#' Oracle local false discovery rate procedure in the beta mixture model
#'
#' @param Ps  Numeric vector of unadjusted p-values.
#' @param pi1s  True probability of being an alternative in the beta mixture model
#' @param mu_alphas Parameter vector of alternative Beta(mu_alpha, 1) distribution of each test.
#' @param alpha Significance level at which to apply method
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#' @export
betamix_oracle_lfdr <- function(Ps, pi1s, mu_alphas, alpha){
  oracle_local_fdr_test(Ps,
                        get_localfdrs_betamix(Ps, pi1s, mu_alphas),
                        alpha)
}



#' Local false discovery rate procedure with local fdrs estimated from the Betamix-model through the EM algorithml
#'
#' @param Ps     Numeric vector of unadjusted p-values.
#' @param Xs     Data frame with features
#' @param alpha  Significance level at which to apply method
#' @param formula_rhs Formula defining the RHS in the fitted GLMs, defaults to ~X1+X2 (used in simulations herein).
#' @param maxiter  Total number of iterations to run the EM algorithm
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#' @export
betamix_datadriven_lfdr <- function(Ps, Xs, alpha, formula_rhs="~X1+X2", maxiter=200,...){
  gamma_glm_fit  <- gamma_glm_basic_em(Ps, Xs, formula_rhs=formula_rhs, maxiter = maxiter, tau_pi0=0.5,...)
  betamix_oracle_lfdr(Ps, gamma_glm_fit$pi1s, gamma_glm_fit$alphas, alpha)
}


#' Fit conditional beta-uniform mixture model using EM algorithm.
#'
#' @param Ps       Numeric vector of p-values
#' @param Xs       Data frame with features
#' @param formula_rhs Formula defining the RHS in the fitted GLMs, defaults to ~X1+X2 (used in simulations herein).
#' @param maxiter  Total number of iterations to run the EM algorithm
#' @param tau_pi0  Number in (0,1) used for the initialization of the EM algorithm through the Boca-Leek approach (default: 0.5).
#' @param pi1_min  Numeric in (0,1),  lower bound on conditionally probability of being an alternative (default: 0.01).
#' @param pi1_max  Numeric in (0,1), upper bound on conditionally probability of being a null (default: 0.9).
#' @param alpha_min  Numeric in (0,1), upper bound on parameter of alternative Beta distribution (default: 0.1).
#' @param alpha_max  Numeric in (0,1), upper bound on parameter of alternative Beta distribution (default: 0.9).

#' @return List with parameters of fitted conditional Beta-Uniform mixture
#' @export
gamma_glm_basic_em <- function(Ps, Xs, formula_rhs="~X1+X2", maxiter = 50, tau_pi0=0.5,
                               pi1_min = 0.01, pi1_max = 0.9, alpha_min=0.1, alpha_max=0.9){
  # basic transform
  Ys <- -log(Ps)
  m <- length(Ys)

  # initialize EM iterations

  ## Use Boca-Leek to initialize pi0s:
  BL_transform <- as.integer(Ps >= tau_pi0)
  glm_bl <- glm(as.formula(paste("BL_transform", formula_rhs)), family=binomial(), data=Xs)
  pi1s_iter <- 1 - predict(glm_bl, type="response")/(1-tau_pi0)
  pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))

  # initial guess for HS is 1-FDR, use BH + our pilot pi0
  Hs_iter <- rep(NA, m)
  Ps_adj_tmp <- p.adjust(Ps, method="BH")
  Hs_iter <- 1 - Ps_adj_tmp*(1-pi1s_iter)

  for (i in 1:maxiter){

    # fit GLM on transformed responses
    glm_mus_iter <- glm( formula(paste("Ys", formula_rhs)), family=Gamma(), weights=Hs_iter, data=Xs )
    # update alphas and Hs
    alphas_iter <- predict(glm_mus_iter, type="link")
    alphas_iter <- pmax(alpha_min, pmin( alphas_iter, alpha_max))
    Hs_iter <- 1-(1 - pi1s_iter)/((1 - pi1s_iter)  + pi1s_iter*dbeta(Ps, alphas_iter, 1))

    # fit logistic GLM
    glm_pi1s_iter <- glm(formula(paste("Hs_iter", formula_rhs)), family=quasibinomial(), data=Xs,
                         mustart = pi1s_iter)
    # get logistic GLM predictions
    pi1s_iter <-  predict(glm_pi1s_iter, type="response")
    pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))
  }
  list(glm_mus = glm_mus_iter, glm_pi1s = glm_pi1s_iter, glm_bl = glm_bl,
       alphas = alphas_iter, pi1s = pi1s_iter)
}



#' Fit conditional beta-uniform mixture model using EM algorithm under p-value tau-censoring.
#'
#' @param censored_Ps       Numeric vector of tau-censored p-values
#' @param Xs       Data frame with features
#' @param tau_censor Numeric, level at which p-valeus have been tau-censored
#' @param formula_rhs Formula defining the RHS in the fitted GLMs, defaults to ~X1+X2 (used in simulations herein).
#' @param maxiter  Total number of iterations to run the EM algorithm.
#' @param tau_pi0  Number in (0,1) used for the initialization of the EM algorithm through the Boca-Leek approach (default: 0.5).
#' @param pi1_min  Numeric in (0,1),  lower bound on conditionally probability of being an alternative (default: 0.01).
#' @param pi1_max  Numeric in (0,1), upper bound on conditionally probability of being a null (default: 0.9).
#' @param alpha_min  Numeric in (0,1), upper bound on parameter of alternative Beta distribution (default: 0.1).
#' @param alpha_max  Numeric in (0,1), upper bound on parameter of alternative Beta distribution (default: 0.9).

#' @return List with parameters of fitted conditional Beta-Uniform mixture
#' @export

gamma_glm_censored_em <- function(censored_Ps, Xs,  tau_censor, formula_rhs="~X1+X2",
                                  maxiter = 50, tau_pi0 = 0.5, alpha_min = 0.1,
                                  pi1_min = 0.01, pi1_max = 0.9, alpha_max=0.9){
  # basic transform
  m <- length(censored_Ps)
  censored_locs <- which(censored_Ps < tau_censor)
  uncensored_locs <- which(censored_Ps >= tau_censor)

  Ys <-  rep(-log(tau_censor/2), length(censored_Ps))
  Ys[uncensored_locs] <- -log(censored_Ps[uncensored_locs])


  ## Use Boca-Leek to initialize pi0s:
  BL_transform <- as.integer(censored_Ps >= tau_pi0)
  glm_bl <- glm(formula(paste("BL_transform", formula_rhs)), family=binomial(), data=Xs)
  pi1s_iter <- 1 - predict(glm_bl, type="response")/(1-tau_pi0)
  pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))

  # initialize EM Hs
  Hs_iter <- rep(NA, m)   #Our best guess for Hs is 1-FDR, use BH + our pilot pi0
  Ps_tmp <- censored_Ps
  Ps_tmp[censored_locs] <- tau_censor
  Ps_adj_tmp <- p.adjust(Ps_tmp, method="BH")
  Hs_iter <- 1 - Ps_adj_tmp*(1-pi1s_iter)

  for (i in 1:maxiter){
    # fit GLM on transformed and imputed responses
    glm_mus_iter <- glm( formula(paste("Ys", formula_rhs)), family=Gamma(), weights=Hs_iter, data=Xs)
    # update alphas and Hs
    alphas_iter <- predict(glm_mus_iter, type="link")
    alphas_iter <- pmax(alpha_min, pmin( alphas_iter, alpha_max))

    # E-step for Hs
    Hs_iter[uncensored_locs] <- 1-(1 - pi1s_iter[uncensored_locs])/((1 - pi1s_iter[uncensored_locs])  +
                                                                      pi1s_iter[uncensored_locs]*dbeta(censored_Ps[uncensored_locs], alphas_iter[uncensored_locs], 1))
    Hs_iter[censored_locs] <-1- (1 - pi1s_iter[censored_locs])*tau_censor/((1 - pi1s_iter[censored_locs])*tau_censor  +
                                                                             pi1s_iter[censored_locs]*pbeta(tau_censor, alphas_iter[censored_locs], 1))

    # E-step for censored Ys
    Ys[censored_locs] <- (1/alphas_iter[censored_locs] - log(tau_censor))#/Hs_iter[censored_locs]

    # fit logistic GLM
    glm_pi1s_iter <- glm(formula(paste("Hs_iter", formula_rhs)), family=quasibinomial(), data=Xs)
    
    # get logistic GLM predictions
    pi1s_iter <-  predict(glm_pi1s_iter, type="response")
    pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))


  }
  list(glm_mus = glm_mus_iter, glm_pi1s = glm_pi1s_iter, glm_bl = glm_bl,
       alphas = alphas_iter, pi1s = pi1s_iter)
}

#' Optimal weights based on a fitted conditional Beta-mixture model
#'
#' @param alpha    Nominal testing level
#' @param pi1s   Numeric vector with probability i-th hypothesis is an alternative
#' @param alphas  Numeric vector of first parameter of alternative Beta distribution (Beta(a_i, 1)).
#' @param numerator_bh  If TRUE (default is FALSE) then replace pi0 in numerator by 1.
#'
#' @return Numeric vector with weights.
#' @export
weights_betamix <- function(alpha, pi1s, alphas, numerator_bh=FALSE){
  wrapped_fun <- function(c) {
    ts <- get_thresholds_betamix(c, pi1s, alphas)
    get_tailfdr_betamix(ts, pi1s, alphas, numerator_bh = numerator_bh) - alpha
  }

  interval_endpts <- c(1e-10, 0.8)
  c_root <- stats::uniroot(wrapped_fun, interval_endpts)$root
  ts_root <- get_thresholds_betamix(c_root, pi1s, alphas)
  ws <- ts_root/sum(ts_root)*length(ts_root)
  ws
}




