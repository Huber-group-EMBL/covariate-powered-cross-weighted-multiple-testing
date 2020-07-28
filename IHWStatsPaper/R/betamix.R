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
#' @return Vector of Tail fdrs.
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



betamix_oracle_lfdr <- function(Ps, pi1s, mu_alphas, alpha){
  oracle_local_fdr_test(Ps,
                        get_localfdrs_betamix(Ps, pi1s, mu_alphas),
                        alpha)
}




#' Fit conditional beta-uniform mixture model using EM algorithm.
#'
#' @param Ps       Numeric vector of p-values
#' @param Xs       Data frame with features
#' @param formula_rhs Formula defining the RHS in the fitted GLMs, defaults to ~X1+X2 (used in simulations herein).
#' @param maxiter  Total number of iterations to run the EM algorithm
#' @param tau_pi0  Number in (0,1) used for the initialization of the EM algorithm through the Boca-Leek approach (default: 0.5).
#' @param pi1_min  Numeric in (0,1),  lower bound on conditionally probability of being an alternative (default: 0.01).
#' @param pi1_min  Numeric in (0,1), upper bound on conditionally probability of being a null (default: 0.9).
#' @param alpha_max  Numeric in (0,1), upper bound on parameter of alternative Beta distribution (default: 0.9).

#' @return List with parameters of fitted GLM
#' @export
gamma_glm_basic_em <- function(Ps, Xs, formula_rhs="~X1+X2", maxiter = 50, tau_pi0=0.5,
                               pi1_min = 0.01, pi1_max = 0.9, alpha_max=0.9){
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
    alphas_iter <- pmin( alphas_iter, alpha_max)
    Hs_iter <- 1-(1 - pi1s_iter)/((1 - pi1s_iter)  + pi1s_iter*dbeta(Ps, alphas_iter, 1))

    # fit logistic GLM
    glm_pi1s_iter <- glm(formula(paste("Hs_iter", formula_rhs)), family=quasibinomial(), data=Xs)
    # get logistic GLM predictions
    pi1s_iter <-  predict(glm_pi1s_iter, type="response")
    pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))
  }
  list(glm_mus = glm_mus_iter, glm_pi1s = glm_pi1s_iter, glm_bl = glm_bl,
       alphas = alphas_iter, pi1s = pi1s_iter)
}

