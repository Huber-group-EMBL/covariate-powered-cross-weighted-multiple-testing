#' Wrapper for the IHW procedure
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Xs   Vector or matrix of covariates
#' @param alpha    Significance level at which to apply method
#' @param Storey   Bool (default: FALSE): is the procedure pi0 adaptive or not?
#' @param pre_bin  Keyword argument that defines the binning procedure
#'
#' @return         Binary vector of rejected/non-rejected hypotheses.
#'
#' @export
ihw_nmeth_wrapper <- function(Ps, Xs, alpha, pre_bin = FALSE, Storey=FALSE){
  if (pre_bin == "2D"){
    binning_idx <- interaction( cut(Xs[,1],5), cut(Xs[,2],5))
  } else if (pre_bin == "1D"){
    binning_idx <- IHW::groups_by_filter(Xs, 10)
  } else {
    binning_idx <- as.factor(Xs)
  }
  if (Storey){
    ihw_nmeth_fit <- ihw(Ps, binning_idx, alpha, lambdas=Inf,
                         null_proportion=TRUE, null_proportion_level=0.5)
  } else {
    ihw_nmeth_fit <- ihw(Ps, binning_idx, alpha, lambdas=Inf)
  }
  rejected_hypotheses(ihw_nmeth_fit)
}


#' The tau-weighted BH multiple testing procedure
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param ws   Numeric vector of multiple testing weights
#' @param tau  Numeric (default = 0.5), the level at which tau-censoring is applied.
#'
#' @return Vector of adjusted p-values
#' @export
tau_weighted_bh <- function(Ps, ws, tau=0.5){
  weighted_pvals <- ifelse(ws == 0, 1, Ps/ws)
  weighted_pvals[weighted_pvals > 1] <- 1
  weighted_pvals[ Ps >= tau] <- 1
  adj_p <- stats::p.adjust(weighted_pvals, method="BH")
  adj_p
}

#' Storey's pi0 estimator for weighted multiple testinng
#'
#' @param pvalues   Numeric vector of unadjusted p-values.
#' @param weights   Numeric vector of multiple testing weights
#' @param tau       Numeric (default = 0.5), the level at which tau-censoring is applied.
#' @param m         Total number of tests (default: `length(pvalues)`)
#'
#' @return Estimated null proportion
#' @export
weighted_storey_pi0 <- function(pvalues, weights, tau=0.5, m = length(pvalues)){
  w_max <- max(weights)
  num <- w_max + sum( weights * (pvalues > tau))
  num/m/(1-tau)
}

#' Main wrapper for general IHW-BH/Storey procedure
#'
#' @param primary_stat  Vector of test statistics, typically pvalues (see `stat_type`)
#' @param Xs       Side-information for each test
#' @param alpha    Significance level at which to apply method
#' @param wt_fitter Function that does the weight learning, for example `grouped_storey_weighter`.
#' @param tau      Numeric (default: 0.5) censoring level
#' @param Storey   Bool (default: FALSE): is the procedure pi0 adaptive or not?
#' @param folds    Factor with assignments of tests to folds, defaults to NULL in which case they are chosen at random.
#' @param kfolds   Integer, number of folds into which to split hypotheses (ignored if `folds` is not NULL).
#' @param stat_type Either "pvalue" (default) or "zscore"
#' @param return_weights Boolean, whether to return weigts (defauls to FALSE)
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#'
#' @export
ihw_bh <- function(primary_stat, Xs, alpha, wt_fitter, tau=0.5,
                   folds=NULL, kfolds=5L, Storey=FALSE, stat_type="pvalue",
                   return_weights=FALSE,...){
  if (stat_type == "pvalue"){
    Ps <- primary_stat
  } else if (stat_type == "zscore") {
    Ps <- 2*(1-pnorm(abs(primary_stat)))
  }
  m <- length(Ps)
  ws <- rep(NA, m)
  if (is.null(folds)){
    folds <- sample(1:kfolds, m, replace=TRUE)
  } else {
    kfolds <- length(unique(folds))
  }
  st_pi0s <- rep(NA, kfolds)
  for (i in 1:kfolds){
    test_idx <- folds == i
    train_idx <- folds != i
    if (is.vector(Xs) ){
      Xs_train <- Xs[train_idx]
      Xs_test <- Xs[test_idx]
    } else {
      Xs_train <- Xs[train_idx,]
      Xs_test <- Xs[test_idx,]
    }
    ws_tmp <- wt_fitter(primary_stat[train_idx], Xs_train, Xs_test, tau, alpha, ...)
    #ws_tmp <- pmin(10, ws_tmp)
    if (all(ws_tmp == 0)){
      ws_tmp <- rep(1, length(ws_tmp))
    } else {
      ws_tmp <- length(ws_tmp)*ws_tmp/sum(ws_tmp)
    }
    st_pi0s[i] <- weighted_storey_pi0(Ps[test_idx], ws_tmp)
    if (Storey){
      ws[test_idx] = ws_tmp/st_pi0s[i]
    } else {
      ws[test_idx] = ws_tmp
    }
  }

  adj_p <- tau_weighted_bh(Ps, ws, tau=tau)
  rjs <- adj_p <= alpha

  if (return_weights){
    return(list(ws=ws, folds=folds, adj_p=adj_p, rjs=rjs, storey_pi0s=st_pi0s))
  } else {
    return(rjs)
  }
}


# Specific instantiations of GBH through specific choice of weighting function

#' Learn multiple testing weights based on the GBH procedure
#'
#' @param Ps        Numeric vector of unadjusted p-values in training set.
#' @param Xs        Factor of groups in training set.
#' @param Xs_new    Factor of groups in test-set (out-of-fold).
#' @param tau       Numeric (default = 0.5), the level at which tau-censoring is applied.
#' @param alpha    Significance level at which to apply method (not used for this weight learning method)
#'
#' @return Numeric vector with weights of test-set hypotheses
#' @export
grouped_storey_weighter <- function(Ps, Xs, Xs_new, tau, alpha){
  both_lvls <- unique(c(Xs,Xs_new))
  train_groups <- factor(Xs, levels = both_lvls)
  test_groups <- factor(Xs_new, levels= both_lvls)

  m_test <- length(test_groups)

  pv_list <- split(Ps, train_groups)
  test_xs_list <-  split(test_groups, test_groups)

  pi0_fun <- function(pv) {
    if (length(pv) == 0){
      pi0 <- 1
    } else {
      pi0 <- min(1, (1+sum( pv >= tau))/length(pv)/(1-tau))
    }
    pi0
  }
  pi0_groups <- sapply(pv_list, pi0_fun)

  # in this case we set all wts to 0
  if (all(pi0_groups ==1 )){
    ws <- rep(1,m_test)
    # do actual work
  } else {
    ws <- unsplit( mapply(function(pi0_g, xs) {(1-pi0_g)/pi0_g}, pi0_groups, test_xs_list,  SIMPLIFY=FALSE),
                   test_groups)
    ws <- ws/sum(ws)*m_test
  }
  ws
}

#' The IHW-GBH multiple testing procedure
#'
#' @param Ps       Numeric vector of unadjusted p-values.
#' @param Xs       Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#' @param tau      Numeric (default: 0.5) level at which to apply Storey's pi0 estimator
#' @param Storey   Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return         Binary vector of rejected/non-rejected hypotheses.
#' @export
ihw_gbh <- function(Ps,Xs, alpha, tau=0.5, Storey=FALSE, ...){
  ihw_bh(Ps,Xs, alpha, grouped_storey_weighter, tau=tau, Storey=Storey, ...)
}

#' Learn multiple testing weights based on EM fitting of the betamix model and local FDRs
#'
#' @param Ps        Numeric vector of unadjusted p-values in training set.
#' @param Xs        Factor of groups in training set.
#' @param Xs_new    Factor of groups in test-set (out-of-fold).
#' @param tau       Numeric (default = 0.5), the level at which tau-censoring is applied.
#' @param alpha    Significance level at which to apply method (not used for this weight learning method)
#' @param formula_rhs Formula defining the RHS in the fitted GLMs, defaults to ~X1+X2 (used in simulations herein).
#' @param maxiter  (Default: 200) Total number of iterations to run the EM algorithm
#' @param pi1_min  Numeric in (0,1),  lower bound on conditionally probability of being an alternative (default: 0.01).
#' @param pi1_max  Numeric in (0,1), upper bound on conditionally probability of being a null (default: 0.9).
#' @param alpha_max  Numeric in (0,1), upper bound on parameter of alternative Beta distribution (default: 0.9).
#' @param numerator_bh Bool (default:TRUE) Define the marginal FDR as sum pi_0(X_i) t(X_i), if FALSE, or sum t(X_i) if TRUE
#' @return Numeric vector with weights of test-set hypotheses
#' @export
censored_betamix_weighter <- function(Ps, Xs, Xs_new, tau, alpha,
                                      formula_rhs="~X1+X2",
                                      pi1_min = 0.01,
                                      pi1_max = 0.9,
                                      alpha_max = 0.9,
                                      maxiter = 200,
                                      numerator_bh = TRUE,
                                      ...){
  Ys <- Ps*(Ps >= tau)
  gamma_glm_censored_fit <- gamma_glm_censored_em(Ys, Xs, tau, maxiter = maxiter, formula_rhs=formula_rhs, tau_pi0=0.5,
                                                  pi1_min = pi1_min,
                                                  pi1_max = pi1_max,
                                                  alpha_max = alpha_max, ...)
  alphas_new <- predict(gamma_glm_censored_fit$glm_mus, newdata=Xs_new, type="link")
  pi1s_new <- predict(gamma_glm_censored_fit$glm_pi1s, newdata=Xs_new, type="response")

  alphas_new <- pmin(alphas_new, alpha_max)
  pi1s_new <- pmax(pi1_min, pmin(pi1_max, pi1s_new))

  weights_betamix(alpha, pi1s_new, alphas_new, numerator_bh = numerator_bh)
}


#' The IHW-BH/Storey-Betamix multiple testing procedure
#'
#' @param Ps       Numeric vector of unadjusted p-values.
#' @param Xs       Numeric vector or matrix with covariates
#' @param alpha    Significance level at which to apply method
#' @param tau      Numeric (default: 0.5) level at which to apply Storey's pi0 estimator
#' @param Storey   Bool (default: FALSE): is the procedure pi0 adaptive or not?
#' @param numerator_bh Bool (default:TRUE) Define the marginal FDR as sum pi_0(X_i) t(X_i), if FALSE, or sum t(X_i) if TRUE
#' @param ...      Other arguments passed to `censored_betamix_weighter`
#'
#' @return         Binary vector of rejected/non-rejected hypotheses.
#' @export
ihw_betamix_censored <- function(Ps,Xs, alpha, tau=0.1, Storey=FALSE, numerator_bh=TRUE,...){
  ihw_bh(Ps,Xs, alpha, censored_betamix_weighter, tau=tau, Storey=Storey, numerator_bh=numerator_bh,...)
}

