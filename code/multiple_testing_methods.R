source("beta_uniform_glm.R")

oracle_local_fdr_test <- function(unadj_p, lfdrs, alpha){

    # now use the rejection rule described in Cai's paper
    
    # Remark:
    # When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
    # we get monotonic adjusted p-values as a function of the p-values
    # This is mainly needed for grenander based lfdrs, with most other
    # lfdr estimation methods lfdr ties are not a problem usually
    
    o <- order(lfdrs, unadj_p)
    lfdrs_sorted <- lfdrs[o]
    fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
    adj_p <- rev(cummin(rev(fdr_estimate)))
    adj_p <- adj_p[order(o)]
    
    rjs <- adj_p <= alpha
    rjs
}

betamix_oracle_lfdr <- function(Ps, pi1s, mu_alphas, alpha){
  oracle_local_fdr_test(Ps,
                      get_localfdrs_betamix(Ps, pi1s, mu_alphas), 
                      alpha)
}


betamix_datadriven_lfdr <- function(Ps, Xs, alpha){
  gamma_glm_fit  <- gamma_glm_basic_em(Ps, Xs, maxiter = 200, tau_pi0=0.5)
  betamix_oracle_lfdr(Ps, gamma_glm_fit$pi1s, gamma_glm_fit$alphas, alpha)
}

# AdaPT
library(adaptMT)

adapt_mtp <- function(Ps, Xs, alpha){
  adapt_glm_fit <- adapt_glm(as.data.frame(Xs), Ps, "~X1+X2", "~X1+X2", alphas=alpha)
  adapt_glm_rjs <- adapt_glm_fit$qvals <= alpha
  adapt_glm_rjs
}


# IHW etc
library(IHW)

ihw_nmeth_wrapper <- function(Ps, Xs, alpha, pre_bin = FALSE, Storey=FALSE){
  if (pre_bin == "2D"){
    binning_idx <- interaction( cut(Xs[,1],5), cut(Xs[,2],5))
  } else if (pre_bin == "1D"){
    binning_idx <- IHW::groups_by_filter(Xs, 10)
  } else {
    binning_idx <- Xs
  }
  if (Storey){
    ihw_nmeth_fit <- ihw(Ps, binning_idx, alpha, lambdas=Inf,
                        null_proportion=TRUE, null_proportion_level=0.5)
  } else {
    ihw_nmeth_fit <- ihw(Ps, binning_idx, alpha, lambdas=Inf)
  }
  rejected_hypotheses(ihw_nmeth_fit)
}

tau_weighted_bh <- function(Ps, ws, tau=0.5){
  weighted_pvals <- ifelse(ws == 0, 1, Ps/ws)
  weighted_pvals[weighted_pvals > 1] <- 1
  weighted_pvals[ Ps >= tau] <- 1
  adj_p <- p.adjust(weighted_pvals, method="BH")
  adj_p
}

weighted_storey_pi0 <- function(pvalues, weights, tau=0.5, m = length(pvalues)){
  w_max <- max(weights)
  num <- w_max + sum( weights * (pvalues > tau))
  num/m/(1-tau)
}

ihw_bh <- function(primary_stat,Xs, alpha, wt_fitter, tau=0.5, 
                   folds=NULL, kfolds=5L, Storey=FALSE, stat_type="pvalue",
                   return_weights=FALSE){
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
    ws_tmp <- wt_fitter(primary_stat[train_idx], Xs_train, Xs_test, tau, alpha)
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

censored_betamix_weighter <- function(Ps, Xs, Xs_new, tau, alpha, 
                                      pi1_min = 0.01,
                                      pi1_max = 0.9,
                                      alpha_max = 0.9){
  Ys <- Ps*(Ps >= tau)
  gamma_glm_censored_fit <- gamma_glm_censored_em(Ys, Xs, tau, maxiter = 200, tau_pi0=0.5)
  alphas_new <- predict(gamma_glm_censored_fit$glm_mus, newdata=Xs_new, type="link")
  pi1s_new <- predict(gamma_glm_censored_fit$glm_pi1s, newdata=Xs_new, type="response")
  
  alphas_new <- pmin(alphas_new, alpha_max)
  pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_new))
  
  weights_betamix(alpha, pi1s_new, alphas_new, numerator_bh = TRUE)
}

uncensored_betamix_weighter <- function(Ps, Xs, Xs_new, tau, alpha,
                                        pi1_min = 0.01,
                                        pi1_max = 0.9,
                                        alpha_max = 0.9){
  gamma_glm_fit <- gamma_glm_basic_em(Ps, Xs, maxiter = 200, tau_pi0=0.5)
  alphas_new <- predict(gamma_glm_fit$glm_mus, newdata=Xs_new, type="link")
  pi1s_new <- predict(gamma_glm_fit$glm_pi1s, newdata=Xs_new, type="response")
  
  alphas_new <- pmin(alphas_new, alpha_max)
  pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_new))
  
  weights_betamix(alpha, pi1s_new, alphas_new, numerator_bh = TRUE)
}

ihw_betamix_censored <- function(Ps,Xs, alpha, tau=0.1, ...){
  ihw_bh(Ps,Xs, alpha, censored_betamix_weighter, tau=tau, ...)
}

ihw_betamix_uncensored <- function(Ps,Xs, alpha, ...){
  ihw_bh(Ps,Xs, alpha, uncensored_betamix_weighter, tau=1, ...)
}

#-----------------------------------------
#------ grouped procedures----------------
#-----------------------------------------

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
  #print(pi0_groups)
 # print(str(test_groups))
 # print(str(train_groups))
  
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

ihw_gbh <- function(Ps,Xs, alpha,tau=0.5, ...){
  ihw_bh(Ps,Xs, alpha, grouped_storey_weighter, tau=tau, ...)
}

gbh_simple <- function(unadj_p, groups, alpha, tau=0.5, Storey=FALSE){
  
  groups <- as.factor(groups)
  pv_list <- split(unadj_p, groups)
  
  pi0_fun <- function(pv) {min(1, (1+sum( pv >= tau))/length(pv)/(1-tau))}

  m        <- length(unadj_p)
  m_groups <- sapply(pv_list, length)
  pi0_groups <- sapply(pv_list, pi0_fun)
  
  # in this case we set all wts to 0
  if (all(pi0_groups ==1 )){
    ws <- rep(1,m)
    # do actual work
  } else if  (all(pi0_groups ==0 )){
    ws <- rep(1,m)
  } else {  
    
    ws <- unsplit( mapply(function(pi0_g, pv) {(1-pi0_g)/pi0_g}, pi0_groups, pv_list,  SIMPLIFY=FALSE), 
                    groups)
    ws <- ws/sum(ws)*m
  }
  if (Storey){
    pi0 <- weighted_storey_pi0(unadj_p, ws, tau=tau)
    unadj_p <- unadj_p * pi0
    
  }
  adj_p <- tau_weighted_bh(unadj_p, ws, tau=tau)
  rjs <- adj_p <= alpha
  #return(list(ws=ws, adj_p=adj_p, rjs=rjs))
  
  rjs
}


stratified_bhq <- function(unadj_p, groups, alpha){
  groups <- as.factor(groups)
  pv_list <- split(unadj_p, groups)
  
  m        <- length(unadj_p)
  m_groups <- sapply(pv_list, length)
  
  adj_pv_list <- lapply(pv_list, function(pv) p.adjust(pv, method="BH"))
  adj_p <- unsplit(adj_pv_list, groups)
  adj_p <= alpha
}
  
library(fdrtool)
stratified_clfdr <- function(unadj_p, groups, alpha){
  lfdr_fun <- function(pv) {fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr}
  # estimate local fdr within each stratum first
  pvals_list <- split(unadj_p, groups)
  lfdr_list <- lapply(pvals_list, lfdr_fun)
  lfdrs <- unsplit(lfdr_list, groups)
  
  # now use the rejection rule described in Cai's paper
  oracle_local_fdr_test(unadj_p, lfdrs, alpha)
    
}

#------------------------------------------
## SABHA
#------------------------------------------
source("groupwise_sabha_pi0.R")

##--------------------------------------------
# CARS and IHW-CARS
#---------------------------------------------
library(CARS)
source("cars_weighting.R")

# Evaluation

fdp_eval <- function(Hs, rjs){
  rjs_total <- sum(rjs)
  pow <- sum(rjs*Hs)/max(1,sum(Hs))
  FDP <- sum(rjs*(1-Hs))/max(1,rjs_total)
  FWER <- sum( (1-Hs)*rjs) > 0
  data.frame(rjs=rjs_total, pow=pow, FDP=FDP, FWER=FWER)
}




