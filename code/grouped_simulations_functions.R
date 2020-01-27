



apply_grouped_methods <- function(sim, alpha, gbh_only=FALSE){
  Xs <- sim$Xs #groups
  Ps <- sim$Ps
  Hs <- sim$Hs
  
  if (gbh_only){
    res <- bind_rows(mutate( fdp_eval(Hs,  gbh_simple(Ps, Xs, alpha)), method="GBH"),
                     mutate( fdp_eval(Hs,  ihw_gbh(Ps, Xs, alpha)), method="IHW-GBH"))
  } else {
    res <- bind_rows(mutate( fdp_eval(Hs,  p.adjust(Ps, method="BH") <= alpha), method="BH"),
                     mutate( fdp_eval(Hs,  gbh_simple(Ps, Xs, alpha)), method="GBH"),
                     mutate( fdp_eval(Hs,  gbh_simple(Ps, Xs, alpha, Storey=TRUE)), method="GBH-Storey"),
                     mutate( fdp_eval(Hs,  ihw_gbh(Ps, Xs, alpha)), method="IHW-GBH"),
                     mutate( fdp_eval(Hs,  ihw_gbh(Ps, Xs, alpha, Storey=TRUE)), method="IHW-GBH-Storey"),
                     mutate( fdp_eval(Hs,  ihw_nmeth_wrapper(Ps, Xs, alpha, pre_bin = FALSE)), method="IHW-Grenander"),
                     mutate( fdp_eval(Hs,  ihw_nmeth_wrapper(Ps, Xs, alpha, pre_bin = FALSE, Storey=TRUE)),
                             method="IHW-Grenander-Storey"),
                     mutate( fdp_eval(Hs,  stratified_clfdr(Ps, Xs, alpha)), method="Clfdr"),
                     mutate( fdp_eval(Hs,  groupwise_sabha(Ps, Xs, alpha)), method="SABHA"),
                     mutate( fdp_eval(Hs,  stratified_bhq(Ps, Xs, alpha)), method="SBH"))
  }
  res
}

#---------------------------------------------------
#------------------- Simulation settings------------
#---------------------------------------------------

# global null
grouped_global_null_sim <- function(m, K){
  m_div_K <- m/K
  Ps <- runif(m)
  Xs <- rep(1:K, m_div_K)
  Hs <- rep(0, m)
  data.frame(Hs=Hs, Ps=Ps, Xs=Xs)
}

eval_grouped_global_null_sim <- function(m, K, seed, alpha=0.1, gbh_only=TRUE){
  print(paste0("seed:",seed," and K:", K))
  sim <- grouped_global_null_sim(m, K)
  mutate(apply_grouped_methods(sim, alpha, gbh_only=gbh_only), m=m, K=K, seed=seed)
}

# global null, varying pi_0

grouped_sim <- function(m, K_coarse, pi0_global, sparsity_multiplier=2,  K=40){
  
  m_div_K <- m/K
  sparsity <- floor(K/sparsity_multiplier)
  # sparsity*pi + (K-sparsity) = K*pi_global => 
  # pi = (K*pi_global - K + sparsity)/sparsity
  #pi0_non_sparse <- K*(pi0_global-1)/sparsity + 1
  pi0_max <- 1
  pi0_left <- sparsity_multiplier * pi0_global - (sparsity_multiplier -1)
  # 2pi0_global <- pi0_min + 1
  pi0_min <- 2*pi0_left - 1
  pi0s_signal <- seq(from = pi0_min, to=pi0_max, length.out = sparsity)
  
  #non_sparse_groups <- sample(1:K, sparsity)
  pi0s <- rep(1, K)
  signal_idx <- seq(1, to=K, by=sparsity_multiplier)
  
  pi0s[signal_idx] <- pi0s_signal
  #pi0s[non_sparse_groups] <- pi0_non_sparse 
  eff_sizes <- rep(0, K)
  eff_sizes[signal_idx] <- seq(2.5, 0.5, length.out=sparsity)
  Hs <- rep(NA,m )
  Xs_tilde <- rep(1:K, m_div_K)
  Zs <- rep(NA,m)
  for (k in 1:K){
    grp_idx <- Xs_tilde == k
    Hs[grp_idx]   <- rbinom(m_div_K, 1, 1-pi0s[k])
    Zs[grp_idx] <- rnorm(m_div_K, Hs[grp_idx]*eff_sizes[k])
  }
  
  Xs <- ceiling(Xs_tilde*K_coarse/K)
  Ps <- 1-pnorm(Zs)
  return(list(signal_idx=signal_idx, pi0s=pi0s))
  data.frame(Hs=Hs, Ps=Ps, Xs=Xs)
}

eval_grouped_sim <- function(m, K_coarse, pi0_global, sparsity_multiplier, seed, alpha=0.1){
  print(paste0("seed:",seed," and K_Coarse:", K_coarse))
  sim <- grouped_sim(m, K_coarse, pi0_global, sparsity_multiplier)
  mutate(apply_grouped_methods(sim, alpha),
         m=m, K_coarse=K_coarse, seed=seed, pi0_global = pi0_global,
         sparsity_multiplier = sparsity_multiplier)
}

