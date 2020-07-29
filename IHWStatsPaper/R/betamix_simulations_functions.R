#' Simulation: Misspecified conditional Beta-uniform mixture model
#
#' @param m Number of hypotheses (default: m=10000)
#' @param mus_slope Numeric (default:1.5) parameter bar(beta) in equation (12)
#' @param one_sided_tests Bool (default:FALSE), if true adds some nulls that are strictly superuniform
#' @param prob_one_sided Numeric (default:0.25) proportion of nulls that are strictly superuniform
#'
#' @return Data frame with columns `Hs` (null or alternative), `Ps` (p-value), `Xs` (side-information),
#'         `alphas` (parameters of alternative distribution), `pi1s` (probabilities of being from the alternative distribution)
#' @export
beta_unif_sim <- function(m=10000, mus_slope=1.5, one_sided_tests=FALSE, prob_one_sided=0.25){
  Xs <- matrix(runif(m*2, 0,1), ncol=2)
  colnames(Xs) <- c("X1", "X2")

  pi1s <- ifelse( Xs[,1]^2 + Xs[,2]^2 <= 1, 0.02, 0.4)
  mus <- pmax(1.3, sqrt(Xs) %*% c(1,1)*mus_slope)

  mu_alphas <- 1/mus

  Hs <- stats::rbinom(m, size=1, prob=pi1s)
  Ps <- stats::runif(m)*(1-Hs) + stats::rbeta(m, mu_alphas, 1)*Hs
  Xs <- data.frame(Xs)
  if (one_sided_tests){
    Hs_alt <-  1- (1-Hs)*stats::rbinom(m, size=1, prob=prob_one_sided)
    Ps[Hs_alt == 0] <- stats::rbeta(sum(Hs_alt == 0), 1, 0.5)
  }
  list(Xs=Xs, Ps=Ps, Hs=Hs, alphas=mu_alphas, pi1s=pi1s)
}


error_fdp_table <- function(x) {
  if (inherits(x, "try-error")){
    x <- data.frame(rjs=NA, pow=NA, FDP=NA, FWER=NA)
  }
  x
}

eval_beta_unif_sim <- function(mu_slope, seed, one_sided_tests=FALSE, prob_one_sided = 0.25, alpha=0.1, m=10000, lfdr_only=FALSE ){
  print(seed)
  sim <- beta_unif_sim(m=m, mus_slope=mu_slope,
                       one_sided_tests=one_sided_tests, prob_one_sided = prob_one_sided)
  Xs <- sim$Xs
  Ps <- sim$Ps
  Hs <- sim$Hs
  mu_alphas <- sim$alphas
  pi1s <- sim$pi1s

  lfdr_oracle_res <- fdp_eval(Hs,  betamix_oracle_lfdr(Ps, pi1s, mu_alphas, alpha))
  lfdr_em_res <- error_fdp_table(try(fdp_eval(Hs,  betamix_datadriven_lfdr(Ps, Xs, alpha))))

  sim_res <-  bind_rows(mutate(lfdr_oracle_res, method="Lfdr-oracle"),
                        mutate(lfdr_em_res, method="Lfdr-EM"))
  if (!lfdr_only){
    bh_res <- fdp_eval(Hs,  p.adjust(Ps, method="BH") <= alpha)
    adapt_res <-  error_fdp_table(try(fdp_eval(Hs, adapt_mtp(Ps, Xs, alpha))))
    ihw_nmeth_res <- fdp_eval(Hs,   ihw_nmeth_wrapper(Ps,
                                                      interaction( cut(Xs[,1],5), cut(Xs[,2],5)), alpha))
    ihw_betamix_res <- error_fdp_table(try(fdp_eval(Hs,  ihw_betamix_censored(Ps, Xs, alpha, tau=0.1, kfolds=5, Storey=TRUE))))

    sim_res <- bind_rows(sim_res,
                         mutate(bh_res, method="BH"),
                         mutate(adapt_res, method="AdaPT"),
                         mutate(ihw_nmeth_res, method="IHW-Grenander"),
                         mutate(ihw_betamix_res, method="IHW-BetaMix"))
  }
  mutate(sim_res,
         seed = seed,
         pi0s = mean(1-Hs),
         mu_slope = mu_slope,
         alpha=alpha,
         m =m,
         one_sided_tests = one_sided_tests,
         prob_one_sided = prob_one_sided)
}
