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

#' stratified_bhq: Stratified Benjamini Hochberg
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#'
#' @references Sun, Lei, et al. "Stratified false discovery control for large-scale hypothesis testing with application to genome-wide
#'    association studies." Genetic epidemiology 30.6 (2006): 519-530.
#' @references Also see similar function in IHWpaper package on Bioconductor.
#' @export
stratified_bhq <- function(unadj_p, groups, alpha){
  groups <- as.factor(groups)
  pv_list <- split(unadj_p, groups)

  m        <- length(unadj_p)
  m_groups <- sapply(pv_list, length)

  adj_pv_list <- lapply(pv_list, function(pv) p.adjust(pv, method="BH"))
  adj_p <- unsplit(adj_pv_list, groups)
  adj_p <= alpha
}


#' stratified clfdr: local fdr based method
#'
#' The two-groups model is estimated using fdrtool.
#'
#' @param unadj_p  Numeric vector of unadjusted p-values.
#' @param groups   Factor to which different hypotheses belong
#' @param alpha    Significance level at which to apply method
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#'
#' @references Cai, T. Tony, and Wenguang Sun. "Simultaneous testing of grouped hypotheses: Finding needles in multiple haystacks."
#'           Journal of the American Statistical Association 104.488 (2009).
#' @references Also see similar function in IHWpaper package on Bioconductor.
#' @export
stratified_clfdr <- function(unadj_p, groups, alpha){
  groups <- as.factor(groups)
  lfdr_fun <- function(pv) {fdrtool::fdrtool(pv, statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr}
  # estimate local fdr within each stratum first
  pvals_list <- split(unadj_p, groups)
  lfdr_list <- lapply(pvals_list, lfdr_fun)
  lfdrs <- unsplit(lfdr_list, groups)

  # now use the rejection rule described in Cai and Sun [2009] paper
  oracle_local_fdr_test(unadj_p, lfdrs, alpha)
}
