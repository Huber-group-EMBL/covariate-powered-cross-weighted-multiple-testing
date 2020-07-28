#' FDR control procedure based on local FDRs
#'
#' @param unadj_p  Numeric vector of unadjusted p-values. This is only used to break ties.
#' @param lfdrs   Numeric vector with Local false discovery rate of each hypothesis
#' @param alpha    Significance level at which to apply method
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#' @export
oracle_local_fdr_test <- function(unadj_p, lfdrs, alpha){

  # now use the rejection rule described in Cai's paper

  # Remark:
  # When sorting lfdrs, we break ties by pvalues so that in the end within each stratum
  # we get monotonic adjusted p-values as a function of the p-values
  # This is mainly needed for grenander based lfdrs, with other
  # lfdr estimation methods lfdr ties are typically not a problem.

  o <- order(lfdrs, unadj_p)
  lfdrs_sorted <- lfdrs[o]
  fdr_estimate <- cumsum(lfdrs_sorted)/(1:length(unadj_p))
  adj_p <- rev(cummin(rev(fdr_estimate)))
  adj_p <- adj_p[order(o)]

  rjs <- adj_p <= alpha
  rjs
}
