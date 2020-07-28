#' Wrapper of AdaPT multiple testing procedure
#'
#' @param Ps       Numeric vector of p-values
#' @param Xs       Data frame with features
#' @param formula_rhs Formula defining the RHS in the fitted GLMs, defaults to ~X1+X2 (used in simulations herein).
#'
#' @param alpha    Nominal testing level
#' @param return_fit  Boolean, whether to return the fitted adapt object (or only the indicator of rejections), defaults to false
#' @return Binary vector of rejected/non-rejected hypotheses.
#'
#' @references AdaptMT CRAN package
#' @export
adapt_mtp <- function(Ps, Xs, alpha, formula_rhs="~X1+X2", return_fit=FALSE){
  adapt_glm_fit <- adaptMT::adapt_glm(as.data.frame(Xs), Ps, formula_rhs, formula_rhs, alphas=alpha)
  adapt_glm_rjs <- adaptMT::adapt_glm_fit$qvals <= alpha
  if (return_fit){
    return(list(rjs=adapt_glm_rjs, fit=adapt_glm_fit))
  } else {
    return(adapt_glm_rjs)
  }
}
