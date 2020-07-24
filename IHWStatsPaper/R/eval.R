#' Evaluate multiple testing procedure
#
#' @param Hs Vector with indicators of alternatives (1) and true nulls (0)
#' @param rjs Vector with indicator of rejected hypotheses
#' @return Data frame with columns `rjs` (total rejections), `pow` (Power), `FDP` (False discovery proportion)
#' @export
fdp_eval <- function(Hs, rjs){
  rjs_total <- sum(rjs)
  pow <- sum(rjs*Hs)/max(1,sum(Hs))
  FDP <- sum(rjs*(1-Hs))/max(1,rjs_total)
  FWER <- sum( (1-Hs)*rjs) > 0
  data.frame(rjs=rjs_total, pow=pow, FDP=FDP, FWER=FWER)
}
