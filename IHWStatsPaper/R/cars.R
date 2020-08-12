#' Modified CARS function
#'
#' Compared to the `CARS::CARS` function, there are three main differences:
#' 1) this function takes input that has already been summarized.
#' 2) it returns additional information (cf. below).
#' 3) we assume that we have Gaussian test-statistics with known variance.

#' @param t_1   Numeric vector of primary z-score statistics (converted to p-values via 2*(1-pnorm(abs(primary_stats)))).
#' @param t_2   Numeric vector of secondary statistic
#' @param alpha    Significance level at which to control FDR
#' @param tau   Threshold for choosing interesting locations for density estimation (default 0.9)
#' @param option Currently only supports 'regular' (and not sparse from the CARS package)
#'
#' @return A list that contains `de`, i.e., the vector with the significant tests and `cars_fun` which is a function that
#'         maps t_1,t_2 to the estimated local fdr.
#'
#' @importFrom np npudensbw
#' @importFrom np npudens
#' @export
CARS_extended <- function(t_1, t_2, alpha, tau=0.9, option='regular'){

  t_1.pval <- 2*stats::pnorm(-abs(t_1));
  m <- length(t_1)

  t_1.p.Est <- CARS::epsest.func(t_1,0,1);
  t_2.p.Est <- CARS::epsest.func(t_2,0,1);

  #Estimate the lfdrs
  t_1.density.Est <- stats::density(t_1,from=min(t_1)-10,to=max(t_1)+10,n=1000);
  t_1.density.Est <- CARS::lin.itp(t_1,t_1.density.Est$x,t_1.density.Est$y); #Ok
  t_1.Lfdr.Est <- (1-t_1.p.Est)*dnorm(t_1)/t_1.density.Est;
  t_1.Lfdr.Est[which(t_1.Lfdr.Est>1)] <- 1;

  t_2.density.Est <- stats::density(t_2,from=min(t_2)-10,to=max(t_2)+10,n=1000);
  t_2.density.Est <- CARS::lin.itp(t_2,t_2.density.Est$x,t_2.density.Est$y);
  t_2.Lfdr.Est <- (1-t_2.p.Est)*dnorm(t_2)/t_2.density.Est;
  t_2.Lfdr.Est[which(t_2.Lfdr.Est>1)] <- 1;
  S <- which(t_1.pval<=0.5);

  bandwidth <- npudensbw(~t_1[S]+t_2[S],bwmethod="normal-reference")$bw;
  hx <- bandwidth[1];
  ht <- bandwidth[2];
  if (option=='regular'){
    # Calculate the estimated CARS statistics based on bivariate density estimation for denominator
    cars.denominator.object <- npudens(~t_1+t_2,bws=bandwidth)
    cars.denominator <- cars.denominator.object$dens
  }

  #Estimate Correction
  sample_null <- stats::rnorm(m);
  t_1.den.Est <- stats::density(t_1,bw=hx,from=min(t_1)-10,to=max(t_1)+10,n=1000);
  sample_lfdr <- (1-t_1.p.Est)*stats::dnorm(sample_null)/CARS::lin.itp(sample_null,t_1.den.Est$x,t_1.den.Est$y);
  correction <- length(which(sample_lfdr>=tau))/m;


  T.tau <- which(t_1.Lfdr.Est>=tau);

  t_2.star.den.Est.object <- stats::density(t_2[T.tau],bw=ht,from=min(t_2[T.tau])-10,to=max(t_2[T.tau])+10,n=1000);
  t_2.star.den.Est <- CARS::lin.itp(t_2,t_2.star.den.Est.object$x, t_2.star.den.Est.object$y);

  cars.numerator <- stats::dnorm(t_1)*length(T.tau)/m*t_2.star.den.Est/correction;

  cars.Est <- cars.numerator/cars.denominator;
  cars.Est[which(cars.Est>=1)] <- 1;

  # function to compute CARS statistic on other data points
  if (option=='regular'){
    cars_fun <- function(t_1, t_2){
      new_df <- data.frame(t_1 = t_1, t_2 = t_2)
      tmp_denom <- predict(cars.denominator.object, newdata=new_df)
      tmp_t_2_star_den <- CARS::lin.itp(t_2, t_2.star.den.Est.object$x, t_2.star.den.Est.object$y);
      tmp_numer <- stats::dnorm(t_1)*length(T.tau)/m*tmp_t_2_star_den/correction;
      pmin(tmp_numer/tmp_denom,1)
    }
  }

  cars.sorted <- sort(cars.Est,decreasing=FALSE,index.return=TRUE);
  cars.sorted.cumsum <- cumsum(cars.sorted$x);

  decision <- rep(0,m);
  threshold <- 0;
  for (i in 1:m){
    if(cars.sorted.cumsum[i]/i <= alpha & cars.sorted.cumsum[i+1]/(i+1) > alpha){
      decision[cars.sorted$ix[1:i]] <- 1;
      threshold <- cars.Est[cars.sorted$ix[i]];
      break;
    }
  }

  res <- list(pvalue = t_1.pval, primary_statistic = t_1, ancillary_statistic = t_2 ,
              de=decision,cars=cars.Est,  th=threshold);
  if  (option=='regular'){
    res$denominator_density <- cars.denominator.object
    res$cars_fun <- cars_fun
  }

  return(res)

}




