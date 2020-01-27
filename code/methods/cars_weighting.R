library(np)

CARS_extended <- function(t_1, t_2, alpha, tau=0.9, option=c('sparse','regular')){
  
  t_1.pval <- 2*pnorm(-abs(t_1));
  m <- length(t_1)
  
  t_1.p.Est <- epsest.func(t_1,0,1);
  t_2.p.Est <- epsest.func(t_2,0,1);
  
  #Estimate the lfdrs
  t_1.density.Est <- density(t_1,from=min(t_1)-10,to=max(t_1)+10,n=1000);
  t_1.density.Est <- lin.itp(t_1,t_1.density.Est$x,t_1.density.Est$y); #Ok
  t_1.Lfdr.Est <- (1-t_1.p.Est)*dnorm(t_1)/t_1.density.Est;
  t_1.Lfdr.Est[which(t_1.Lfdr.Est>1)] <- 1;
  
  t_2.density.Est <- density(t_2,from=min(t_2)-10,to=max(t_2)+10,n=1000);
  t_2.density.Est <- lin.itp(t_2,t_2.density.Est$x,t_2.density.Est$y);
  t_2.Lfdr.Est <- (1-t_2.p.Est)*dnorm(t_2)/t_2.density.Est;
  t_2.Lfdr.Est[which(t_2.Lfdr.Est>1)] <- 1;
  S <- which(t_1.pval<=0.5);
  
  bandwidth <- np::npudensbw(~t_1[S]+t_2[S],bwmethod="normal-reference")$bw;
  hx <- bandwidth[1];
  ht <- bandwidth[2];
  if (option=='regular'){
    # Calculate the estimated CARS statistics based on bivariate density estimation for denominator
    cars.denominator.object <- np::npudens(~t_1+t_2,bws=bandwidth)
    cars.denominator <- cars.denominator.object$dens
  }else if(option=='sparse'){
    # Calculate the estimated CARS statistics based on decomposed bivariate density estimation for denominator (null + non-null on t_2)
    t_2.pval <- 2*pnorm(-abs(t_2));
    S_t2 <- which(t_2.pval<=0.2);
    screened.den.Est <- density(t_1[S_t2],from=min(t_1[S_t2])-10,to=max(t_2[S_t2])+10,n=1000);
    screened.den.Est <- lin.itp(t_1,screened.den.Est$x,screened.den.Est$y);
    #Sparse option, decompose the denominator into two components
    cars.denominator <- (1-t_2.p.Est)*dnorm(t_2)*dnorm(t_1)+(1-t_2.Lfdr.Est)*t_2.density.Est*screened.den.Est;
  }
  
  #Estimate Correction 
  sample_null <- rnorm(m);
  t_1.den.Est <- density(t_1,bw=hx,from=min(t_1)-10,to=max(t_1)+10,n=1000);
  sample_lfdr <- (1-t_1.p.Est)*dnorm(sample_null)/lin.itp(sample_null,t_1.den.Est$x,t_1.den.Est$y);
  correction <- length(which(sample_lfdr>=tau))/m;  
  
  
  T.tau <- which(t_1.Lfdr.Est>=tau);
  
  t_2.star.den.Est.object <- density(t_2[T.tau],bw=ht,from=min(t_2[T.tau])-10,to=max(t_2[T.tau])+10,n=1000);
  t_2.star.den.Est <- lin.itp(t_2,t_2.star.den.Est.object$x, t_2.star.den.Est.object$y);
  
  cars.numerator <- dnorm(t_1)*length(T.tau)/m*t_2.star.den.Est/correction;
  
  cars.Est <- cars.numerator/cars.denominator;
  cars.Est[which(cars.Est>=1)] <- 1;
  
  # function to compute CARS statistic on other data points
  if (option=='regular'){
    cars_fun <- function(t_1, t_2){
      new_df <- data.frame(t_1 = t_1, t_2 = t_2) 
      tmp_denom <- predict(cars.denominator.object, newdata=new_df)
      tmp_t_2_star_den <- lin.itp(t_2, t_2.star.den.Est.object$x, t_2.star.den.Est.object$y);
      tmp_numer <- dnorm(t_1)*length(T.tau)/m*tmp_t_2_star_den/correction;
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



cars_weighter <- function(primary_stats, Xs, Xs_new, tau, alpha){
  cars_res <- CARS_extended(primary_stats, Xs, alpha, option="regular")
  Xs_quantiles <- quantile(Xs_new, seq(0,to=1, length.out=100))
  cars_get_threshold <- function(x, cars_fun, threshold){
    zs <- seq(1, 5, length=200)
    idx_min <- which.min(abs( cars_fun(zs,x) - threshold))
    if (idx_min == 1){
      return(4)
    } else {
      return(zs[idx_min])
    }
  }
  ts <- sapply(Xs_quantiles, function(x)  cars_get_threshold(x, cars_res$cars_fun, cars_res$th))
  ws <- approx(x=Xs_quantiles,y=2*(1-pnorm(ts)), xout=Xs_new)$y
  ws <- ws*length(ws)/sum(ws)
  ws
}


ihw_cars <- function(primary_stats, Xs, alpha, ...){
  ihw_bh(primary_stats, Xs, alpha, cars_weighter, tau=1, stat_type="zscore",...)
}

