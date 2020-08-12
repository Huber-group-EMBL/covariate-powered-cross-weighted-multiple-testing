#' Simulation: Two-sample t-test
#
#' @param m Number of hypotheses (default: m=10000)
#' @param k Number of alternatives (default: k=10)
#' @param n_x Number of samples for X (first group) in each test (default: n_x=50)
#' @param n_y Number of samples for Y (second group) in each test (default: n_y=50)
#' @param mu_x_1 Signal strength for first locations of X
#' @param mu_y_1 Signal strength for first locations of Y (default: 0)
#' @param mu_x_2 Signal strength for locations (k+1) to 2k of X
#' @param mu_y_2 Signal strength for locations (k+1) to 2k of Y
#' @param sd_x   Standard deviation of each measurement from X
#' @param sd_y   Standard deviation of each measurement from Y
#'
#' @return  List with entries `H` (0/1 vector with null or alternative), `x` (Matrix of dimension m * n_x with X data),
#'          `y` (Matrix with `Y` data) and `var_mat` (Matrix with 2 columns containing the measurement variance for each test * group)
#' @references The code here is a modification of the example code in the CARS package vignette.
#' @export
cars_ttest_data <- function(m=10000, k=10, n_x=50, n_y=50,
                            mu_x_1=3/sqrt(n_x),
                            mu_y_1=0,
                            mu_x_2 = 1/sqrt(n_x),
                            mu_y_2 = 1/sqrt(n_y),
                            sd_x=1,
                            sd_y=1){

  theta <- rep(0,m); # the vector for true states of nature

  # initialize matrices for x and y
  x <- matrix(rnorm(m*n_x,0,sd_x),m,n_x);
  y <- matrix(rnorm(m*n_y,0,sd_y),m,n_y);

  if (k > 0){
    theta[1:k] <- 1  # signal locations
    x[c(1:k),] <- x[c(1:k),] + mu_x_1;
    y[c(1:k),] <- y[c(1:k),] + mu_y_1;
    x[c((k+1):(2*k)),] <- x[c((k+1):(2*k)),] + mu_x_2;
    y[c((k+1):(2*k)),] <- y[c((k+1):(2*k)),] + mu_y_2;
  }
  var_mat <-   matrix(c(sd_x,sd_y),  ncol=2, nrow=m, byrow=TRUE)^2
  list(x=x, y=y, H=theta, var_mat=var_mat)
}



#' Data preprocessing for CARS to be used with p-value methods
#'
#' @param X,Y,variance X,Y matrices with data for two groups used for two-sample multiple testing and per-test per-sample variance matrix, as returned
#'                     by the `cars_ttest_data` function.
#' @return Data.frame with columns `pvalue` and `t_2` (ancillary --independent under the null-- testing statistic)
#' @export
#' @examples
#' cars_example_data <- cars_ttest_data(m=10000, k=2000, n_x=50, n_y=50, mu_x_1 = 0.5,
#                                     mu_x_2 = 0.25, mu_y_2 = 0.25)
#' cars_example_preprocess <- with(cars_example_data, CARS_preprocess(x, y, var_mat))
CARS_preprocess <- function(X,Y, variance){
  n_x <- ncol(X); m <- nrow(X);
  stopifnot(nrow(Y)==m)
  n_y <- ncol(Y);
  n <- n_x+n_y;

  #Calculate the mean of X, Y for each location
  X.mean <- base::rowMeans(X);
  Y.mean <- base::rowMeans(Y);

  #Assign corresponding variances for each location
  X.var <- variance[,1];
  Y.var <- variance[,2];

  #Calculate kappa and pooled variance
  kappa <- n_y*X.var/(n_x*Y.var);
  pool.sd <- sqrt(n_y/n*X.var+n_x/n*Y.var);

  #Calculate primary and auxiliary statistics
  t_1 <- (X.mean-Y.mean)/pool.sd*sqrt(n_x*n_y/n);
  t_2 <- (X.mean+kappa*Y.mean)/(sqrt(kappa)*pool.sd)*sqrt(n_x*n_y/n);


  t_1.pval <- 2*pnorm(-abs(t_1));
  deg <- (X.var/n_x+Y.var/n_y)^2/((X.var/n_x)^2/(n_x-1)+(Y.var/n_y)^2/(n_y-1));

  preprocessed_df <- data.frame( pvalue = t_1.pval, t_1 = t_1, t_2=t_2,
                                 X_var = variance[,1], Y_var = variance[,2],
                                 deg = deg)

  attr(preprocessed_df, "n_x") <- n_x
  attr(preprocessed_df, "n_y") <- n_y

  preprocessed_df
}


apply_cars_methods <- function(cars_data, alpha){
  Hs <- cars_data$H
  cars_prepr <- CARS_preprocess(cars_data$x, cars_data$y, cars_data$var_mat )

  # run methods
  cars_res <- CARS::CARS(cars_data$x, cars_data$y, alpha, tau=0.9, cars_data$var_mat, option="regular")$de
  cars_sparse_res <- CARS::CARS(cars_data$x, cars_data$y, alpha, tau=0.9, cars_data$var_mat, option="sparse")$de
  BH_res <- stats::p.adjust(cars_prepr$pvalue, method="BH") <= alpha

  ihw_cars_res <-  ihw_cars(cars_prepr$t_1,cars_prepr$t_2, alpha, Storey=TRUE)
  ihw_nmeth_res <- ihw_nmeth_wrapper(cars_prepr$pvalue,
                                     IHW::groups_by_filter(cars_prepr$t_2, 10),
                                     alpha, pre_bin=FALSE, Storey=TRUE) #with Storey

  # collect results
  res <- bind_rows(mutate( fdp_eval(Hs, BH_res), method="BH"),
                   mutate( fdp_eval(Hs, cars_res), method="CARS"),
                   mutate( fdp_eval(Hs, cars_sparse_res), method="CARS-sparse"),
                   mutate( fdp_eval(Hs, ihw_cars_res), method="IHW-Storey-CARS"),
                   mutate( fdp_eval(Hs, ihw_nmeth_res), method="IHW-Storey-Grenander"))

  res
}

#' Apply multiple testing methods to the simultaneous two-sample testing simulation.
#
#' @param m Number of hypotheses (default: m=10000)
#' @param k Number of alternatives (default: k=10)
#' @param mu_x_1 Signal strength for first locations of X
#' @param mu_x_2 Signal strength for locations (k+1) to 2k of X
#' @param seed     Integer; used for printing which simulation it running (does not set an actual RNG seed)
#' @param n_x Number of samples for X (first group) in each test (default: n_x=50)
#' @param n_y Number of samples for Y (second group) in each test (default: n_y=50)
#' @param alpha Numeric (default: 0.1), nominal significance level at which to apply methods
#'
#' @return Data frame with FDP and Power of different methods on this simulation.
#' @export
eval_cars_sim <- function(m, k, mu_x_1, mu_x_2 , seed, n_x=50, n_y=50, alpha=0.1){
  print(paste0("seed:",seed," and k:", k))

  cars_data <- cars_ttest_data(m=m, k=k, n_x=n_x, n_y=n_y, mu_x_1 = mu_x_1,
                               mu_x_2 = mu_x_2, mu_y_2 = mu_x_2)

  mutate(apply_cars_methods(cars_data, alpha),
         m = m, k = k, seed = seed, n_x = n_x, n_y = n_y,
         mu_x_1 = mu_x_1, mu_x_2 = mu_x_2, mu_y_2 = mu_x_2)
}

