  #----------------------------------------------------------------
  #  Data generation: Simulation function
  #----------------------------------------------------------------
  

cars_ttest_data <- function(m=10000, k=10, n_x=50, n_y=50, 
                            mu_x_1=3/sqrt(n_x),  # signal strength for first locations of x
                            mu_y_1=0,  # signal strength for first locations of y
                            mu_x_2 = 1/sqrt(n_x), # signal strength for location 501 to 1000 of x
                            mu_y_2 = 1/sqrt(n_y), # signal strength for location 501 to 1000 of y
                            sd_x=1,
                            sd_y=1){ 
  
  theta <- rep(0,m); # the vector for true states of nature
  
  # initialize matrices for x and y
  x <- matrix(rnorm(m*n_x,0,sd_x),m,n_x);
  y <- matrix(rnorm(m*n_y,0,sd_y),m,n_y);
  
  if (k > 0){
    theta[1:k] <- 1; # signal locations: first 500
    # fill in data matrices according to Setting 1 (k=500)
    x[c(1:k),] <- x[c(1:k),] + mu_x_1;
    y[c(1:k),] <- y[c(1:k),] + mu_y_1;
    x[c((k+1):(2*k)),] <- x[c((k+1):(2*k)),] + mu_x_2;
    y[c((k+1):(2*k)),] <- y[c((k+1):(2*k)),] + mu_y_2;
  }
  var_mat <-   matrix(c(sd_x,sd_y),  ncol=2, nrow=m, byrow=TRUE)^2
  list(x=x, y=y, H=theta, var_mat=var_mat)
}

#----------------------------------------------------------------
#  Data preprocessing for CARS to be used with p-value methods
#----------------------------------------------------------------

CARS_preprocess <- function(X,Y, variance){
#Validate input types.
  if (is.data.frame(X)){
    X.names=names(X)
    X = as.matrix(X,rownames.force=F)
  } else if (is.matrix(X))
    X.names=colnames(X)
  else
    stop('Input X must be a matrix or data frame')

  if (is.data.frame(Y)){
    Y.names=names(Y)
    Y = as.matrix(Y,rownames.force=F)
  } else if (is.matrix(Y))
    Y.names=colnames(Y)
  else
    stop('Input Y must be a matrix or data frame')

  #validate input dimensions.
  n_x <- ncol(X); m <- nrow(X);
  stopifnot(nrow(Y)==m)
  n_y <- ncol(Y);
  n <- n_x+n_y;

  #Calculate the mean of X, Y for each location
  X.mean <- rowMeans(X);
  Y.mean <- rowMeans(Y);

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
  cars_res <- CARS(cars_data$x, cars_data$y, alpha, tau=0.9, cars_data$var_mat, option="regular")$de
  cars_sparse_res <- CARS(cars_data$x, cars_data$y, alpha, tau=0.9, cars_data$var_mat, option="sparse")$de
  BH_res <- p.adjust(cars_prepr$pvalue, method="BH") <= alpha
  
  ihw_cars_res <- ihw_cars(cars_prepr$t_1,cars_prepr$t_2, alpha, Storey=TRUE)
  ihw_nmeth_res <- ihw_nmeth_wrapper(cars_prepr$pvalue,IHW::
                                       groups_by_filter(cars_prepr$t_2, 10),
                                     alpha, pre_bin=FALSE, Storey=TRUE) #with Storey
  
  # collect results
  res <- bind_rows(mutate( fdp_eval(Hs, BH_res), method="BH"),
                     mutate( fdp_eval(Hs,  cars_res), method="CARS"),
                     mutate( fdp_eval(Hs,  cars_sparse_res), method="CARS-sparse"),
                     mutate( fdp_eval(Hs,  ihw_cars_res), method="IHW-CARS"),
                     mutate( fdp_eval(Hs,  ihw_nmeth_res), method="IHW-Grenander"))
  
  res
}


eval_cars_sim <- function(m, k, mu_x_1, mu_x_2 , seed, n_x=50, n_y=50, alpha=0.1){
  print(paste0("seed:",seed," and k:", k))
  
  cars_data <- cars_ttest_data(m=m, k=k, n_x=n_x, n_y=n_y, mu_x_1 = mu_x_1, 
                               mu_x_2 = mu_x_2, mu_y_2 = mu_x_2)
  
  mutate(apply_cars_methods(cars_data, alpha),
         m = m, k = k, seed = seed, n_x = n_x, n_y = n_y,
         mu_x_1 = mu_x_1, mu_x_2 = mu_x_2)
}

