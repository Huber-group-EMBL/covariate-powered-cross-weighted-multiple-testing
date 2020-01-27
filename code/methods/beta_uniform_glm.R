gamma_glm_basic_em <- function(Ps, Xs, maxiter = 50, tau_pi0=0.5){
  # basic transform
  Ys <- -log(Ps)
  m <- length(Ys)
  pi1_min <- 0.01
  pi1_max <- 0.9
  alpha_max <- 0.9
  # initialize EM iterations

  ## Use Boca-Leek to initialize pi0s:
  BL_transform <- as.integer(Ps >= tau_pi0)
  glm_bl <- glm(BL_transform ~ X1 + X2, family=binomial(), data=Xs)
  pi1s_iter <- 1 - predict(glm_bl, type="response")/(1-tau_pi0)
  pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))
  
  # initial guess for HS is 1-FDR, use BH + our pilot pi0
  Hs_iter <- rep(NA, m)
  Ps_adj_tmp <- p.adjust(Ps, method="BH")
  Hs_iter <- 1 - Ps_adj_tmp*(1-pi1s_iter)
  
  for (i in 1:maxiter){
    
    # fit GLM on transformed responses
    glm_mus_iter <- glm( Ys ~ X1 + X2, family=Gamma(), weights=Hs_iter, data=Xs )
    # update alphas and Hs
    alphas_iter <- predict(glm_mus_iter, type="link")
    alphas_iter <- pmin( alphas_iter, alpha_max)
    Hs_iter <- 1-(1 - pi1s_iter)/((1 - pi1s_iter)  + pi1s_iter*dbeta(Ps, alphas_iter, 1))
    
    # fit logistic GLM
    glm_pi1s_iter <- glm(Hs_iter ~ X1 + X2, family=quasibinomial(), data=Xs)
    # get logistic GLM predictions
    pi1s_iter <-  predict(glm_pi1s_iter, type="response")
    pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))
  }
  list(glm_mus = glm_mus_iter, glm_pi1s = glm_pi1s_iter, glm_bl = glm_bl, 
       alphas = alphas_iter, pi1s = pi1s_iter)
}


gamma_glm_censored_em <- function(censored_Ps, Xs,  tau_censor, maxiter = 50, tau_pi0=0.5){
  # basic transform
  m <- length(censored_Ps)
  censored_locs <- which(censored_Ps == 0)
  uncensored_locs <- which(censored_Ps >= tau_censor)
  #Ys <- -log(runif(length(censored_Ps) , min=0, max=tau_censor)) # initialize like this
  Ys <-  rep(-log(tau_censor/2), length(censored_Ps))
  Ys[uncensored_locs] <- -log(censored_Ps[uncensored_locs])
  
  pi1_min <- 0.01
  pi1_max <- 0.9
  alpha_max <- 0.9
  
 
  ## Use Boca-Leek to initialize pi0s:
  BL_transform <- as.integer(censored_Ps >= tau_pi0)
  glm_bl <- glm(BL_transform ~ X1+X2, family=binomial(), data=Xs)
  pi1s_iter <- 1 - predict(glm_bl, type="response")/(1-tau_pi0)
  pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))
  
  # initialize EM Hs
  ## Hs: just set to 0.5
  Hs_iter <- rep(NA, m)   #Our best guess for Hs is 1-FDR, use BH + our pilot pi0
  Ps_tmp <- censored_Ps
  Ps_tmp[censored_locs] <- tau_censor
  Ps_adj_tmp <- p.adjust(Ps_tmp, method="BH")
  Hs_iter <- 1 - Ps_adj_tmp*(1-pi1s_iter)
    
  for (i in 1:maxiter){
    # fit GLM on transformed and imputed responses
    glm_mus_iter <- glm( Ys ~ X1+X2, family=Gamma(), weights=Hs_iter, data=Xs)
    # update alphas and Hs
    alphas_iter <- predict(glm_mus_iter, type="link")
    alphas_iter <- pmin( alphas_iter, alpha_max)
    
    # E-step for Hs
    Hs_iter[uncensored_locs] <- 1-(1 - pi1s_iter[uncensored_locs])/((1 - pi1s_iter[uncensored_locs])  + 
                                                                      pi1s_iter[uncensored_locs]*dbeta(censored_Ps[uncensored_locs], alphas_iter[uncensored_locs], 1))
    Hs_iter[censored_locs] <-1- (1 - pi1s_iter[censored_locs])*tau_censor/((1 - pi1s_iter[censored_locs])*tau_censor  + 
                                                                             pi1s_iter[censored_locs]*pbeta(tau_censor, alphas_iter[censored_locs], 1))
    
    # E-step for censored Ys
    Ys[censored_locs] <- (1/alphas_iter[censored_locs] - log(tau_censor))#/Hs_iter[censored_locs]
    
    # fit logistic GLM
    glm_pi1s_iter <- glm(Hs_iter ~ X1+X2, family=quasibinomial(), data=Xs)
    # get logistic GLM predictions
    pi1s_iter <-  predict(glm_pi1s_iter, type="response")
    pi1s_iter <- pmax(pi1_min, pmin(pi1_max, pi1s_iter))
    
    
  }
  list(glm_mus = glm_mus_iter, glm_pi1s = glm_pi1s_iter, glm_bl = glm_bl,
       alphas = alphas_iter, pi1s = pi1s_iter)
}

## get thresholds + weights
get_thresholds_betamix <- function(c, pi1s, alphas){
  pi0s <- 1 - pi1s
  ((pi0s/c - pi0s)/(alphas*pi1s))^(1/(alphas-1))
}

get_localfdrs_betamix <- function(ts, pi1s, alphas){
  pi0s <- 1 - pi1s
  pi0s/(pi0s + pi1s*dbeta(ts, alphas, 1))
}

get_tailfdr_betamix <- function(ts, pi1s, alphas, numerator_bh=FALSE){
  m <- length(ts)
  pi0s <- 1 - pi1s
  if (numerator_bh){
    pi0s_num <- rep(1,m)
  } else {
    pi0s_num <- pi0s
  }
  num <- sum(pi0s_num * ts)
  denom <- sum(pi0s*ts + pi1s*pbeta(ts, alphas, 1))
  num/denom
}


weights_betamix <- function(alpha, pi1s, alphas, numerator_bh=FALSE){
  wrapped_fun <- function(c) {
    ts <- get_thresholds_betamix(c, pi1s, alphas)
    get_tailfdr_betamix(ts, pi1s, alphas, numerator_bh = numerator_bh) - alpha
  }

  interval_endpts <- c(1e-10, 0.8)
  c_root <- uniroot(wrapped_fun, interval_endpts)$root
  ts_root <- get_thresholds_betamix(c_root, pi1s, alphas)
  ws <- ts_root/sum(ts_root)*length(ts_root)
  ws
}