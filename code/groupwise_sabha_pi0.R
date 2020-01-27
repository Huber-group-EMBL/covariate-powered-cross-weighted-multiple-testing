groupwise_sabha <- function(pvals, Xs, alpha, tau=0.5, eps=0.1, return_fit=FALSE,
                            beta = 10^3, max_iters=400){
  # min -sum_i (B[i]*log((1-tau) q[i]) + (1-B[i])*log(1-(1-tau) q[i]))
  # subject to (1) q \in Qset (characterized by M*q \in Mset)
  # and (2) sum_i B[i]/q[i] <= gamma and (3) eps<=q<=1
  # introduce auxiliary variables x, y under the constraint Mq = x, q = y
  # ADMM optimization:
  # minimize -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + <u, Mq-x> + <v, q-y> + alpha/2 ||Mq-x||^2 + beta/2 ||q-y||^2 + alpha/2 (q-qt)'(eta I - M'M)(q-qt)
  # where qt is the previous iteration's q value
  
  # ADMM_params are: beta, max_iters
  Xs <- as.factor(Xs)

  n = length(pvals)
  m_groups <- length(unique(Xs))
  
  #B = (Pvals > tau) 
  pval_list <- split(pvals, Xs)
  q_init <- sapply(pval_list, function(pv) {sum( pv > tau)/length(pv)/(1-tau)})
  
  Ls = sapply(pval_list, function(ps) {sum(ps > tau)})     #split (Pvals > tau)
  Rs = sapply(pval_list, function(ps) {sum(ps <= tau)})     #split  (Pvals <= tau)
  gamma = n*(1-tau) # bound on sum_i (Pvals[i]>tau) / q[i]*(1-tau)
  # check if we need to return unit weights
  
  if (sum(Ls) > gamma){
    q <- rep(1, m_groups)
  } else {
    # initialize q
    q <- q_init
    # check if maybe we are done
    if (!(all(q >= eps) & all(q <= 1))){
       # continue with ADMM
      y = q 
      v = rep(0, m_groups)
  
      stop = FALSE
      iter = 0
      while(!stop){
        #print(iter)
        iter = iter+1
        q_old = q;  y_old = y;  v_old = v
        q = q_update(Ls, Rs, tau, eps, y, v, beta)
        y = y_update(Ls, q, y, v, beta, gamma)
        #return(y)
        v = v_update(v, q, y, beta)
        if(iter>=max_iters){stop=TRUE}
      }
    }
  }
  m_each_group <- lapply(pval_list, length)
  
  q_list <- lapply(1:m_groups, function(i) rep(q[i], m_each_group[i]))
  q_all <- unsplit(q_list, Xs)
  
  q_init_list <- lapply(1:m_groups, function(i) rep(q_init[i], m_each_group[i]))
  q_init_all <- unsplit(q_init_list, Xs)
  
  
  pvals[pvals>tau] <- Inf
  
  khat=max(c(0,which(sort(q_all*pvals)<=alpha*(1:length(pvals))/length(pvals))))
  
  rejected_idx = q_all*pvals <= alpha*khat/length(pvals)
  
  
  if (return_fit){
    res_list <-  list(q=q, q_all=q_all, q_init=q_init, q_init_all=q_init_all, rejected_idx=rejected_idx)
    return(res_list)
  } else {
    return(rejected_idx)
  }
}



q_update = function(L, R, tau, eps, y, v, beta){
  # minimize -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + <u, Mq-x> + <v, q-y> + alpha/2 ||Mq-x||^2 + beta/2 ||q-y||^2 + alpha/2 (q-qt)'(eta I - M'M)(q-qt)
  # where qt is the previous iteration's q value
  # equivalently, -sum_i (B_i*log((1-tau) q_i)+(1-B_i)*log(1-(1-tau) q_i)) + (alpha eta + beta)/2 * ||q-w||_2^2
  # where w = - (M'(ut + alpha (M qt - xt)) + (vt - beta yt - alpha eta qt))/(alpha eta + beta)
  m_groups <- length(L)
  q <- rep(1, m_groups)
  for (i in 1:m_groups){
    #print(paste0("L:",L[i]))
    #print(paste0("R:",R[i]))
    #print(paste0("v:",v[i]))
    #print(paste0("y:",y[i]))
    
    uniroot_fun <- function(c) {-L[i]/c + R[i]*(1-tau)/(1-c*(1-tau)) + v[i] + beta*(c - y[i])}
    q[i] = uniroot(uniroot_fun, c(0.0001, 1.5))$root
  }
  q[q<eps] = eps
  q[q>1] = 1
  q
}


y_update = function(L,q, y, v, beta, gamma){
  # Prof_B (q + v/beta)
  # where B = {sum_i B[i]/y[i]<= gamma}
  y = q + v/beta
  y = inverse_sum_prox(L, y, gamma)
  y
}



v_update = function(v, q, y, beta){
  v = v + beta * (q-y)
  v
}


# inverse_sum_prox solves: min{1/2 ||x-y||^2 : x_i>0, sum_i 1/x_i <= bound}
# Used in y-update step of ADMM
inverse_sum_prox = function(Ls, y,bound){
  
  y = pmax(0,y) # the solution will have all positive x_i's now
  # and we can now ignore the constraint x_i>0
  
  if(sum(Ls/y)<= bound){
    x=y
  }else{ # use Lagrange multipliers
    
    # we should have - lambda * d/dx_j (sum_i 1/x_i) = d/dx_j (1/2 ||x-y||^2)
    # for all j, for some single lambda>0
    # in other words, lambda / x^2 = x-y (this holds elementwise)
    # rearranging, lambda = x^3 - x^2*y
    # let c = log(lambda) so that it's real-valued
    # we need to solve x^3 - x^2*y - exp(c) = 0 (elementwise)
    
    cuberoot = function(c){ # this solves the cubic equation x^3-x^2*y-exp(c)=0
      temp1 = ((y/3)^3 + exp(c)/2 + (exp(c)*(y/3)^3 + exp(c)^2/4)^0.5)
      temp2 = ((y/3)^3 + exp(c)/2 - (exp(c)*(y/3)^3 + exp(c)^2/4)^0.5)
      x = sign(temp1)*abs(temp1)^(1/3) + sign(temp2)*abs(temp2)^(1/3) + (y/3)
      x
    }
    
    # now we need to choose c, i.e. choose the lagrange multiplier lambda=exp(c)
    # the right value of c is the one that produces an x satisfying sum_i 1/x_i = bound
    uniroot_fun <- function(c){sum(Ls/cuberoot(c+log(Ls)))-bound}
    #return(uniroot_fun)
    c = uniroot(uniroot_fun,c(-100,100))$root
    x = cuberoot(c+log(Ls))
  }
  x
}




