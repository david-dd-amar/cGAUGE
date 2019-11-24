# Modified code from znormix (pi0 was depricated from CRAN on 2018)
modified_znormix <- function (p, start.pi0=0.85, eps = 1e-03, 
                              niter = 1000, verbose = FALSE,min_mean_z1=1,theoretical_null=F){
  z = as.matrix(qnorm(1 - p))
  z[is.infinite(z) & z < 0] = min(z[is.finite(z)])
  z[is.infinite(z) & z > 0] = max(z[is.finite(z)])
  G = length(z)
  stopifnot(G >= 4)
  if (start.pi0 <= 0) 
    start.pi0 = 0.001
  if (start.pi0 >= 1) 
    start.pi0 = 1 - 0.001
  zcut = quantile(z, start.pi0)
  z0.idx = which(z < zcut)
  last.par = c(start.pi0, mean(z[z0.idx]), sd(drop(z[z0.idx])), 
               mean(z[-z0.idx]), sd(drop(z[-z0.idx])))
  iter = 1
  new.par = last.par
  repeat {
    f0 = last.par[1] * dnorm(z, last.par[2], last.par[3])
    ppee = pmin(1 - 1e-06, pmax(1e-06, f0/(f0 + (1 - last.par[1]) * 
                                             dnorm(z, last.par[4], last.par[5]))))
    new.par[1] = mean(ppee)
    sum.ppee = sum(ppee)
    new.par[4] = crossprod(z, 1 - ppee)/(G - sum.ppee)
    new.par[5] = sqrt(crossprod((z - new.par[4])^2, 1 - ppee)/(G - 
                                                                 sum.ppee))
    
    new.par[2] = crossprod(z, ppee)/sum.ppee
    new.par[3] = sqrt(crossprod((z - new.par[2])^2, ppee)/sum.ppee)
    if (abs(new.par[2]) > abs(new.par[4])) {
      tmp = new.par[2]
      new.par[2] = new.par[4]
      new.par[4] = tmp
      tmp = new.par[3]
      new.par[3] = new.par[5]
      new.par[5] = tmp
    }
    new.par[2] = 0
    new.par[3] = max(new.par[3],1)
    if(theoretical_null){
      new.par[3] = 1
    }
    new.par[4] = max(min_mean_z1,new.par[4])
    if (isTRUE(verbose)) 
      cat("iter", iter, "\tparameters=", new.par, "\tmax.diff=", 
          max(abs(new.par - last.par)), fill = TRUE)
    if (iter >= niter || max(abs(new.par - last.par)) < eps) 
      break
    last.par = new.par
    iter = iter + 1
  }
  ord = order(ppee)
  fdr = numeric(G)
  fdr[ord] = cumsum(ppee[ord])/(1:G)
  names(new.par) = c("pi0", "mean.z0", "sd.z0", "mean.z1", 
                     "sd.z1")
  attr(new.par, "converged") = iter < niter
  attr(new.par, "iter") = iter
  attr(new.par, "call") = match.call()
  attr(new.par, "lfdr") = ppee
  attr(new.par, "fdr") = fdr
  class(new.par) = "znormix"
  new.par
}


# x1 and x2 are lists or vectors or matrices
# x1 are the constaints - non NA vals
fixConstraints<-function(x1,x2){
  for(i in 1:length(x1)){
    inds = !is.na(x1[[i]])
    x2[[i]][inds] = x1[[i]][inds]
  }
  return(x2)
}

# Base code taken from mixtools
library(mixtools)
mvnormalmixConstrainedEM = 
  function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, 
            epsilon = 1e-06, maxit = 10000, verb = FALSE) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  tmp <- mvnormalmix.init(x = x, lambda = lambda, mu = mu, k = k)
  if(is.null(mu)){
    initial_mu = tmp$mu
    mu_constraints = F
  }
  else{
    initial_mu = mu
    mu_constraints = T
  }
  if(is.null(sigma)){
    initial_sigma = tmp$sigma
    sigma_constraints = F
  }
  else{
    initial_sigma = sigma
    sigma_constraints = T
  }
  lambda <- tmp$lambda
  lambda = lambda/sum(lambda)
  k = tmp$k
  diff <- 1
  iter <- 0
  mu = tmp$mu
  sigma = tmp$sigma
  comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x, mu[[i]], sigma[[i]]))
  comp <- sapply(comp, cbind)
  compsum <- apply(comp, 1, sum)
  obsloglik <- sum(log(compsum))
  ll <- obsloglik
  restarts <- 0
  while (diff > epsilon & iter < maxit) {
    # print(all(sapply(sigma,isSymmetric)))
      z = matrix(nrow = n, ncol = k)
      for (i in 1:n) {
        for (j in 1:k) {
          z.denom = c()
          for (m in 1:k) {
            z.denom = c(z.denom, lambda[m]/lambda[j] * 
                          (det(sigma[[j]])/det(sigma[[m]]))^(0.5) * 
                          exp(-0.5 * ((x[i, ] - mu[[m]]) %*% solve(sigma[[m]]) %*% 
                                        t(t(x[i, ] - mu[[m]])) - (x[i, ] - mu[[j]]) %*% 
                                        solve(sigma[[j]]) %*% t(t(x[i, ] - mu[[j]])))))
          }
          z[i, j] = 1/sum(z.denom)
        }
      }
      # print(table(is.na(z)))
      # print(table(is.infinite(z)))
      z = z/apply(z, 1, sum)
      sing <- sum(is.nan(z))
      lambda.new <- apply(z, 2, mean)
      # print(lambda.new)
      if (sum(lambda.new < 1e-08) > 0 || is.na(sum(lambda.new))) {
        lambda.new[is.na(lambda.new) | lambda.new < 1e-08] = 1e-6
      }
      mu.new <- lapply(1:k, function(j) sapply(1:p, 
                        function(i) apply(z * x[, i], 2, sum))[j, ]/sum(z[, j]))
      sigma.new <- lapply(1:k, 
                            function(j) matrix(apply(sapply(1:n, 
                                    function(i) z[i, j] * (x[i, ] - mu.new[[j]]) %*% 
                                                 t(x[i, ] - mu.new[[j]])), 1, sum), p, p)/sum(z[,j]))
      if(mu_constraints){
          mu.new = fixConstraints(initial_mu,mu.new)
      } 
      if(sigma_constraints){
          sigma.new = fixConstraints(initial_sigma,sigma.new)
      } 
      # print(sigma.new)
      
      lambda <- lambda.new
      mu <- mu.new
      # for(j in length(sigma)){
      #   if(det(sigma.new[[j]])>1e-06){
      #     sigma[[j]] = sigma.new[[j]]
      #   }
      # }
      sigma <- sigma.new
      comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,mu[[i]], sigma[[i]]))
      comp <- sapply(comp, cbind)
      compsum <- apply(comp, 1, sum)
      newobsloglik <- sum(log(compsum))
      
      diff <- newobsloglik - obsloglik
      obsloglik <- newobsloglik
      ll <- c(ll, obsloglik)
      iter <- iter + 1
      if (verb) {
        cat("iteration=", iter, "diff=", diff, "log-likelihood", 
            obsloglik, "\n")
      }
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
  # cat("number of iterations=", iter, "\n")
  a = list(x = x, lambda = lambda, mu = mu, sigma = sigma, 
           loglik = obsloglik, posterior = z, all.loglik = ll, restarts = restarts, 
           ft = "mvnormalmixEM")
  class(a) = "mixEM"
  a
}
