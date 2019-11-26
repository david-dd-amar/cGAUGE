# Modified code from znormix (pi0 was depricated from CRAN on 2018)
modified_znormix <- function (p, z = NULL,start.pi0=0.85, eps = 1e-03, 
                              niter = 1000, verbose = FALSE,min_mean_z1=1,theoretical_null=F){
  if(!is.null(p)){
    z = as.matrix(qnorm(1 - p))
    z[is.infinite(z) & z < 0] = min(z[is.finite(z)])
    z[is.infinite(z) & z > 0] = max(z[is.finite(z)])    
  }

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
# Code from mixtools (copied and simplified to avoid messages and unneeded options)
library(mixtools)
library(mixtools)
univar_mixtools_em<-function(p1,p2){
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
  z1_m = modified_znormix(p1,theoretical_null = T)[1:5]
  z2_m = modified_znormix(p2,theoretical_null = T)[1:5]
  pval = 1
  try({
    zz = z1-z2
    params1 = list(
      pro = c(0.8,0.05,0.05,0.1),
      mean = c(0,-z2_m[4],z1_m[4],z1_m[4]-z2_m[4]),
      variance = c(1,1,1,1)
    )
    emEst1 = normalmixEM(zz,lambda = params1$pro,mean.constr = params1$mean,
                         sd.constr = c("a","a","a","a"),k=4,epsilon = 1e-4)
    
    params0 = list(
      pro = c(0.8,0.1,0.1),
      mean = c(0,-z2_m[4],z1_m[4]-z2_m[4]),
      variance = c(1,1,1)
    )
    emEst0 = normalmixEM(zz,lambda = params0$pro,mean.constr = params0$mean,
                         sd.constr = c("a","a","a"),k=3,epsilon = 1e-4)
    
    l1 = emEst1$loglik
    l0 = emEst0$loglik
    dfs1 = 3 + 4 + 1
    dfs2 = 2 + 3 + 1
    chi = 2*(l1 - l0)
    dof = dfs1-dfs2
    pval = pchisq(chi,dof,lower.tail = F)
    
    for(j in 1:10){
      emEst0_samp = normalmixEM(zz,lambda = params0$pro,mean.constr = params0$mean,
                                sd.constr = c("a","a","a"),k=3,epsilon = 1e-4)
      emEst1_samp = normalmixEM(zz,lambda = params1$pro,mean.constr = params1$mean,
                                sd.constr = c("a","a","a","a"),k=4,epsilon = 1e-4)
      
      l1 = emEst1_samp$loglik
      l0 = emEst0_samp$loglik
      chi = 2*(l1 - l0)
      newp = pchisq(chi,dof,lower.tail = F)
      print(pval-newp)
      pval = max(pval,newp)
    }
  })
  pval
}


