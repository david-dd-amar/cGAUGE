try({library(mixtools)})

# Modified code from znormix (pi0 was depricated from CRAN on 2018)
# Alternatives: spEMsymlocN01, locfdr, fdrtool
fix_inf_z <- function(z){
  z[is.infinite(z) & z < 0] = min(z[is.finite(z)])
  z[is.infinite(z) & z > 0] = max(z[is.finite(z)])  
  return(z)
}

modified_znormix <- function (p, z = NULL,start.pi0=0.85, eps = 1e-04, 
            niter = 1000, verbose = FALSE,min_mean_z1=1,theoretical_null=T,
            min_sd1 = 0.25,max_sd1 = -1){
  if(!is.null(p)){z = -qnorm(p)}
  z = fix_inf_z(z)

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
  if(theoretical_null){
    last.par[2] = 0
    last.par[3] = 1
  }
  iter = 1
  new.par = last.par
  repeat {
    f0 = last.par[1] * dnorm(z, last.par[2], last.par[3])
    f1 = (1 - last.par[1]) * dnorm(z, last.par[4], last.par[5])
    ppee = pmin(1 - 1e-10,pmax(1e-10, f0/(f0 + f1)))
    new.par[1] = mean(ppee)
    sum.ppee = sum(ppee)
    new.par[4] = crossprod(z, 1 - ppee)/(G - sum.ppee)
    new.par[5] = sqrt(crossprod((z - new.par[4])^2, 1 - ppee)/(G - sum.ppee))
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
    
    if(min_mean_z1>0){
      new.par[4] = max(min_mean_z1,new.par[4])
    }
    if(min_sd1 > 0){
      new.par[5] = max(min_sd1,new.par[5])
    }
    if(max_sd1 > 0){
      new.par[5] = min(max_sd1,new.par[5])
    }
    
    if(theoretical_null){
      new.par[2] = 0
      new.par[3] = 1
    }
    else{
      new.par[3] = max(new.par[3],1)
    }

    if (isTRUE(verbose)) 
      cat("iter", iter, "\tparameters=", new.par, "\tmax.diff=", 
          max(abs(new.par - last.par)), fill = TRUE)
    if (iter >= niter || max(abs(new.par - last.par)) < eps) 
      break
    
    last.par = new.par
    iter = iter + 1
  }
  
  # ord = order(ppee)
  # fdr = numeric(G)
  # fdr[ord] = cumsum(ppee[ord])/(1:G)
  
  F0 = new.par[1] * pnorm(z, new.par[2], new.par[3],lower.tail = F)
  F1 = (1 - new.par[1]) * pnorm(z, new.par[4], new.par[5],lower.tail = F)
  fdr = pmin(1 - 1e-10,pmax(1e-10, F0/(F0 + F1)))
  names(new.par) = c("pi0", "mean.z0", "sd.z0", "mean.z1", "sd.z1")
  attr(new.par, "converged") = iter < niter
  attr(new.par, "iter") = iter
  attr(new.par, "call") = match.call()
  attr(new.par, "lfdr") = ppee
  attr(new.par, "fdr") = fdr
  class(new.par) = "znormix"
  new.par
}

znormix_wrapper<-function(p, z = NULL,
          start.pi0s = seq(0.85,0.88,0.01),reps=10,...){
  znormix_vals = c()
  for(start.pi0 in start.pi0s){
    for(i in 1:reps){
      currz = sample(z,replace = T)
      currp = sample(p,replace = T)
      m = modified_znormix(p=currp,z=currz,start.pi0 = start.pi0,...)
      znormix_vals = rbind(znormix_vals,m[1:5])
    }
  }
  new.par = colMeans(znormix_vals)
  if(!is.null(p)){z = -qnorm(p)}
  z[is.infinite(z) & z < 0] = min(z[is.finite(z)])
  z[is.infinite(z) & z > 0] = max(z[is.finite(z)])  
  f0 = new.par[1] * dnorm(z, new.par[2], new.par[3])
  f1 = (1 - new.par[1]) * dnorm(z, new.par[4], new.par[5])
  ppee = pmin(1 - 1e-10,pmax(1e-10, f0/(f0 + f1)))
  F0 = new.par[1] * pnorm(z, new.par[2], new.par[3],lower.tail = F)
  F1 = (1 - new.par[1]) * pnorm(z, new.par[4], new.par[5],lower.tail = F)
  fdr = pmin(1 - 1e-10,pmax(1e-10, F0/(F0 + F1)))
  names(new.par) = c("pi0", "mean.z0", "sd.z0", "mean.z1", "sd.z1")
  attr(new.par, "call") = match.call()
  attr(new.par, "lfdr") = ppee
  attr(new.par, "fdr") = fdr
  class(new.par) = "znormix"
  new.par
}

znormix_compute_fdrs<-function(z,new.params){
  f0 = new.params[1] * dnorm(z, new.params[2], new.params[3])
  f1 = (1 - new.params[1]) * dnorm(z, new.params[4], new.params[5])
  ppee = pmin(1 - 1e-10,pmax(1e-10, f0/(f0 + f1)))
  F0 = new.params[1] * pnorm(z, new.params[2], new.params[3],lower.tail = F)
  F1 = (1 - new.params[1]) * pnorm(z, new.params[4], new.params[5],lower.tail = F)
  fdr = pmin(1 - 1e-10,pmax(1e-10, F0/(F0 + F1)))
  names(new.params) = c("pi0", "mean.z0", "sd.z0", "mean.z1", "sd.z1")
  attr(new.params, "lfdr") = ppee
  attr(new.params, "fdr") = fdr
  class(new.params) = "znormix"
  new.params
}

normalmixEM_wrapper<-function(zz,reps=100,k,...){
  if(k==1){
    lambda=1
    mu = mean(zz,na.rm = T)
    s = sd(zz,na.rm = T)
    loglik = sum(dnorm(zz,mean=mu,sd=d,log = T))
    return(list(
      lambda=lambda,mu=mu,sd=s,loglik=loglik
    ))
  }
  bestModel = NULL
  for(j in 1:reps){
    try({
      currModel = normalmixEM(zz,k=k,...)
      if(is.null(bestModel) || currModel$loglik > bestModel$loglik){
        bestModel = currModel
      }
    })
  }
  return(bestModel)
}

#' Two groups estimation of the tendency of non-null observations from p1 to be null in p2.
#' 
#' We estimate the Bayes Fdr values of each point to infer two fdr vectors: fdr1 and fdr2.
#' Under the null hypothesis of no such events (Pr(h1=1,h2=0) = 0), fdr values from p1 can only
#' become lower in p2. We therefore use the Kolmogorov-Smirnov test on the pairwise fdr values: 
#' ks.test(fdr1,fdr2,alternative = "g").
simple_lfdr_test<-function(p1,p2,zthr = 10){
  z1 = c(-qnorm(p1))
  z1 = fix_inf_z(z1)
  z1[z1 > zthr] = zthr
  z1[z1 < -zthr] = -zthr
  z2 = c(-qnorm(p2))
  z2 = fix_inf_z(z2)
  z2[z2>zthr]  = zthr
  z2[z2<-zthr]  = -zthr
  z1_m = modified_znormix(p1,theoretical_null = T,min_mean_z1 = 1,min_sd1 = 0.5)
  z2_m = modified_znormix(p2,theoretical_null = T,min_mean_z1 = 1,min_sd1 = 0.5)
  # print(cbind(z1_m[1:5],z2_m[1:5]))
  fdr1 = attr(z1_m,"fdr")
  fdr2 = attr(z2_m,"fdr")
  
  # Null hypothesis: CDF of fdr1 lies below that of fdr2 - thus, values
  # in fdr1 are greater than those in fdr2. 
  # Alternative - we have non-null p1 cases that become null in p2: low fdrs in
  # p1 with high fdr in p2.
  pv = ks.test(fdr1,fdr2,alternative = "g",paired=T)$p.value
  return(pv)
}


# library(EnvStats)
# simple_z1_variance_test<-function(p1,p2,ltdr_val = 0.4){
#   z1 = c(-qnorm(p1))
#   z1 = fix_inf_z(z1)
#   z2 = c(-qnorm(p2))
#   z2 = fix_inf_z(z2)
#   z2_m = znormix_wrapper(p2,theoretical_null = T,
#             min_mean_z1 = 1,min_sd1 = 0.25,reps=2)
#   tdr2 = 1-attr(z2_m,"lfdr")
#   
#   inds = tdr2 < ltdr_val
#   print(table(inds))
#   p = varTest(z1[inds],alternative = "greater", sigma.squared =  1)$p.value
#   return(p)
# }


# # Check the tdr plots
# order_tdr_test<-function(p1,p2,reps=20,zthr = 10){
#   z1 = c(-qnorm(p1))
#   z1 = fix_inf_z(z1)
#   z1[z1 > zthr] = zthr
#   z2 = c(-qnorm(p2))
#   z2 = fix_inf_z(z2)
#   z2[z2>zthr]  = zthr
#   z1_m = znormix_wrapper(p1,theoretical_null = T,
#                          min_mean_z1 = 1,min_sd1 = 0.25,reps=2)
#   z2_m = znormix_wrapper(p2,theoretical_null = T,
#                          min_mean_z1 = 1,min_sd1 = 0.25,reps=2)
#   # compute our statistics
#   fdr1 = attr(z1_m,"fdr")
#   fdr2 = attr(z2_m,"fdr")
#   
#   tdr1 = 1-fdr1
#   tdr2 = 1-fdr2
#   ord1 = order(z1,decreasing = T)
#   ord2 = order(z2,decreasing = T)
#   F0 = cumsum(tdr1[ord1])/sum(tdr1)
#   F1 = cumsum(tdr1[ord2])/sum(tdr1)
#   plot(F0,type="l",ylab="tdr",xlab="")
#   lines(F1,type="l",col="blue")
#   diffs = F0-F1
#   lines(diffs)
#   abline(0,0,xpd=F)
#   STATISTIC = max(F0-F1)
# }

# # Too slow, under powered
# # # Code from mixtools (copied and simplified to avoid messages and unneeded options)
# # library(mixtools)
# # library(mixtools)
univar_mixtools_em<-function(p1,p2,zthr = 10,...){
  z1 = c(-qnorm(p1))
  z1 = fix_inf_z(z1)
  z1[z1 > zthr] = zthr
  z1[z1 < -zthr] = -zthr
  z2 = c(-qnorm(p2))
  z2 = fix_inf_z(z2)
  z2[z2>zthr]  = zthr
  z2[z2<-zthr]  = -zthr
  z1_m = znormix_wrapper(p1,theoretical_null = T,
                          min_mean_z1 = 1,min_sd1 = 0.25)[1:5]
  z2_m = znormix_wrapper(p2,theoretical_null = T,
                          min_mean_z1 = 1,min_sd1 = 0.25)[1:5]

  # define the differences
  zz = z1-z2
  try({
    possible_means = c(0,-z2_m[4],z1_m[4]-z2_m[4],z1_m[4])
    null_prior = c(0.8,rep(0.2/2,2))
    null_m = normalmixEM_wrapper(zz,lambda = null_prior,
                                 mean.constr = possible_means[1:3],k=3)
    alt_prior = c(0.8,rep(0.2/3,3))
    alt_m = normalmixEM_wrapper(zz,lambda = alt_prior,
                                mean.constr = possible_means[1:4],k=4)
    l_diff = alt_m$loglik - null_m$loglik
    l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
    return(l_diff_p)
  })

  return(1)
}


#' library(MASS)
#' rmvn_mix<-function(n,lambda,mu,sigma){
#'   samp = c()
#'   clustSamples = rmultinom(1,n,lambda)
#'   for(i in 1:nrow(clustSamples)){
#'     if(clustSamples[i,1]==0){next}
#'     currsamp = mvrnorm(clustSamples[i,1],mu[[i]],sigma[[i]])
#'     samp = rbind(samp,currsamp)
#'   }
#'   return(samp)
#' }
#' get_correl_from_diff_sd<-function(sd_diff,sd1,sd2){
#'   var_diff = sd_diff^2
#'   correl = (sd1^2 + sd2^2 - var_diff)/(2*sd1*sd2)
#'   return(correl)
#' }
#' 
#' #' Parametric bootstrap test for the existance of edge sep points.
#' #' 
#' #' This test transforms the input p-value vectors to z-scores (-qnorm(p)). 
#' #' For each z-score vector we fit a two groups model as a mixture of two gaussians: 
#' #' standard normal and a distribution of significant non-null z-scores.
#' #' We assume that this estimation is correct and take the two groups models as fixed.
#' #' Thus, each vector has three parameters: \pi (probability of null), \mu mean of discoveries
#' #' and \sigma - standard deviation of discoveries (assumed to be > 0.25). 
#' #' 
#' #' As we are interested in cases with "high" z1 and "insignificant" z2, we use the following test
#' #' statistic: \sum_i{tdr^1({z^1}_i)fdr^2({z^2}_i)}. Note that this statistic is not equivalent 
#' #' to summing over Pr(non-null|{z^1}_i)Pr(null|{z^2}_i) due to the dependence structure of z1 and z2.
#' #' However, as we explain below, we can simulate data points from the null distribution and compute
#' #' the significance of our statistic.
#' #' 
#' #' Assuming that z1 and z2 each follows a two groups model, 
#' #' their the joint distribution is a mixture of up to four Guassians:
#' #' 1. (null,null) - a mixture around (0,0) with known marginal variance and unknown covariance.
#' #' 2. (non-null,non-null) - a mixture around (\mu_1,\mu_2) with known marginal variance and unknown covariance.
#' #' 3. (null,non-null) - a mixture around (0,\mu_2)
#' #' 4. (non-null,null) - a mixture around (\mu_1,0)
#' #' Assuming that if cluster 3 exists in the data then its contribution to the statistic
#' #' is zero (because Pr(non-null) for low zscores is zero) then simulating without this component
#' #' can only increase the statistic under the null. We thus simulate from a mixture of 1 and 2 only.
#' #' We have two unknown parameters: the covariances within these mixtures. We estimate these using a four
#' #' mixture model for the difference zz = z1-z2. 
#' #' We fit three possible models: all with the component for cluster 4 and:
#' #' with cluster 1 only, with 1 and 2 and with 1-3. 
#' #' Model selection for the null is based on BIC.
#' #' Given the selected model, we simulate data from the joint distribution with
#' #' clusters 1 and 2 only. 
#' #' We then compute empirical p-values using the normal approximation (requires assuming
#' #' that z-scores are independent.)
#' param_bootstrap_test<-function(p1,p2,reps=20,zthr = 10,BICdiff = 10){
#'   z1 = c(-qnorm(p1))
#'   z1 = fix_inf_z(z1)
#'   z1[z1 > zthr] = zthr
#'   z2 = c(-qnorm(p2))
#'   z2 = fix_inf_z(z2)
#'   z2[z2>zthr]  = zthr
#'   z1_m = znormix_wrapper(p1,theoretical_null = T,
#'                          min_mean_z1 = 1,min_sd1 = 0.25,reps=2)
#'   z2_m = znormix_wrapper(p2,theoretical_null = T,
#'                          min_mean_z1 = 1,min_sd1 = 0.25,reps=2)
#'   # compute our statistics
#'   fdr1 = attr(z1_m,"lfdr")
#'   fdr2 = attr(z2_m,"lfdr")
#'   stat = sum((1-fdr1)*fdr2)
#'   
#'   # create the null distribution
#'   minp = min(z1_m[1],z2_m[1])
#'   lambda = c(minp,1-minp)
#'   mu = list(
#'     c(0,0),
#'     c(z1_m[4],z2_m[4])
#'   )
#'   
#'   # Look at the diff distribution
#'   zz = z1-z2
#'   possible_means = c(0,-z2_m[4],z1_m[4]-z2_m[4],z1_m[4])
#'   newk = 4
#'   new_prior = c(0.8,rep(0.2/newk,newk))
#'   
#'   # Fit the EM models for the mixtures
#'   # curr_EM2 = list(lambda=1,mu=0,sd=sd(zz,na.rm = T),
#'   #       loglik=sum(dnorm(zz,0,sd(zz,na.rm = T),log = T)))
#'   # curr_EM3 = normalmixEM_wrapper(zz,lambda = new_prior,
#'   #       mean.constr = possible_means[c(1,2)],k=2)
#'   # curr_EM4 = normalmixEM_wrapper(zz,lambda = new_prior,
#'   #       mean.constr = possible_means[c(1,3)],k=2)
#'   curr_EM5 = normalmixEM_wrapper(zz,lambda = new_prior,
#'                                  mean.constr = possible_means[1:3],k=3)
#'   # curr_EM2$BIC = log(length(zz))*1-2*curr_EM2$loglik
#'   # curr_EM3$BIC = log(length(zz))*3-2*curr_EM3$loglik
#'   # curr_EM4$BIC = log(length(zz))*3-2*curr_EM4$loglik
#'   # curr_EM5$BIC = log(length(zz))*5-2*curr_EM5$loglik
#'   curr_EM = curr_EM5
#'   selected_model = 5
#'   # if(curr_EM4$BIC < curr_EM$BIC+BICdiff){
#'   #   curr_EM = curr_EM4
#'   #   selected_model = 4
#'   # }
#'   # if(curr_EM3$BIC < curr_EM$BIC+BICdiff){
#'   #   curr_EM = curr_EM3
#'   #   selected_model = 3
#'   # }
#'   # if(curr_EM2$BIC < curr_EM$BIC+BICdiff){
#'   #   curr_EM = curr_EM5
#'   #   selected_model = 2
#'   # }
#'   
#'   if(selected_model == 5){
#'     zcor = get_correl_from_diff_sd(curr_EM$sigma[1],1,1)
#'     zcor1 = get_correl_from_diff_sd(curr_EM$sigma[3],z1_m[5],z2_m[5])
#'     zcor1 = min(zcor1,1)
#'     zcor = min(zcor,1)
#'     zcor1 = max(zcor1,0)
#'     zcor = max(zcor,0)
#'   }
#'   if(selected_model == 4){
#'     zcor = get_correl_from_diff_sd(curr_EM$sigma[1],1,1)
#'     zcor1 = get_correl_from_diff_sd(curr_EM$sigma[2],z1_m[5],z2_m[5])
#'     zcor1 = min(zcor1,1)
#'     zcor = min(zcor,1)
#'   }
#'   if(selected_model == 3 || selected_model == 2){
#'     zcor = get_correl_from_diff_sd(curr_EM$sigma[1],1,1)
#'     zcor = min(zcor,1)
#'     zcor1 = zcor
#'   }
#'   
#'   zcov1 = zcor1*z1_m[5]*z2_m[5]
#'   corrmat0 = matrix(c(1,zcor,zcor,1),2,2)
#'   corrmat1 = matrix(c(z1_m[5]^2,zcov1,zcov1,z2_m[5]^2),2,2)
#'   sigma = list(corrmat0,corrmat1)
#'   print(sigma)
#'   # # check the sigmas
#'   # for(i in 1:length(sigma)){
#'   #   print(sigma[[i]])
#'   #   print(det(sigma[[i]]))
#'   # }
#'   
#'   # parametric bootstrap
#'   b_scores = c()
#'   for(j in 1:reps){
#'     print(j)
#'     zzsamp = rmvn_mix(length(z1),lambda,mu,sigma)
#'     currz1 = zzsamp[,1]
#'     currz2 = zzsamp[,2]
#'     plot(currz1,currz2,xlim=c(-5,15),ylim=c(-5,15),
#'          main=paste(format(zcor1,digits = 3),
#'                     format(zcor,digits = 3),sep=","));abline(0,1)
#'     plot(z1,z2,xlim=c(-5,15),ylim=c(-5,15));abline(0,1)
#'     currz1_m = znormix_compute_fdrs(currz1,z1_m[1:5])
#'     currz2_m = znormix_compute_fdrs(currz2,z2_m[1:5])
#'     
#'     currfdr1 = attr(currz1_m,"lfdr")
#'     currfdr2 = attr(currz2_m,"lfdr")
#'     currstat = sum((1-currfdr1)*currfdr2)
#'     b_scores = c(b_scores,currstat)
#'     # if(currstat>stat){break}
#'   }
#'   
#'   # use normal approximation for p
#'   pval = min((1+sum(b_scores>=stat))/length(b_scores),
#'              pnorm(stat,mean(b_scores),sd(b_scores),lower.tail = F))
#'   return(pval)
#' }
#
# get_univar_diff_null_model<-function(zz,z1_m,z2_m,...){
#   # infer the correct null model
#   null_k1_likelihoods = dnorm(zz,mean(zz,na.rm = T),
#                               sd(zz,na.rm = T),log = T)
#   null_m =  list(mu = mean(zz,na.rm=T),
#                  sigma = sd(zz,na.rm=T),lambda=1,
#                  loglik = sum(null_k1_likelihoods)
#   )
#   
#   # try adding a the negative component
#   try({
#     curr_EM = normalmixEM_wrapper(zz,lambda = c(0.9,0.1),
#                                   mean.constr = c(0,-z2_m[4]),
#                                   k=2)
#     l_diff = curr_EM$loglik - null_m$loglik
#     l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
#     if(l_diff_p < 0.001){
#       null_m = curr_EM
#     }
#   })
#   
#   # try adding a the z1 z2 non-null shift component
#   try({
#     newk = length(null_m$lambda)+1
#     new_prior = c(0.8,rep(0.2/newk,newk))
#     new_means = c(0,-z2_m[4],z1_m[4]-z2_m[4])
#     if(newk == 2){
#       new_means = c(0,z1_m[4]-z2_m[4])
#     }
#     curr_EM = normalmixEM_wrapper(zz,lambda = new_prior,
#                                   mean.constr = new_means,
#                                   k=newk)
#     l_diff = curr_EM$loglik - null_m$loglik
#     l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
#     if(l_diff_p < 0.001){
#       null_m = curr_EM
#     }
#   })
#   return(null_m)
# }

# # Too slow:
# bivar_mixtools_em<-function(p1,p2,...){
#   z1 = c(-qnorm(p1))
#   z1 = fix_inf_z(z1)
#   z2 = c(-qnorm(p2))
#   z2 = fix_inf_z(z2)
#   z1_m = znormix_wrapper(p1,theoretical_null = T,
#                          min_mean_z1 = 1,min_sd1 = 0.25)[1:5]
#   z2_m = znormix_wrapper(p2,theoretical_null = T,
#                          min_mean_z1 = 1,min_sd1 = 0.25)[1:5]
#   
#   # define the differences
#   zz = cbind(z1,z2)
#   
#   # infer the correct null model
#   # no mixture
#   null_k1_likelihoods = log(dmvnorm(zz,colMeans(zz),cov(zz)))
#   null_m =  list(mu = mean(zz,na.rm=T),
#                  sigma = sd(zz,na.rm=T),lambda=1,
#                  loglik = sum(null_k1_likelihoods))
#   
#   # Two components - additional one for non-nulls
#   try({
#     curr_EM = mvnormalmixEM(zz,lambda = c(0.9,0.1),
#                  mu = list(c(0,0),c(z1_m[4],z2_m[4])),k=2,verb=T,epsilon = 1e-02)
#     l_diff = curr_EM$loglik - null_m$loglik
#     l_diff_p = pchisq(2*l_diff,4,lower.tail = F)
#     if(l_diff_p < 0.001){
#       null_m = curr_EM
#     }
#   })
#   
#   # try adding a the z1 z2 non-null shift component
#   try({
#     newk = length(null_m$lambda)+1
#     new_prior = c(0.8,rep(0.2/(newk-1),newk-1))
#     new_means = list(
#       c(0,0),
#       c(z1_m[4],z2_m[4]),
#       c(0,z2_m[4])
#     )
#     if(newk == 2){
#       new_means = list(
#         c(0,0),
#         c(0,z2_m[4])
#       )
#     }
#     curr_EM = mvnormalmixEM(zz,lambda = new_prior, mu = new_means,
#                                   k=newk,epsilon = 1e-02)
#     l_diff = curr_EM$loglik - null_m$loglik
#     l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
#     if(l_diff_p < 1e-05){
#       null_m = curr_EM
#     }
#   })
#   
#   # At this point we have the null model, now we try adding the new component
#   try({
#     newk = length(null_m$lambda)+1
#     new_prior = c(0.8,rep(0.2/(newk-1),newk-1))
#     new_means = list(
#       c(0,0),
#       c(z1_m[4],z2_m[4]),
#       c(0,z2_m[4]),
#       c(z1_m[4],0)
#     )
#     if(newk == 2){
#       new_means = list(
#         c(0,0),
#         c(0,z2_m[4]),
#         c(z1_m[4],0)
#       )
#     }
#     curr_EM = mvnormalmixEM(zz,lambda = new_prior, mu = new_means,
#                             k=newk,epsilon = 1e-02,verb = T)
#     l_diff = curr_EM$loglik - null_m$loglik
#     l_diff_p = pchisq(2*l_diff,4,lower.tail = F)
#     return(l_diff_p)
#   })
#   return(1)
# }
# 
# 
# univar_mixtools_em(
#   c(runif(10000),runif(1000)/1000000000),
#   c(runif(10000),runif(1000)/1000000)
# )
# 
# p1 = c(runif(10000),runif(1000)/10000)
# p2 = c(runif(10000),runif(950)/10000,runif(50))
# z1 = c(-qnorm(p1))
# z2 = c(-qnorm(p2))
# plot(z1,z2)
# univar_mixtools_em(p1,p2)

# # Use rstan for a mixture model
# # code based on
# # http://rpubs.com/kaz_yos/fmm2
# models = "
#   
#  data { 
#      // Hyperparameters for the hyper-priors
#      real alpha; 
#      real beta; 
#      real m; 
#      real s_squared; 
#      real<lower=0> dirichlet_alpha; 
#   
#      // Define variables in data 
#      // Number of observations (an integer) 
#      int<lower=0> n; 
#      // Outcome (a real vector of length n) 
#      real y[n]; 
#      // Number of latent clusters 
#      int<lower=1> H; 
#   
#      // Grid evaluation 
#      real grid_max; 
#      real grid_min; 
#      int<lower=1> grid_length; 
#  } 
#   
#  transformed data { 
#      real s; 
#      real grid_step; 
#   
#      s = sqrt(s_squared); 
#      grid_step = (grid_max - grid_min) / (grid_length - 1); 
#  } 
#   
#  parameters { 
#      // Define parameters to estimate 
#      // Population mean (a real number) 
#      vector[H] mu; 
#      // Population variance (a positive real number) 
#      real<lower=0> sigma_squared[H]; 
#      // Cluster probability 
#      simplex[H] Pi; 
#  } 
#   
#  transformed parameters { 
#      // Population standard deviation (a positive real number) 
#      real<lower=0> sigma[H]; 
#      // Standard deviation (derived from variance) 
#      sigma = sqrt(sigma_squared); 
#  } 
#   
#  model { 
#      // Temporary vector for loop use. Need to come first before priors. 
#      real contributions[H]; 
#   
#      // Prior part of Bayesian inference 
#      // All vectorized 
#      // Mean 
#      mu ~ normal(m, s); 
#      // sigma^2 has inverse gamma (alpha = 1, beta = 1) prior 
#      sigma_squared ~ inv_gamma(alpha, beta); 
#      // cluster probability vector 
#      Pi ~ dirichlet(rep_vector(dirichlet_alpha / H, H)); 
#   
#      // Likelihood part of Bayesian inference 
#      // Outcome model N(mu, sigma^2) (use SD rather than Var) 
#      for (i in 1:n) { 
#          // Loop over individuals 
#   
#            for (h in 1:H) { 
#                // Loop over clusters within each individual 
#                // Log likelihood contributions log(Pi[h] * N(y[i] | mu[h],sigma[h])) 
#                contributions[h] = log(Pi[h]) + normal_lpdf(y[i] | mu[h], sigma[h]); 
#            } 
#   
#            // log(sum(exp(contribution element))) 
#            target += log_sum_exp(contributions); 
#   
#      } 
#  } 
#   
#  generated quantities { 
#   
#      real log_f[grid_length]; 
#   
#      for (g in 1:grid_length) { 
#          // Definiting here avoids reporting of these intermediates. 
#          real contributions[H]; 
#          real grid_value; 
#   
#          grid_value = grid_min + grid_step * (g - 1); 
#          for (h in 1:H) { 
#              contributions[h] = log(Pi[h]) + normal_lpdf(grid_value | mu[h], sigma[h]); 
#          } 
#   
#          log_f[g] = log_sum_exp(contributions); 
#      } 
#   
#  } 
#  
# "
# 
# # range of values to examine estimated density.
# grid_max <- 20
# grid_min <- -20
# grid_length <- 100
# data = list(alpha = 10^(-3), beta = 10^(-3),
#             m = 0, s_squared = 10^(3),
#             n = length(zz),
#             y = zz,
#             dirichlet_alpha = 1,
#             H = 4,
#             grid_max = grid_max, grid_min = grid_min, grid_length = grid_length)
# fit = stan(model_code=models,data=data,chains=4)
# 
# 
# models4 = "
#   
# data { 
#   // Define variables in data 
#   // Number of observations (an integer) 
#   int<lower=0> n; 
#   // Outcome (a real vector of length n) 
#   real y[n];
# } 
# 
# parameters { 
#   // Define parameters to estimate 
#   // Population mean (a real number) 
#   real mu00;
#   real mu10;
#   real mu01;
#   // real mu11;
#   // Population variance (a positive real number) 
#   real<lower=0> sigma_squared00; 
#   real<lower=0> sigma_squared10; 
#   real<lower=0> sigma_squared01; 
#   // real<lower=0> sigma_squared11; 
#   // Cluster probability 
#   simplex[3] Pi;
# } 
# 
# transformed parameters { 
#   // Population standard deviation (a positive real number) 
#   real<lower=0> sigma00;
#   real<lower=0> sigma10;
#   real<lower=0> sigma01;
#   // real<lower=0> sigma11;
#   // Standard deviation (derived from variance) 
#   sigma00 = sqrt(sigma_squared00);
#   sigma10 = sqrt(sigma_squared10);
#   sigma01 = sqrt(sigma_squared01);
#   // sigma11 = sqrt(sigma_squared11);
# } 
# 
# model { 
#   // Temporary vector for loop use. Need to come first before priors. 
#   real contributions[3]; 
# 
#   // Prior part of Bayesian inference 
#   // All vectorized 
#   // Mean 
#   mu00 ~ normal(0, 0.2);
#   mu10 ~ normal(4, 1);
#   mu01 ~ normal(-4, 1);
#   //mu11 ~ normal(1, 0.5);
#   // sigma^2 has inverse gamma (alpha = 1, beta = 1) prior 
#   sigma_squared00 ~ inv_gamma(2, 1); 
#   sigma_squared10 ~ inv_gamma(3, 1); 
#   sigma_squared01 ~ inv_gamma(3, 1); 
#   //sigma_squared11 ~ inv_gamma(2, 1); 
#   // cluster probability vector 
#   Pi ~ dirichlet(rep_vector(1, 3)); 
# 
#   // Likelihood part of Bayesian inference 
#   // Outcome model N(mu, sigma^2) (use SD rather than Var) 
#   for (i in 1:n) {
#     // Log likelihood contributions log(Pi[h] * N(y[i] | mu[h],sigma[h])) 
#     contributions[1] = log(Pi[1]) + normal_lpdf(y[i] | mu00, sigma00);
#     contributions[2] = log(Pi[2]) + normal_lpdf(y[i] | mu10, sigma10);
#     contributions[3] = log(Pi[3]) + normal_lpdf(y[i] | mu01, sigma01);
#     // contributions[4] = log(Pi[4]) + normal_lpdf(y[i] | mu11, sigma11);
#     // log(sum(exp(contribution element))) 
#    target += log_sum_exp(contributions); 
#   } 
# } 
# 
# generated quantities { 
#   vector[3] contributions;
#   vector[n] log_lik;
# 
#   for (i in 1:n) {
#     // Log likelihood contributions log(Pi[h] * N(y[i] | mu[h],sigma[h])) 
#     contributions[1] = log(Pi[1]) + normal_lpdf(y[i] | mu00, sigma00);
#     contributions[2] = log(Pi[2]) + normal_lpdf(y[i] | mu10, sigma10);
#     contributions[3] = log(Pi[3]) + normal_lpdf(y[i] | mu01, sigma01);
#     // contributions[4] = log(Pi[4]) + normal_lpdf(y[i] | mu11, sigma11);
#     log_lik[i] = log_sum_exp(contributions);
#   }
# 
# } 
# 
# "
# 
# # range of values to examine estimated density.
# data = list(n = length(zz),y = zz)
# fit = stan(model_code=models4,data=data,chains=2,iter=100,warmup = 20)
# 
#
# zdiff = z1-z2
# loglik0 = sum(dnorm(zdiff,mean(zdiff),sd(zdiff),log = T))
# emEst1 = normalmixEM(zdiff,lambda = c(0.5,0.5),mean.constr = c(0,"a"),k=2,epsilon = 1e-8)
# m2 = spEMsymlocN01(z1,maxiter = 20)                      
#                          
