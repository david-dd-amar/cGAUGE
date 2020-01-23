try({library(mixtools)})
try({library(mclust)})

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

znormix_wrapper<-function(p, z = NULL, Nsamp = 5000,
          start.pi0s = seq(0.85,0.88,0.01),reps=10,...){
  znormix_vals = c()
  Nsamp = min(Nsamp,max(length(p),length(z)))
  for(start.pi0 in start.pi0s){
    for(i in 1:reps){
      currz = sample(z,replace = T)[1:Nsamp]
      currp = sample(p,replace = T)[1:Nsamp]
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
    loglik = sum(dnorm(zz,mean=mu,sd=s,log = T))
    return(list(
      lambda=lambda,mu=mu,sd=s,loglik=loglik
    ))
  }
  bestModel = NULL
  logliks = c()
  for(j in 1:reps){
    try({
      currModel = normalmixEM(sample(zz),k=k,...)
      logliks[j] = currModel$loglik
      if(is.null(bestModel) || currModel$loglik > bestModel$loglik){
        bestModel = currModel
      }
    })
  }
  bestModel$logliks = logliks
  return(bestModel)
}

#' Two groups estimation of the tendency of non-null observations from p1 to be null in p2.
#' 
#' We estimate the Bayes Fdr values of each point to infer two fdr vectors: fdr1 and fdr2.
#' Under the null hypothesis of no such events (Pr(h1=1,h2=0) = 0), fdr values from p1 can only
#' become lower in p2. We therefore use the Kolmogorov-Smirnov test on the pairwise fdr values: 
#' ks.test(fdr1,fdr2,alternative = "g").
simple_lfdr_test<-function(p1,p2,zthr = 10){
  
  inds = !is.na(p1) & !is.na(p2)
  p1 = p1[inds];p2 = p2[inds]
  
  z1 = c(-qnorm(p1))
  z1 = fix_inf_z(z1)
  z1[z1 > zthr] = zthr
  z1[z1 < -zthr] = -zthr
  z2 = c(-qnorm(p2))
  z2 = fix_inf_z(z2)
  z2[z2 > zthr]  = zthr
  z2[z2 < -zthr]  = -zthr
  
  z1_m = znormix_wrapper(p=NULL,z=z1,theoretical_null = T,
                         min_mean_z1 = 1,min_sd1 = 0.25,reps = 50)
  z2_m = znormix_wrapper(p=NULL,z=z2,theoretical_null = T,
                         min_mean_z1 = 1,min_sd1 = 0.25,reps = 50)
  # print(cbind(z1_m[1:5],z2_m[1:5]))
  fdr1 = attr(z1_m,"lfdr")
  fdr2 = attr(z2_m,"lfdr")
  
  tdr1 = 1-fdr1
  tdr2 = 1-fdr2
  
  # Null hypothesis: under the null, z2 cannot have less discoveries
  # thus, its tdr CDF is equal or greater than that of z1
  pv = ks.test(tdr2,tdr1,alternative = "g",paired=T)$p.value
  return(pv)
}

# # Too slow, under powered
# # # Code from mixtools (copied and simplified to avoid messages and unneeded options)
# # library(mixtools)
univar_mixtools_em<-function(p1,p2,zthr = 10,return_models=F,...){
  
  inds = !is.na(p1) & !is.na(p2)
  p1 = p1[inds];p2 = p2[inds]
  z1 = c(-qnorm(p1))
  z1 = fix_inf_z(z1)
  z1[z1 > zthr] = zthr
  z1[z1 < -zthr] = -zthr
  z2 = c(-qnorm(p2))
  z2 = fix_inf_z(z2)
  z2[z2 > zthr]  = zthr
  z2[z2 < -zthr]  = -zthr
  
  z1_m = znormix_wrapper(p=NULL,z=z1,theoretical_null = F,
                         min_mean_z1 = 1,min_sd1 = 0.25,reps = 50)
  z2_m = znormix_wrapper(p=NULL,z=z2,theoretical_null = F,
                         min_mean_z1 = 1,min_sd1 = 0.25,reps = 50)

  # define the differences
  zz = z1-z2
  try({
    possible_means = c(z1_m[2]-z2_m[2],z1_m[2]-z2_m[4],z1_m[4]-z2_m[4],z1_m[4]-z2_m[2])
    alt_prior = c(0.8,rep(0.2/3,3))
    alt_m = normalmixEM_wrapper(zz,lambda = alt_prior,
                                mean.constr = possible_means[1:4],k=4,...)
    alt_m = emFixedMeans(zz,mu=possible_means)
    null_prior = c(0.8,rep(0.2/2,2))
    null_m = normalmixEM_wrapper(zz,lambda = null_prior,
                                 mean.constr = possible_means[1:3],k=3,...)
    l_diff = alt_m$loglik - null_m$loglik
    l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
    if(return_models){
      return(list(
        null_m=null_m,alt_m=alt_m,l_diff_p=l_diff_p,possible_means=possible_means)
      )
    }
    return(l_diff_p)
  })

  return(1)
}

# # Some tests
# setwd("~/Desktop/causal_inference_projects/ms3/edge_sep_em/")
# load("./Lipoprotein_A_LDL_direct_input.RData")
# load("./Lipoprotein_A_LDL_direct_edgesep_em_output.RData")
# prlist = read.delim("./plink.prune.in",stringsAsFactors = F)[,1]
# prlist = intersect(rownames(ps),prlist)
# p1 = ps[prlist,1];p2 = ps[prlist,2]
# grid_ms_test(p1,p2)

# allfiles = list.files("./")
# res_vs_sum = c()
# for(f in allfiles){
#   if(!grepl("_input.RData",f)){next}
#   currn = gsub("_input.RData","",f)
#   load(f)
#   load(paste(currn,"_edgesep_em_output.RData",sep=""))
#   p1 = ps[prlist,1];p2 = ps[prlist,2]
#   currs = sum(p1 < 1e-04 & p2>0.01,na.rm=T)
#   inds = !is.na(p1) & !is.na(p2)
#   p1 = p1[inds];p2 = p2[inds]
#   res_vs_sum = rbind(res_vs_sum,c(res,currs))
#   if(length(res)==0){next}
#   if(res<1e-10){
#     if(currs == 0){
#       print(f)
#     }
#   }
# }

# inds = !is.na(p1) & !is.na(p2)
# p1 = p1[inds];p2 = p2[inds]
# z1 = c(-qnorm(p1))
# z1 = fix_inf_z(z1)
# z2 = c(-qnorm(p2))
# z2 = fix_inf_z(z2)
# plot(z1,z2)
# univar_mixtools_em(p1,p2,reps=10)

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

mvnormix_e_step <- function(x, mu,covs, lambda) {
  comp_prods = c()
  for(j in 1:nrow(mu)){
    comp_prods = cbind(comp_prods,
                       mclust::dmvnorm(x, mu[j,], covs[[j]]) * lambda[j])
  }
  sum.of.comps <- rowSums(comp_prods)
  comp_posts = c()
  for(j in 1:nrow(mu)){
    comp_posts = cbind(comp_posts,
                       comp_prods[,j]/sum.of.comps)
  }
  
  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  
  list("loglik" = sum.of.comps.ln.sum,
       "posterior" = comp_posts)
}

mvnormix_lambda_m_step <- function(x, posterior.df) {
  return(colSums(posterior.df) / nrow(x))
}

lambda_em<-function(x,mu,cov,initial,maxiter=100,eps = 1e-03){
  e_step = mvnormix_e_step(x,mu,cov,initial)
  newl = mvnormix_lambda_m_step(x,e_step$posterior)
  mem = initial
  delta = sum((newl-mem)^2)
  iter = 1
  while(delta > eps && iter < maxiter){
    mem = newl 
    e_step = mvnormix_e_step(x,mu,cov,newl)
    newl = mvnormix_lambda_m_step(x,e_step$posterior)
    delta = sum((newl-mem)^2)
    iter = iter + 1
  }
  return(newl)
}

grid_bivar_normix_fixed_marginals<-function(z1,z2,z1_m,z2_m,cor_ranges = seq(0.5,0.9,0.1)){
  # define the mus
  zz = cbind(z1,z2)
  mu = rbind(
    c(z1_m[2],z2_m[2]),
    c(z1_m[2],z2_m[4]),
    c(z1_m[4],z2_m[4]),
    c(z1_m[4],z2_m[2])
  )
  covs = list(
    rbind(c(z1_m[3],NA),c(NA,z2_m[3])),
    rbind(c(z1_m[3],NA),c(NA,z2_m[5])),
    rbind(c(z1_m[5],NA),c(NA,z2_m[5])),
    rbind(c(z1_m[5],NA),c(NA,z2_m[3]))
  )
  
  null_models_loglik = c()
  alt_models_loglik = c()
  for(r00 in cor_ranges){
    cov00 = r00 * z1_m[3]*z2_m[3]
    for(r01 in cor_ranges){
      cov01 = r01 *z1_m[3]*z2_m[5]
      for(r11 in cor_ranges){
        cov11 = r11 *z1_m[5]*z2_m[5]
        curr_null_name = paste(r00,r01,r11,sep=",")
        curr_covs = list(
          rbind(c(z1_m[3],cov00),c(cov00,z2_m[3])),
          rbind(c(z1_m[3],cov01),c(cov01,z2_m[5])),
          rbind(c(z1_m[5],cov11),c(cov11,z2_m[5]))
        )
        null_lambda = c(0.8,0.1,0.1)
        est_lambda = lambda_em(zz,mu[-4,],curr_covs,null_lambda)
        curr_null_estep = mvnormix_e_step(x,mu[-4,],curr_covs,est_lambda)
        null_models_loglik[curr_null_name] = curr_null_estep$loglik
        
        for(r10 in cor_ranges){
          cov10 = r10 *z1_m[5]*z2_m[3]
          curr_alt_name = paste(r00,r01,r11,r10,sep=",")
          curr_alt_covs = list(
            rbind(c(z1_m[3],cov00),c(cov00,z2_m[3])),
            rbind(c(z1_m[3],cov01),c(cov01,z2_m[5])),
            rbind(c(z1_m[5],cov11),c(cov11,z2_m[5])),
            rbind(c(z1_m[5],cov10),c(cov10,z2_m[3]))
          )
          alt_initial_lambda = c(0.8,rep(0.2/3,3))
          alt_est_lambda = lambda_em(zz,mu,curr_alt_covs,alt_initial_lambda)
          curr_alt_estep = mvnormix_e_step(x,mu,curr_alt_covs,alt_est_lambda)
          alt_models_loglik[curr_alt_name] = curr_alt_estep$loglik
        }
      }
    }
  }
  
  liks = list(
    null_models_loglik = null_models_loglik,
    alt_models_loglik = alt_models_loglik
  )
  
  return(liks)
}


grid_ms_test<-function(p1,p2,zthr=10,marginal_em_reps = 20){
  inds = !is.na(p1) & !is.na(p2)
  p1 = p1[inds];p2 = p2[inds]
  z1 = c(-qnorm(p1))
  z1 = fix_inf_z(z1)
  z1[z1 > zthr] = zthr
  z1[z1 < -zthr] = -zthr
  z2 = c(-qnorm(p2))
  z2 = fix_inf_z(z2)
  z2[z2 > zthr]  = zthr
  z2[z2 < -zthr]  = -zthr
  z1_m = znormix_wrapper(p=NULL,z=z1,theoretical_null = F,
                         min_mean_z1 = 1,min_sd1 = 0.25,reps = marginal_em_reps)
  z2_m = znormix_wrapper(p=NULL,z=z2,theoretical_null = F,
                         min_mean_z1 = 1,min_sd1 = 0.25,reps = marginal_em_reps)
  
  liks = grid_bivar_normix_fixed_marginals(z1,z2,z1_m,z2_m)
  l_diff = max(liks$alt_models_loglik,na.rm = T) - 
    max(liks$null_models_loglik,na.rm = T)
  l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
  # print(l_diff_p)
  return(l_diff_p)
}

tr1 = "T4"
tr2 = "T7"
p1 = GWAS_Ps[,tr2]
p2 = trait_pair_pvals[[tr2]][[tr1]][,1]

# library(MASS)
# rmvn_mix<-function(n,lambda,mu,sigma){
#   samp = c()
#   clustSamples = rmultinom(1,n,lambda)
#   for(i in 1:nrow(clustSamples)){
#     if(clustSamples[i,1]==0){next}
#     currsamp = mvrnorm(clustSamples[i,1],mu[[i]],sigma[[i]])
#     samp = rbind(samp,currsamp)
#   }
#   return(samp)
# }

# get_correl_from_diff_sd<-function(sd_diff,sd1,sd2){
#   var_diff = sd_diff^2
#   correl = (sd1^2 + sd2^2 - var_diff)/(2*sd1*sd2)
#   return(correl)
# }

# Parametric bootstrap test for the existance of edge sep points.
# 
# This test transforms the input p-value vectors to z-scores (-qnorm(p)).
# For each z-score vector we fit a two groups model as a mixture of two gaussians:
# standard normal and a distribution of significant non-null z-scores.
# We assume that this estimation is correct and take the two groups models as fixed.
# Thus, each vector has three parameters: \pi (probability of null), \mu mean of discoveries
# and \sigma - standard deviation of discoveries (assumed to be > 0.25).
# 
# As we are interested in cases with "high" z1 and "insignificant" z2, we use the following test
# statistic: \sum_i{tdr^1({z^1}_i)fdr^2({z^2}_i)}. Note that this statistic is not equivalent
# to summing over Pr(non-null|{z^1}_i)Pr(null|{z^2}_i) due to the dependence structure of z1 and z2.
# However, as we explain below, we can simulate data points from the null distribution and compute
# the significance of our statistic.
# 
# Assuming that z1 and z2 each follows a two groups model,
# their the joint distribution is a mixture of up to four Guassians:
# 1. (null,null) - a mixture around (0,0) with known marginal variance and unknown covariance.
# 2. (non-null,non-null) - a mixture around (\mu_1,\mu_2) with known marginal variance and unknown covariance.
# 3. (null,non-null) - a mixture around (0,\mu_2)
# 4. (non-null,null) - a mixture around (\mu_1,0)
# Assuming that if cluster 3 exists in the data then its contribution to the statistic
# is zero (because Pr(non-null) for low zscores is zero) then simulating without this component
# can only increase the statistic under the null. We thus simulate from a mixture of 1 and 2 only.
# We have two unknown parameters: the covariances within these mixtures. We estimate these using a four
# mixture model for the difference zz = z1-z2.
# We fit three possible models: all with the component for cluster 4 and:
# with cluster 1 only, with 1 and 2 and with 1-3.
# Model selection for the null is based on BIC.
# Given the selected model, we simulate data from the joint distribution with
# clusters 1 and 2 only.
# We then compute empirical p-values using the normal approximation (requires assuming
# that z-scores are independent.)
# param_bootstrap_test<-function(p1,p2,reps=500,zthr = 10){
#   inds = !is.na(p1) & !is.na(p2)
#   p1 = p1[inds];p2 = p2[inds]
# 
#   z1 = c(-qnorm(p1))
#   z1 = fix_inf_z(z1)
#   z1[z1 > zthr] = zthr
#   z1[z1 < -zthr] = -zthr
#   z2 = c(-qnorm(p2))
#   z2 = fix_inf_z(z2)
#   z2[z2 > zthr]  = zthr
#   z2[z2 < -zthr]  = -zthr
# 
#   z1_m = znormix_wrapper(p=NULL,z=z1,theoretical_null = T,
#                          min_mean_z1 = 1,min_sd1 = 0.25,reps = 5)
#   z2_m = znormix_wrapper(p=NULL,z=z2,theoretical_null = T,
#                          min_mean_z1 = 1,min_sd1 = 0.25,reps = 5)
#   # compute our statistics
#   fdr1 = attr(z1_m,"lfdr")
#   fdr2 = attr(z2_m,"lfdr")
# 
#   stat = sum((1-fdr1)*fdr2)
# 
#   try({
#     # create the null distribution
#     zz = z1-z2
#     # define the differences
#     possible_means = c(0,-z2_m[4],z1_m[4]-z2_m[4],z1_m[4])
#     null_prior = c(0.8,rep(0.2/3,3))
#     null_m = normalmixEM_wrapper(zz,lambda = null_prior,
#                                mean.constr = possible_means,k=4,reps=5)
#     zcor = get_correl_from_diff_sd(null_m$sigma[1],1,1)
#     zcor1 = get_correl_from_diff_sd(null_m$sigma[3],z1_m[5],z2_m[5])
#     zcor2 = get_correl_from_diff_sd(null_m$sigma[2],1,z2_m[5])
#     zcor1 = min(zcor1,1);zcor = min(zcor,1);zcor2 = min(zcor2,1)
#     zcor1 = max(zcor1,0);zcor = max(zcor,0);zcor2 = max(zcor2,0)
#     zcov1 = zcor1*z1_m[5]*z2_m[5]
#     zcov2 = zcor2*z2_m[5]
#     corrmat0 = matrix(c(1,zcor,zcor,1),2,2)
#     corrmat1 = matrix(c(z1_m[5]^2,zcov1,zcov1,z2_m[5]^2),2,2)
#     corrmat2 = matrix(c(1,zcov2,zcov2,z2_m[5]^2),2,2)
#     sigma = list(corrmat0,corrmat2,corrmat1)
#     lambda = null_m$lambda[1:3]
#     lambda = lambda/sum(lambda)
#     mu = list(c(0,0),c(0,z2_m[4]),c(z1_m[4],z2_m[4]))
# 
#     # pi0 = (z1_m[1] + z2_m[1])*0.5
#     # pi0 = (m1$fp0[3,3] + m2$fp0[3,3])*0.5
#     # lambda = c(pi0,1-pi0)
#     # corrmat0 = matrix(c(1,0,0,1),2,2)
#     # corrmat1 = matrix(c(z1_m[5]^2,0,0,z2_m[5]^2),2,2)
#     # sigma = list(corrmat0,corrmat1)
#     # mu = list(c(0,0),c(z1_m[4],z2_m[4]))
#     # parametric bootstrap
#     b_scores = c();sum_b_greater = 0
#     for(j in 1:reps){
#       if(j %% 20 == 0){print(j)}
#       zzsamp = rmvn_mix(length(z1),lambda,mu,sigma)
#       currz1 = zzsamp[,1]
#       currz2 = zzsamp[,2]
#       # plot(currz1,currz2,xlim=c(-5,15),ylim=c(-5,15),
#       #    main=paste(format(zcor1,digits = 3),
#       #               format(zcor,digits = 3),sep=","));abline(0,1)
#       # plot(z1,z2,xlim=c(-5,15),ylim=c(-5,15));abline(0,1)
#       currz1_m = znormix_compute_fdrs(currz1,z1_m[1:5])
#       currz2_m = znormix_compute_fdrs(currz2,z2_m[1:5])
#       currfdr1 = attr(currz1_m,"lfdr")
#       currfdr2 = attr(currz2_m,"lfdr")
# 
#       currstat = sum((1-currfdr1)*currfdr2)
#       b_scores = c(b_scores,currstat)
#       sum_b_greater = sum_b_greater + as.numeric(currstat > stat)
#       emp_p = (1+sum_b_greater)/(1+length(b_scores))
#       if(j > 20 && emp_p > 0.1){break}
#     }
#     return(emp_p)
#   })
#   return(1)
# }
# 


#' #' Expectation Step of the EM Algorithm
#' #'
#' #' Calculate the posterior probabilities (soft labels) that each component
#' #' has to each data point.
#' #'
#' #' @param sd.vector Vector containing the standard deviations of each component
#' #' @param sd.vector Vector containing the mean of each component
#' #' @param alpha.vector Vector containing the mixing weights  of each component
#' #' @return Named list containing the loglik and posterior.df
#' e_step <- function(x, mu.vector, sd.vector, alpha.vector) {
#'   comp_prods = c()
#'   for(j in 1:length(mu.vector)){
#'     comp_prods = cbind(comp_prods,
#'                        dnorm(x, mu.vector[j], sd.vector[j]) * alpha.vector[j])
#'   }
#'   sum.of.comps <- rowSums(comp_prods)
#'   comp_posts = c()
#'   for(j in 1:length(mu.vector)){
#'     comp_posts = cbind(comp_posts,
#'                        comp_prods[,j]/sum.of.comps)
#'   }
#'   
#'   sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
#'   sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
#'   
#'   list("loglik" = sum.of.comps.ln.sum,
#'        "posterior.df" = comp_posts)
#' }
#' 
#' #' Maximization Step of the EM Algorithm
#' #'
#' #' Update the Component Parameters
#' #'
#' #' @param x Input data.
#' #' @param posterior.df Posterior probability data.frame.
#' #' @return Named list containing the mean (mu), variance (var), and mixing
#' #'   weights (alpha) for each component.
#' m_step <- function(x, posterior.df) {
#'   comp.n = c()
#'   comp.mu = c()
#'   comp.var = c()
#'   comp.alpha = c()
#'   for(j in 1:ncol(posterior.df)){
#'     comp.n = c(comp.n,sum(posterior.df[,j]))
#'     comp.mu = c(comp.mu,1/comp.n[j]* sum(posterior.df[, j] * x))
#'     comp.var = c(comp.var,sum(posterior.df[, j] * (x - comp.mu[j])^2) * 1/comp.n[j])
#'     comp.alpha = c(comp.alpha,comp.n[j] / length(x))
#'   }
#'   # comp1.n <- sum(posterior.df[, 1])
#'   # comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
#'   # comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
#'   # comp1.alpha <- comp1.n / length(x)
#'   
#'   list("mu" = comp.mu,
#'        "var" = comp.var,
#'        "alpha" = comp.alpha)
#' }
#' 
#' emFixedMeans<-function(x,maxiter = 100,mu,eps=1e-06){
#'   k = length(mu)
#'   for (i in 1:maxiter) {
#'     if (i == 1) {
#'       # Initialization via hard clustering
#'       initial_diffs = c()
#'       for(j in 1:k){
#'         initial_diffs = cbind(initial_diffs,abs(x - mu[j]))
#'       }
#'       initial_clusters = apply(initial_diffs,1,function(y)which(y==min(y))[1])
#'       # make sure each cluster is there
#'       initial_sds = c();initial_lambdas=c()
#'       for(j in 1:k){
#'         currinds = initial_clusters==j
#'         if(sum(currinds)<5){
#'           currsamp = sample(1:length(initial_clusters))[1:5]
#'           initial_clusters[currsamp]=j
#'           currinds[currsamp] = T
#'         }
#'         initial_sds[j] = sd(x[currinds])
#'         initial_lambdas[j] = sum(currinds)/length(currinds)
#'       }
#'       initial_lambdas = initial_lambdas/sum(initial_lambdas)
#'       e.step <- e_step(x, mu, initial_sds,initial_lambdas)
#'       m.step <- m_step(x, e.step[["posterior.df"]])
#'       cur.loglik <- e.step[["loglik"]]
#'       loglik.vector <- e.step[["loglik"]]
#'     } else {
#'       # Repeat E and M steps till convergence
#'       e.step <- e_step(x, mu, sqrt(m.step[["var"]]), 
#'                        m.step[["alpha"]])
#'       m.step <- m_step(x, e.step[["posterior.df"]])
#'       loglik.vector <- c(loglik.vector, e.step[["loglik"]])
#'       
#'       loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
#'       if(loglik.diff < eps) {
#'         break
#'       } else {
#'         cur.loglik <- e.step[["loglik"]]
#'       }
#'     }
#'   }
#'   return(list(
#'     loglik.vector = loglik.vector,
#'     loglik = cur.loglik,
#'     mu=mu,sd = sqrt(m.step[["var"]]),lambda = m.step[["alpha"]],
#'     posterior = e.step[["posterior.df"]]
#'   ))
#' }

