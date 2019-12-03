# Modified code from znormix (pi0 was depricated from CRAN on 2018)
# Alternatives: spEMsymlocN01, locfdr, fdrtool
modified_znormix <- function (p, z = NULL,start.pi0=0.85, eps = 1e-04, 
            niter = 1000, verbose = FALSE,min_mean_z1=1,theoretical_null=T,
            min_sd1 = 0.5){
  if(!is.null(p)){z = -qnorm(p)}
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

znormix_wrapper<-function(p, z = NULL,start.pi0s = seq(0.85,0.9,0.01),reps=10,...){
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
  attr(new.par, "converged") = iter < niter
  attr(new.par, "iter") = iter
  attr(new.par, "call") = match.call()
  attr(new.par, "lfdr") = ppee
  attr(new.par, "fdr") = fdr
  class(new.par) = "znormix"
  new.par
}

# library(EnvStats)
# simple_hg_test<-function(p1,p2){
#   z1 = c(-qnorm(p1))
#   z2 = c(-qnorm(p2))
#   z1_m = modified_znormix(p1,theoretical_null = T,min_mean_z1 = 2)
#   z2_m = modified_znormix(p2,theoretical_null = T,min_mean_z1 = 2)
#   fdr1 = attr(z1_m,"fdr")
#   fdr2 = attr(z2_m,"fdr")
#   selected1 = fdr1 < 0.1
#   z2_s1 = z2[selected1]
#   p = varTest(z2_s1,alternative = "greater", sigma.squared =  z2_m[5]^2)$p.value
#   return(p)
# }

#' Two groups estimation of the tendency of non-null observations from p1 to be null in p2.
#' 
#' We estimate the Bayes Fdr values of each point to infer two fdr vectors: fdr1 and fdr2.
#' Under the null hypothesis of no such events (Pr(h1=1,h2=0) = 0), fdr values from p1 can only
#' become lower in p2. We therefore use the Kolmogorov-Smirnov test on the pairwise fdr values: 
#' ks.test(fdr1,fdr2,alternative = "g").
simple_lfdr_test<-function(p1,p2){
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
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

simple_locfdr_test<-function(p1,p2){
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
  # z1_m = locfdr(z1,nulltype = 0,plot=0)
  # z2_m = locfdr(z2,nulltype = 0,plot=0)
  # pi0Est1 = z1_m$fp0[1,3]
  # pi0Est2 = z2_m$fp0[1,3]
  # pi0Est2SD = z2_m$fp0[2,3]
  # pval = pnorm(pi0Est1,mean = pi0Est2,sd = pi0Est2SD)
  try({
    zz = z1-z2
    zz_m = locfdr(zz[zz>0],nulltype = 0,plot=0)
    zscore = (1-zz_m$fp0[1,3]) / zz_m$fp0[2,3]
    return(pnorm(zscore,lower.tail = F))
  })
  return(1)
}
# # tests
simple_lfdr_test(
  c(runif(10000),runif(1000)/1000000000),
  c(runif(10000),runif(1000)/1000000)
)
simple_lfdr_test(c(runif(10000),runif(1000)/100000),
                 c(runif(10000),runif(1000)/10))

zdiff_test<-function(p1,p2){
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
  z1_m = modified_znormix(p1,theoretical_null = T)
  z2_m = modified_znormix(p2,theoretical_null = T)
  zdiff = pmax(z1-z2,0)
  zdiff = z1-z2
  z_diff_m = modified_znormix(z=zdiff,p=NULL,theoretical_null = F,min_mean_z1 = 0)
  pv = ks.test(fdr1,fdr2,alternative = "g",paired=T)$p.value
  return(pv)
}

# # Code from mixtools (copied and simplified to avoid messages and unneeded options)
# library(mixtools)
# library(mixtools)
univar_mixtools_em<-function(p1,p2,B=NULL){
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
  z1_m = znormix_wrapper(p1,theoretical_null = T,
                          min_mean_z1 = 1,min_sd1 = 0.25)[1:5]
  z2_m = znormix_wrapper(p2,theoretical_null = T,
                          min_mean_z1 = 1,min_sd1 = 0.25)[1:5]
  
  # define the differences
  zz = z1-z2
  
  # infer the correct null model
  null_k1_likelihoods = dnorm(zz,mean(zz,na.rm = T),sd(zz,na.rm = T),log = T)
  null_m =  list(mu = mean(zz,na.rm=T),
                 sigma = sd(zz,na.rm=T),lambda=1,
                 loglik = sum(null_k1_likelihoods)
  )
  
  # try adding a the negative component
  try({
    curr_EM = normalmixEM(zz,lambda = c(0.9,0.1),
                          mean.constr = c(0,-z2_m[4]),
                          k=2,epsilon = 1e-6)
    l_diff = curr_EM$loglik - null_m$loglik
    l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
    if(l_diff_p < 1e-05){
      null_m = curr_EM
    }
  })
  
  # try adding a the z1 z2 non-null shift component
  try({
    newk = length(null_m$lambda)+1
    new_prior = c(0.8,rep(0.2/newk,newk))
    new_means = c(0,-z2_m[4],z1_m[4]-z2_m[4])
    if(newk == 2){
      new_means = c(0,z1_m[4]-z2_m[4])
    }
    curr_EM = normalmixEM(zz,lambda = new_prior,
                          mean.constr = new_means,
                          k=newk,epsilon = 1e-6)
    l_diff = curr_EM$loglik - null_m$loglik
    l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
    if(l_diff_p < 1e-05){
      null_m = curr_EM
    }
  })
  
  # At this point we have the null model, now we try adding the new component
  try({
    newk = length(null_m$lambda)+1
    new_prior = c(0.8,rep(0.2/newk,newk))
    new_means = c(null_m$mu,z1_m[4])
    if(newk == 2){
      new_means = c(0,z2_m[4])
    }
    curr_EM = normalmixEM(zz,lambda = new_prior,
                          mean.constr = new_means,
                          k=newk,epsilon = 1e-6)
    l_diff = curr_EM$loglik - null_m$loglik
    l_diff_p = pchisq(2*l_diff,2,lower.tail = F)
    return(l_diff_p)
  })
  return(1)
}

univar_mixtools_em(
  c(runif(10000),runif(1000)/1000000000),
  c(runif(10000),runif(1000)/1000000)
)

p1 = c(runif(10000),runif(1000)/1000)
p2 = c(runif(10000),runif(950)/500,runif(50))
z1 = c(-qnorm(p1))
z2 = c(-qnorm(p2))
plot(z1,z2)
univar_mixtools_em(p1,p2)

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
