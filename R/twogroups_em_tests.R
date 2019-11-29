# Modified code from znormix (pi0 was depricated from CRAN on 2018)
modified_znormix <- function (p, z = NULL,start.pi0=0.85, eps = 1e-03, 
            niter = 1000, verbose = FALSE,min_mean_z1=1,theoretical_null=F,
            sd1GreaterThanSd2 = F){
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
    if(sd1GreaterThanSd2){
      new.par[5] = max(new.par[3],new.par[5])
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

library(EnvStats)
simple_hg_test<-function(p1,p2){
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
  z1_m = modified_znormix(p1,theoretical_null = T,min_mean_z1 = 2)
  z2_m = modified_znormix(p2,theoretical_null = T,min_mean_z1 = 2)
  fdr1 = attr(z1_m,"fdr")
  fdr2 = attr(z2_m,"fdr")
  selected1 = fdr1 < 0.1
  z2_s1 = z2[selected1]
  p = varTest(z2_s1,alternative = "greater", sigma.squared =  z2_m[5]^2)$p.value
  return(p)
}

library(locfdr)
#' Two groups estimation of the tendency of non-null observations from p1 to be null in p2.
#' 
#' We estimate the Bayes Fdr values of each point to infer two fdr vectors: fdr1 and fdr2.
#' Under the null hypothesis of no such events (Pr(h1=1,h2=0) = 0), fdr values from p1 can only
#' become lower in p2. We therefore use the Kolmogorov-Smirnov test on the pairwise fdr values: 
#' ks.test(fdr1,fdr2,alternative = "g").
simple_lfdr_test<-function(p1,p2){
  # p1 = round(p1,digits = 5)
  # p1[p1==0] = 1e-06
  # p1[p1==1] = 1-1e-06
  # p2 = round(p2,digits = 5)
  # p2[p2==0] = 1e-06
  # p2[p2==1] = 1-1e-06
  z1 = c(-qnorm(p1))
  z2 = c(-qnorm(p2))
  z1_m = modified_znormix(p1,theoretical_null = T)
  z2_m = modified_znormix(p2,theoretical_null = T)
  fdr1 = attr(z1_m,"fdr")
  fdr2 = attr(z2_m,"fdr")
  # try({
  #   locfdr1 = locfdr(z1,nulltype = 0,plot=0)
  #   locfdr2 = locfdr(z2,nulltype = 0,plot=0)
  #   fdr1 = locfdr1$fdr
  #   fdr2 = locfdr2$fdr
  # })

  # Null hypothesis: CDF of fdr1 lies below that of fdr2 - thus, values
  # in fdr1 are greater than those in fdr2. 
  # Alternative - we have non-null p1 cases that become null in p2: low fdrs in
  # p1 with high fdr in p2.
  pv = ks.test(fdr1,fdr2,alternative = "g",paired=T)$p.value
  return(pv)
}
# # tests
# simple_lfdr_test(c(runif(1000),runif(100)/1000000000),
#                  c(runif(1000),runif(100)/1000000000))
# simple_lfdr_test(c(runif(1000),runif(100)/100),
#                  c(runif(1000),runif(100)/100))
# simple_hg_test(c(runif(1000),runif(100)/1000000),
#                  c(runif(1000),runif(100)/100000))
# simple_hg_test(c(runif(1000),runif(100)/100),
#                  c(runif(1000),runif(100)/100))


# # Code from mixtools (copied and simplified to avoid messages and unneeded options)
# library(mixtools)
# library(mixtools)
# univar_mixtools_em<-function(p1,p2,B=NULL){
#   z1 = c(-qnorm(p1))
#   z2 = c(-qnorm(p2))
#   z1_m = modified_znormix(p1,theoretical_null = T)[1:5]
#   z2_m = modified_znormix(p2,theoretical_null = T)[1:5]
#   pval = 1
#   try({
#     zz = z1-z2
#     
#     params0 = list(
#       pro = c(0.8,0.1,0.1),
#       mean = c(0,-z2_m[4],z1_m[4]-z2_m[4])
#     )
#     emEst0 = normalmixEM(zz,lambda = params0$pro,mean.constr = params0$mean,
#                          sd.constr = c("a","a","a"),k=3,epsilon = 1e-4)
#     
#     params1 = list(
#       pro = c(0.8,0.05,0.1,0.05),
#       mean = c(0,-z2_m[4],z1_m[4]-z2_m[4],z1_m[4])
#     )
#     emEst1 = normalmixEM(zz,lambda = params1$pro,mean.constr = params1$mean,
#                          sd.constr = rep(emEst0$sigma[1],4),k=4,epsilon = 1e-4)
#     
#     l1 = emEst1$loglik
#     l0 = emEst0$loglik
#     dfs1 = 3 + 4 + 1
#     dfs2 = 2 + 3 + 1
#     chi = 2*(l1 - l0)
#     dof = dfs1-dfs2
#     pval = pchisq(chi,dof,lower.tail = F)
#     if(is.null(B) || B < 2){return(pval)}
#     l0_b = c()
#     for(j in 1:B){
#       emEst0_samp = normalmixEM(sample(zz,replace=T),
#                                 lambda = params0$pro,mean.constr = params0$mean,
#                                 sd.constr = c("a","a","a"),k=3,epsilon = 1e-4)
#       l0_b[j] = emEst0_samp$loglik
#     }
#     pval = sum(l1 > l0_b) / length(l0_b)
#   })
#   return(pval)
# }

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
# emEst1 = normalmixEM(z1,lambda = c(0.8,0.2),mean.constr = c(0,"a"),
#                          sd.constr = c("b","b"),k=2,epsilon = 1e-8)
                         
                         
