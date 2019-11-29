###################################################################################

# Set the session
required_libs = c("igraph","bnlearn","MRPRESSO",
                  "optparse","limma","MendelianRandomization")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
  })
}
# Add the cGAUGE functions and auxiliary functions for MR
source("https://raw.githubusercontent.com/david-dd-amar/cGAUGE/master/R/cGAUGE.R")
source("https://raw.githubusercontent.com/david-dd-amar/cGAUGE/master/R/twogroups_em_tests.R")
print("Completed loading libraries and code")

###################################################################################
# Simulation parameters

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list( 
  make_option(c("--N"), action="store", default=2000,type="integer",
              help="Sample size for simulated data"),
  make_option(c("--p"), action="store", default=15,type="integer",
              help="Number of phenotypes"),
  make_option(c("--deg"), action="store", default=1,type="double",
              help="Expected in/out degree in the causal graph, greater values mean more cycles"),
  make_option(c("--minBeta"), action="store", default=0.1,type="double",
              help="Min absolue value for causal effects (beta coefficients)"),
  make_option(c("--maxBeta"), action="store", default=0.95,type="double",
              help="Max absolue value for causal effects (beta coefficients)"),
  make_option(c("--minIVs"), action="store", default=10,type="integer",
              help="Min number of instruments per trait"),
  make_option(c("--maxIVs"), action="store", default=20,type="integer",
              help="Max number of instruments per trait"),
  make_option(c( "--minPleio"), action="store", default=1,type="integer",
              help="When applying pleiotropy, this is the min number of IV-trait links to add (at random)"),
  make_option(c( "--maxPleio"), action="store", default=10,type="integer",
              help="When applying pleiotropy, this is the max number of IV-trait links to add (at random)"),
  make_option(c( "--probPleio"), action="store", default=0.1,type="double",
              help="Probability that a variant is pleiotropic"),
  make_option(c( "--minMAF"), action="store", default=0.05,type="double",
              help="Minimal MAF - these are sampled with U(minMAF,maxMAF)"),
  make_option(c( "--maxMAF"), action="store", default=0.4,type="double",
              help="Maximal MAF - these are sampled with U(minMAF,maxMAF)"),
  make_option(c( "--p1"), action="store", default=0.001,type="double",
              help="P-value threshold for dependence - i.e., if p<p1"),
  make_option(c( "--p2"), action="store", default=0.01,type="double",
              help="P-value threshold for independence - i.e., if p>p2"),
  make_option(c( "--out"), action="store", default="simulation_result.RData",type="character",
              help="Output RData file with the simulated data and the analysis results"),
  make_option(c( "--edgeSepRun"), action="store", default="0",type="character",
              help="0: run edgeSep with the other methods; 1: run edge sep only, ignore p2")
  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

# sample size
N = opt$N
p = opt$p
deg = opt$deg
minBeta = opt$minBeta
maxBeta = opt$maxBeta
minIVs = opt$minIVs
maxIVs = opt$maxIVs
pleio_levelMin = opt$minPleio
pleio_levelMax=opt$maxPleio
prob_pleio = opt$probPleio
minMAF = opt$minMAF
maxMAF = opt$maxMAF
p1 = opt$p1
p2 = opt$p2
outfile = opt$out
edgeSepRun = opt$edgeSepRun

print("Completed parsing input parameters, starting the simulation")

###################################################################################
# Useful functions for the analysis below

# Functions for simulating data and obtaining the cGAUGR/MR results 
# Get the effect size, se, and p-value for x~y|z 
run_lm<-function(x,y,z,df){
  if(is.null(z)){
    df = data.frame(x=df[,x],y=df[,y])
  }
  else{
    df = data.frame(x=df[,x],y=df[,y],df[,z])
  }
  model = lm(x~.,data=df)
  coefs = summary(model)$coefficients
  return(coefs[2,])
}
# get distances (ancestry graph)
igraph_directed_distances<-function(Bg){
  distances = igraph::distances(Bg) # undirected distances
  for(i in 1:nrow(distances)){
    for(j in 1:ncol(distances)){
      distances[i,j] = length(igraph::shortest_paths(Bg,from=j,to=i)$vpath[[1]])-1
    }
  }
  return(distances)
}
# a function to add the known distances to the results
add_distances<-function(m,dists,newcolname = "KnownDistance"){
  m[,1] = as.character(m[,1])
  m[,2] = as.character(m[,2])
  v = c()
  for(i in 1:nrow(m)){
    v = c(v,dists[m[i,2],m[i,1]]) # add dist from exposure to outcome
  }
  m = cbind(m,v)
  colnames(m)[ncol(m)] = newcolname
  return(m)
}

###################################################################################
# Method for running LCV
source("https://raw.githubusercontent.com/lukejoconnor/LCV/master/R/RunLCV.R")
source("https://raw.githubusercontent.com/lukejoconnor/LCV/master/R/MomentFunctions.R")
exec_lcv_on_pair<-function(x1,x2,N,ell=NULL){
  inds = !is.na(x1) & !is.na(x2)
  x1=x1[inds];x2=x2[inds]
  M = length(x1)
  ## LCV parameters
  # Whether to use intercept in cross-trait LDSC regression. 
  # Should be 0 if there is no LD (generally, recommend setting this to 1 in real-data analyses)
  crosstrait.intercept <- 0
  # Whether to use intercept in univariate LDSC regression. 
  # Should be 0 if there is no LD (generally, recommend setting this to 1 in real-data analyses)
  ldsc.intercept <- 0 
  # Significance threshold: throw out SNPs about this threshold times mean
  #   chisq for computing LDSC intercept (all SNPs are included in other
  #   parts of the analysis)
  sig.threshold<-30
  # Number of jackknife blocks
  no.blocks<-100
  # Value of crosstrait LDSC intercept. Only applicable if crosstrait.intercept==1.
  # Should be zero for non-overlapping cohorts.
  intercept.12<-0;
  if(is.null(ell)){ell<-rep(1,M)} # LD scores
  weights <- rep(1,M) # regression weights
  N.1=N 
  N.2=N
  LCV<- RunLCV(ell,x1,x2,no.blocks,crosstrait.intercept,
               ldsc.intercept,weights,sig.threshold,N.1,N.2,intercept.12)
  logpval.1tail<-pt(LCV$zscore,no.blocks-2,log.p=1)
  return(c(LCV,logpval.1tail))
}
# A copy of LCV from Feb 2019 - We use this code to avoid the internal
# call of RunLCV to "source(Estimate4.R)"
# 
# RunLCV runs LCV on summary statistics for two traits.
#	Input vectors are sorted lists of LD scores and signed summary statistics. They
#		must be sorted by genomic position, as LCV uses a block-jackknife procedure
#		to compute standard errors; if consecutive SNPs are not approximately
#		contiguous, standard errors will be underestimated.
#
#   INPUT VARIBLES: ell, Mx1 vector of LD scores; 
#	z.1, Mx1 vector of estimated marginal per-normalized-genotype effects on trait 1
#   	(or Z scores; invariant to scaling);  
#	z.2, Mx2 vector of effects on trait2;  
#	crosstrait.intercept, 0 if cohorts are disjoint or overlap is known, 1 if cohorts 
#   	are possibly nondisjoint and necessary correction is unknown
#   ldsc.intercept, 0 if intercept should be fixed at 1 and 1 otherwise;
#   weights, Mx1 vector of regression weights; sig.threshold, threshold
#   	above which to discard chisq statistics for the purpose of estimating
#   	the LDSC intercept if they are above sig_threshold*mean(chisq);
#   no.blocks, number of jackknife blocks.
#   n1, sample size for trait 1, only needed if ldsc.intercept=1; 
#   n2, sample size for trait 2;  
#	intercept.12, covariance between sampling errors for Z1 and Z2, only needed if
#   	crosstrait_intercept=0.
#
#   OUTPUT VARIABLES: lcv.output, a data frame with columns:
#   "zscore", Z score for partial genetic causality; 
#   	to obtain 2-tailed p-value from zscore, compute:
#   	x<-pt(lcv.output[["zscore"]],no.blocks-2); pval<-2*min(x,1-x).
#   	Significantly positive zscore implies trait 1 partially genetically causal for trait 2.
#   "gcp.pm", posterior mean gcp (positive: trait 1 -> trait 2); 
#   "gcp.pse", posterior standard error for gcp; 
#   "rho.est", estimated genetic correlation; 
#   "rho.err", standard error of rho estimate;
#   "pval.fullycausal.1", p-value for null that gcp=-1; 
#   "pval.fullycausal.2", p-value for null that gcp=1; 
#   "h2.zscore.1", z score for trait 1 being significantly heritable;
#   	we recommend reporting results for h2.zscore.k > 7 (a very stringent threshold).
#   "h2.zscore.2", z score for trait 2 being significantly heritable.

RunLCV <- function(ell,z.1,z.2,no.blocks=100,crosstrait.intercept=1,ldsc.intercept=1,weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0,nargout=5){
  grid<- (-100:100)/100;
  mm=length(z.1)
  
  # Jackknife estimates of moments
  size.blocks=floor(mm/no.blocks)
  jackknife=matrix(0,no.blocks,8)
  for(jk in 1:no.blocks){
    if(jk==1)
    {ind<-(size.blocks+1):mm}
    else if(jk==no.blocks)
    {ind <- 1:(size.blocks*(jk-1))}
    else
    {ind<-c(1:((jk-1)*size.blocks), (jk*size.blocks+1):mm)}
    jackknife[jk,] <- EstimateK4(ell[ind],z.1[ind],z.2[ind],crosstrait.intercept,ldsc.intercept,weights[ind],sig.threshold,n.1,n.2,intercept.12,8)
  }
  
  jackknife[,2:3]<-jackknife[,2:3]-3*jackknife[,1]
  
  estimate<-mean(jackknife)
  error<-sd(jackknife)*sqrt(no.blocks+1)
  
  # Likelihood of each gcp value
  gcp.likelihood=grid;gcp.likelihood[]=0
  for(kk in 1:length(grid)){
    xx<-grid[kk]
    fx<-abs(jackknife[,1])^(-xx)
    numer<-jackknife[,2]/fx-fx*jackknife[,3]
    denom=pmax(1/abs(jackknife[,1]),sqrt(jackknife[,2]^2/fx^2 + jackknife[,3]^2*fx^2 ))
    pct.diff<-numer/denom
    
    gcp.likelihood[kk]<-dt(mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    
    if(xx==-1){
      pval.fullycausal.2<-pt(mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==1){
      pval.fullycausal.1<-pt(-mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==0){
      zscore<- mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1)
    }
  }
  
  rho.est<-mean(jackknife[,1])
  rho.err=sd(jackknife[,1])*sqrt(no.blocks+1)
  
  gcp.pm<-WeightedMean(grid,gcp.likelihood)
  gcp.pse<-sqrt(WeightedMean(grid^2,gcp.likelihood)-gcp.pm^2)
  
  h2.zscore.1<-mean(jackknife[,5])/sd(jackknife[,5])/sqrt(no.blocks+1)
  h2.zscore.2<-mean(jackknife[,6])/sd(jackknife[,6])/sqrt(no.blocks+1)
  
  lcv.output<-data.frame(zscore,gcp.pm,gcp.pse,rho.est,rho.err,
                         pval.fullycausal.1,pval.fullycausal.2,h2.zscore.1,h2.zscore.2)
  
  return(lcv.output)
}


###################################################################################
# create the basic pheno graph: edges are column to row
# j->i means B[i,j] != 0
B = matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    if(i==j){next}
    edgesign = (-1)^rbinom(1,1,0.5)
    edgeval = runif(1,minBeta,maxBeta)
    B[i,j] = edgesign*edgeval*rbinom(1,1,prob=deg/p)
  }
}
Bg = igraph::graph_from_adjacency_matrix(t(abs(B)>0)[1:p,1:p])
# plot(igraph::simplify(Bg))
# is.dag(Bg)
print("Completed creating the traits graph")
print("In-degrees:")
print(rowSums(abs(B)>0))
print("Out-degrees:")
print(colSums(abs(B)>0))
print("Is DAG?")
print(is_dag(Bg))

# Compute distances in B, among phenotypes
B_distances = igraph_directed_distances(Bg)
# sanity check - edges mean distance==1
# all((B_distances==1) == (abs(B)>0))

print("Adding instruments")
# sample variants per phenotype
for(i in 1:p){
  # add IVs
  curr_num_ivs = round(runif(1,min = minIVs,max=maxIVs))
  for(j in 1:curr_num_ivs){
    v = rep(0,nrow(B))
    edgesign = (-1)^rbinom(1,1,0.5)
    edgeval = runif(1,minBeta,maxBeta)
    v[i] = edgesign*edgeval
    # add pleiotropy
    is_pleio = rbinom(1,1,prob_pleio)
    if(is_pleio>0){
      pleio_level = round(runif(1,pleio_levelMin,pleio_levelMax))
      pleio_inds = sample(setdiff(1:p,i))[1:pleio_level]
      for(another_i in pleio_inds){
        edgesign = (-1)^rbinom(1,1,0.5)
        edgeval = runif(1,minBeta,maxBeta)
        v[another_i] = edgesign*edgeval
      }
    }
    # Expand B
    B = cbind(B,v)
    B = rbind(B,rep(0,ncol(B)))
  }
}
# rename dimensions
phenos = paste("T",1:p,sep="")
colnames(B)[1:p] = phenos
num_ivs = ncol(B) - p
ivs = paste("IV",1:num_ivs,sep="")
colnames(B)[(p+1):ncol(B)] = ivs
rownames(B) = colnames(B)
colnames(B_distances) = phenos
rownames(B_distances) = phenos
print("Done, out-degrees of instruments:")
print(colSums(abs(B[,ivs])>0))
print("Simulating the dataset using the graph")

# Simulate MAFs
mafs = runif(num_ivs,minMAF,maxMAF)
# parameters of the current SEM
I = diag(ncol(B))
A = solve(I-B)
theoretical_cov = A%*%t(A)
# Use the simulated B to generate data, use Fisher's fixed point method
simulated_data = c()
for(s in 1:N){
  v_t = rnorm(p)
  v_iv = rep(0,num_ivs)
  for(j in 1:num_ivs){v_iv[j] = rbinom(1,1,mafs[j])}
  e = c(v_t,v_iv)
  # # to make sure that the theoretical cov fits the calculation -
  # # this is a linear normal noise
  # e = rnorm(p+num_ivs) 
  v = (A%*%e)[,1]
  # This is not required for linear models:
  # for(kk in 1:100){
  #   v_prev = v
  #   v_new = (B%*%v_prev)[,1] + e
  #   v = v_new
  #   if(max(abs(v_prev-v_new))<1e-10){
  #     if(kk>1){print(kk)}
  #     break
  #   }
  # }
  simulated_data = rbind(simulated_data,v)
}
colnames(simulated_data) = colnames(B)
print("Done, computing associations and summary statistics")

###################################################################################
# Create the input for cGAUGE and MR
df = data.frame(simulated_data)
# Get all IV-phenotype associations
GWAS_Ps = matrix(1,num_ivs,p,dimnames = list(ivs,phenos))
GWAS_effects = matrix(0,num_ivs,p,dimnames = list(ivs,phenos))
GWAS_ses = matrix(0,num_ivs,p,dimnames = list(ivs,phenos))
GWAS_Zs = matrix(0,num_ivs,p,dimnames = list(ivs,phenos))
for(pheno in phenos){
  print(pheno)
  gwas_res = sapply(ivs,run_lm,x=pheno,z=NULL,df = df)
  GWAS_Ps[,pheno] = gwas_res[4,]
  GWAS_effects[,pheno] = gwas_res[1,]
  GWAS_ses[,pheno] = gwas_res[2,]
  GWAS_Zs[,pheno] = gwas_res[3,]
}

G_it = GWAS_Ps < p1
print("Done, number of associations with p<p1 per instrument:")
print(colSums(G_it))
print("Running standard MR analysis")

###########################
# cGAUGE starts here:
###########################

print("Starting the cGAUGE CI analysis")
# Skeleton learning
# G_t
print("Computing the trait skeleton matrix")
p_thr = p1
skeleton_pmax = matrix(-1,p,p,dimnames=list(phenos,phenos))
sepsets = list()
for(tr1 in phenos){
  sepsets[[tr1]] = list()
  for(tr2 in phenos){
    sepsets[[tr1]][[tr2]] = list()
  }
}
# Go over singletons
for(tr1 in phenos){
  for(tr2 in phenos){
    if(tr1==tr2){break}
    skeleton_pmax[tr1,tr2] = run_lm(tr1,tr2,NULL,simulated_data)[4]
    if(skeleton_pmax[tr1,tr2]>p_thr){
      skeleton_pmax[tr2,tr1] = skeleton_pmax[tr1,tr2]
      next
    }
    # go over singletons
    for(tr3 in phenos){
      if (tr3 %in% c(tr1,tr2)){next}
      currp = run_lm(tr1,tr2,tr3,simulated_data)[4]
      skeleton_pmax[tr1,tr2] = max(skeleton_pmax[tr1,tr2],currp)
      if(currp > p1){sepsets[[tr1]][[tr2]][[tr3]] = list(p=currp,sep=tr3)}
    }
    sepsets[[tr2]][[tr1]] = sepsets[[tr1]][[tr2]]
    skeleton_pmax[tr2,tr1] = skeleton_pmax[tr1,tr2]
  }
}
# Go over pairs
for(tr1 in phenos){
  for(tr2 in phenos){
    if(tr1==tr2){break}
    if(skeleton_pmax[tr1,tr2]>p_thr){next}
    for(tr3 in phenos){
      if (tr3 %in% c(tr1,tr2)){next}
      for(tr4 in phenos){
        if (tr4 %in% c(tr1,tr2,tr3)){next}
        currp = run_lm(tr1,tr2,c(tr3,tr4),simulated_data)[4]
        skeleton_pmax[tr1,tr2] = max(skeleton_pmax[tr1,tr2],currp)
        if(currp > p1){
          sepsets[[tr1]][[tr2]][[paste(tr1,tr2,sep=";")]] = 
            list(p=currp,sep=c(tr4,tr3))
        }
      }
    }
    sepsets[[tr2]][[tr1]] = sepsets[[tr1]][[tr2]]
    skeleton_pmax[tr2,tr1] = skeleton_pmax[tr1,tr2]
  }
}
G_t = skeleton_pmax < p1
print("Done, node degrees:")
print(colSums(G_t))

# Merge the sepsets
merged_sepsets = list()
for(tr1 in phenos){
  merged_sepsets[[tr1]] = list()
  for(tr2 in phenos){
    l = sepsets[[tr1]][[tr2]]
    if(length(l)==0){next}
    l = l[sapply(l,function(x)x$p)>p1]
    merged_sepsets[[tr1]][[tr2]] = unique(unlist(sapply(l,function(x)x$sep)))
  }
}

# G_vt
print("Computing all instrument vs trait pair CI tests:")
trait_pair_pvals = list()
for(pheno1 in phenos){
  trait_pair_pvals[[pheno1]] = list()
  for(pheno2 in phenos){
    if(pheno1==pheno2){next}
    gwas_res = t(sapply(ivs,run_lm,x=pheno1,z=pheno2,df = df))
    gwas_res = gwas_res[,4:1]
    trait_pair_pvals[[pheno1]][[pheno2]] = gwas_res
  }
}
G_vt = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,P1=p1,P2=p2,test_columns = 1)[[1]]
real_G_vt = abs(t(B[phenos,ivs])>0)

if(edgeSepRun=="1"){
  edge_sep_results_statTest = EdgeSepTest2(GWAS_Ps,G_t,trait_pair_pvals,text_col_name=1)
  edge_sep_results_statTest = add_distances(edge_sep_results_statTest,
                                            B_distances,newcolname = "KnownDistance")
  
  # boxplot(-log10(edge_sep_results_statTest$`pval:trait1->trait2`)~edge_sep_results_statTest$KnownDistance)
  edge_sep_results_statTest = edge_sep_results_statTest[
    p.adjust(edge_sep_results_statTest$`pval:trait1->trait2`)<0.01,]
  print("EdgeSep, Bonf correction (0.1), FDR and num discoveries:")
  print(paste(
    sum(edge_sep_results_statTest$KnownDistance==-1)/nrow(edge_sep_results_statTest),
    nrow(edge_sep_results_statTest)))
  
  save(
    opt, # input parameters
    B,Bg,simulated_data,B_distances, # simulated data
    edge_sep_results_statTest, # EdgeSepStatTest
    G_it,G_vt,G_t, iv_sets, # Skeletons
    file = outfile
  )
  
  q(save = "no",status = 0)
}


# Get new instrument sets after the cGAUGE filter
iv_sets = list()
for(tr1 in phenos){
  iv_sets[[tr1]] = list()
  for(tr2 in phenos){
    iv_sets[[tr1]][[tr2]] = rownames(GWAS_Ps)[GWAS_Ps[,tr1]<p1]
    currseps = merged_sepsets[[tr1]][[tr2]]
    # remove IVs into separating variables
    for(sep in currseps){
      curr_sep_ivs = rownames(G_vt)[G_vt[,sep]]
      iv_sets[[tr1]][[tr2]] = setdiff(iv_sets[[tr1]][[tr2]],curr_sep_ivs)
    }
    print(paste("before:",sum(G_it[,tr1]),"after:",length(iv_sets[[tr1]][[tr2]]),"sepNodes:",length(currseps)))
  }
}

# Run the MR
print("Done, rerunning MR")
# Pleio size is set to 1, to satisfy the conditions of Theorem 3.1
cgauge_mr_anal_res = list(
  "Egger" = run_pairwise_mr_analyses_with_iv_sets(GWAS_effects,GWAS_ses,iv_sets,
                                     func=mr_egger,robust=T),
  "IVW" = run_pairwise_mr_analyses_with_iv_sets(GWAS_effects,GWAS_ses,iv_sets,
                                   func=mr_ivw,robust=T)
)
# Add the known distances
cgauge_egger_res = cgauge_mr_anal_res$Egger
if(!is.null(dim(cgauge_mr_anal_res$Egger))){
  cgauge_egger_res = add_distances(cgauge_mr_anal_res$Egger,B_distances)
  cgauge_egger_res = add_distances(cgauge_egger_res,1-G_t,newcolname = "PleioProperty")
}

cgauge_ivw_res = cgauge_mr_anal_res$IVW
if(!is.null(dim(cgauge_mr_anal_res$IVW))){
  cgauge_ivw_res = add_distances(cgauge_ivw_res,B_distances)
  cgauge_ivw_res = add_distances(cgauge_ivw_res,1-G_t,newcolname = "PleioProperty")
}

print("Done with Egger and IVW, adding MRPRESSO")
cgauge_mrpresso_res = c()
try({
  # Add MRPRESSO
  for(tr1 in phenos){
    for(tr2 in phenos){
      if(tr1==tr2){next}
      currivs = iv_sets[[tr1]][[tr2]]
      if(length(currivs)<5){next}
      X = data.frame(E1b=GWAS_effects[currivs,tr1],O1b=GWAS_effects[currivs,tr2],
                     E1sd=GWAS_ses[currivs,tr1],O1sd=GWAS_ses[currivs,tr2])
      try({
        res = mr_presso(BetaOutcome = "O1b", BetaExposure = "E1b", 
                        SdOutcome = "O1sd", SdExposure = "E1sd",data=X,
                        OUTLIERtest=T,
                        DISTORTIONtest = T,
                        NbDistribution = 100,SignifThreshold = 0.1)
        if(is.na(res$`Main MR results`[2,"P-value"])){
          cgauge_mrpresso_res = rbind(cgauge_mrpresso_res,
                                      c(tr1,tr2,unlist(res$`Main MR results`[1,])))
        }
        else{
          cgauge_mrpresso_res = rbind(cgauge_mrpresso_res,
                                      c(tr1,tr2,unlist(res$`Main MR results`[2,])))
        }
      })
    }
  }
  if(!is.null(dim(cgauge_mrpresso_res))){
    cgauge_mrpresso_res = as.data.frame(cgauge_mrpresso_res)
    for(j in 3:ncol(cgauge_mrpresso_res)){
      cgauge_mrpresso_res[[j]] = as.numeric(as.character(cgauge_mrpresso_res[[j]]))
    }
  }
  # Add distances and the nonpleio property
  cgauge_mrpresso_res = add_distances(cgauge_mrpresso_res,B_distances)
  cgauge_mrpresso_res = add_distances(cgauge_mrpresso_res,1-G_t,"PleioProperty")
  
})

cgauge_mr_results = list()
try({cgauge_mr_results[["Egger"]] = cgauge_egger_res})
try({cgauge_mr_results[["IVW"]] = cgauge_ivw_res})
try({cgauge_mr_results[["MRPRESSO"]] = cgauge_mrpresso_res})

# EdgeSep
print("Done, running the skeletong edge separation analysis")
edge_sep_results = c()
try({
  # dummy_g_t = matrix(1,ncol(G_t),nrow(G_t),dimnames=dimnames(G_t))
  edge_sep_res = EdgeSep(GWAS_Ps,G_t,trait_pair_pvals,p1=p1,p2=p2,
                         text_col_name=1,pheno_names=NULL)
  for(nn in names(edge_sep_res)){
    arr = strsplit(nn,split=";")[[1]][c(1,3)]
    arr = c(arr,length(edge_sep_res[[nn]]$variants),
            edge_sep_res[[nn]]$num_tests,B_distances[arr[2],arr[1]])
    names(arr) = c("Exposure","Outcome","num_edgesep",
                   "num_exposure_variants","real_distance")
    edge_sep_results = rbind(edge_sep_results,arr)
  }
  if(!is.null(dim(edge_sep_results))){
    edge_sep_results = as.data.frame(edge_sep_results)
    for(j in 3:ncol(edge_sep_results)){
      edge_sep_results[[j]] = as.numeric(as.character(edge_sep_results[[j]]))
    }
  }
})

edge_sep_results_statTest = EdgeSepTest2(GWAS_Ps,G_t,trait_pair_pvals,text_col_name=1)
edge_sep_results_statTest = add_distances(edge_sep_results_statTest,
                                 B_distances,newcolname = "KnownDistance")

# boxplot(-log10(edge_sep_results_statTest$`pval:trait1->trait2`)~edge_sep_results_statTest$KnownDistance)
edge_sep_results_statTest = edge_sep_results_statTest[
  p.adjust(edge_sep_results_statTest$`pval:trait1->trait2`)<0.01,]
print("EdgeSep, Bonf correction (0.1), FDR and num discoveries:")
print(paste(
  sum(edge_sep_results_statTest$KnownDistance==-1)/nrow(edge_sep_results_statTest),
  nrow(edge_sep_results_statTest)))

#############################################################################
# Run MR methods
# pleio size is set to 100 - no filtering of variants (a standard MR analysis)
mr_anal_res = list(
  "Egger" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
                                     pleio_size=100,pruned_lists=NULL,func=mr_egger,robust=T),
  "IVW" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
                                   pleio_size=100,pruned_lists=NULL,func=mr_ivw,robust=T)
)
# Add MRPRESSO
print("Done with Egger and IVW, running MRPRESSO")
mrpresso_res = c()
try({
  for(tr1 in phenos){
    currivs = rownames(GWAS_Ps)[G_it[,tr1]]
    for(tr2 in phenos){
      if(tr1==tr2){next}
      X = data.frame(E1b=GWAS_effects[currivs,tr1],O1b=GWAS_effects[currivs,tr2],
                     E1sd=GWAS_ses[currivs,tr1],O1sd=GWAS_ses[currivs,tr2])
      try({
        res = mr_presso(BetaOutcome = "O1b", BetaExposure = "E1b", 
                        SdOutcome = "O1sd", SdExposure = "E1sd",data=X,
                        OUTLIERtest=T,
                        DISTORTIONtest = T,
                        NbDistribution = 100,SignifThreshold = 0.1)
        if(is.na(res$`Main MR results`[2,"P-value"])){
          mrpresso_res = rbind(mrpresso_res,
                               c(tr1,tr2,unlist(res$`Main MR results`[1,])))
        }
        else{
          mrpresso_res = rbind(mrpresso_res,
                               c(tr1,tr2,unlist(res$`Main MR results`[2,])))
        }
      })
    }
  }
  if(!is.null(dim(mrpresso_res))){
    mrpresso_res = as.data.frame(mrpresso_res)
    for(j in 3:ncol(mrpresso_res)){
      mrpresso_res[[j]] = as.numeric(as.character(mrpresso_res[[j]]))
    }  
  }
})

print("Done, adding LCV")
lcv_res = c()
try({
  # Add LCV
  for(tr1 in phenos){
    for(tr2 in phenos){
      if(tr1==tr2){next}
      curr_lcv = exec_lcv_on_pair(GWAS_Zs[,tr1],GWAS_Zs[,tr2],N=N)
      curr_p = pnorm(curr_lcv$zscore,lower.tail = F)
      lcv_res = rbind(lcv_res,
                      c(tr1,tr2,curr_lcv[[1]],curr_p,curr_lcv$gcp.pm,curr_lcv$gcp.pse,
                        curr_lcv$rho.est,curr_lcv$rho.err))
    }
  }
  colnames(lcv_res) = c("Trait1","Trait2","Z","P","gcp","gcpse","rho","rhose")
  if(!is.null(dim(lcv_res))){
    lcv_res = as.data.frame(lcv_res)
    for(j in 3:ncol(lcv_res)){
      lcv_res[[j]] = as.numeric(as.character(lcv_res[[j]]))
    }
  }
})
print("Done with LCV")

print("Adding the known trait distances to the MR data frames:")
# Save the MR results in a list
standard_mr_results = list()
try({standard_mr_results[["Egger"]] = add_distances(mr_anal_res$Egger,B_distances)})
try({standard_mr_results[["IVW"]] = add_distances(mr_anal_res$IVW,B_distances)})
try({standard_mr_results[["MRPRESSO"]] = add_distances(mrpresso_res,B_distances)})
try({standard_mr_results[["LCV"]] = add_distances(lcv_res,B_distances)})

#############################################################################
# Save the results in an RData file
save(
  opt, # input parameters
  B,Bg,simulated_data,B_distances, # simulated data
  cgauge_mr_results,standard_mr_results, # MR results
  edge_sep_results, # EdgeSep
  edge_sep_results_statTest, # EdgeSepStatTest
  G_it,G_vt,G_t, iv_sets, # Skeletons
  file = outfile
)


############################################################################
# explore the results (commented out, but can be used locally)
is_causal<-function(dists){
  return(dists>0 )
}

xx = standard_mr_results$MRPRESSO
# boxplot(-log10(xx$`P-value`)~xx$KnownDistance,main="MRPRESSO",las=2)
xx = xx[p.adjust(xx$`P-value`)< 0.1,]
print("MRPRESSO, Bonf correction (0.1), FDR and num discoveries:")
print(paste(sum(!is_causal(xx$KnownDistance))/nrow(xx),nrow(xx)))
xx = cgauge_mr_results$MRPRESSO
# boxplot(-log10(xx$`P-value`)~xx$KnownDistance,main="MRPRESSO",las=2)
xx = xx[p.adjust(xx$`P-value`)< 0.1,]
print("MRPRESSO+cGAUGE, Bonf correction (0.1), FDR and num discoveries:")
print(paste(sum(!is_causal(xx$KnownDistance))/nrow(xx),nrow(xx)))

par(mfrow=c(1,2))
xx = standard_mr_results$Egger
# cor.test(xx$p,xx$NumIVs,method = "spearman")$p.value
# boxplot(-log10(xx$p)~xx$KnownDistance,main="Egger",las=2)
xx = xx[p.adjust(xx$p)<0.1 & !is.na(xx$p),]
print("Egger, Bonf correction (0.1), FDR and num discoveries:")
print(paste(sum(!is_causal(xx$KnownDistance),na.rm = T)/nrow(xx),nrow(xx)))
xx = cgauge_mr_results$Egger
# cor.test(xx$p,xx$NumIVs,method = "spearman")$p.value
# boxplot(-log10(xx$p)~xx$KnownDistance,main="Egger+cGAUGE", las=2)
xx = xx[p.adjust(xx$p)<0.1 & !is.na(xx$p),]
print("Egger+cGAUGE, Bonf correction (0.1), FDR and num discoveries:")
print(paste(sum(!is_causal(xx$KnownDistance),na.rm = T)/nrow(xx),nrow(xx)))

par(mfrow=c(1,2))
xx = standard_mr_results$IVW
# cor.test(xx$p,xx$NumIVs,method = "spearman")$p.value
# boxplot(-log10(xx$p)~xx$KnownDistance,main="IVW",las=2)
xx = xx[p.adjust(xx$p)<0.1,]
print("IVW, Bonf correction (0.1), FDR and num discoveries:")
print(paste(sum(!is_causal(xx$KnownDistance),na.rm = T)/nrow(xx),nrow(xx)))

xx = cgauge_mr_results$IVW
# cor.test(xx$p,xx$NumIVs,method = "spearman")$p.value
# boxplot(-log10(xx$p)~xx$KnownDistance,main="IVW+cGAUGE", las=2)
xx = xx[p.adjust(xx$p)<0.1,]
print("IVW+cGAUGE, Bonf correction (0.1), FDR and num discoveries:")
print(paste(sum(!is_causal(xx$KnownDistance),na.rm = T)/nrow(xx),nrow(xx)))

# boxplot(-log10(edge_sep_results_statTest$`pval:trait1->trait2`)~edge_sep_results_statTest$KnownDistance)
edge_sep_results_statTest = edge_sep_results_statTest[
  p.adjust(edge_sep_results_statTest$`pval:trait1->trait2`)<0.01,]
print("EdgeSep, Bonf correction (0.1), FDR and num discoveries:")
print(paste(
  sum(edge_sep_results_statTest$KnownDistance==-1)/nrow(edge_sep_results_statTest),
  nrow(edge_sep_results_statTest)))

edge_sep_results = edge_sep_results[edge_sep_results$num_edgesep>2,]
print(paste(sum(edge_sep_results$real_distance==-1)/nrow(edge_sep_results),
            nrow(edge_sep_results)))
# 
# edge_sep_results_statTest = edge_sep_results_statTest[edge_sep_results_statTest$`pval:trait1->trait2` > 1,]
# print("EdgeSep, Bonf correction (0.1), FDR and num discoveries:")
# print(paste(sum(edge_sep_results_statTest$KnownDistance==-1)/nrow(edge_sep_results_statTest),nrow(edge_sep_results_statTest)))


# dev.off()
# plot(Bg)
# plot(simplify(igraph::graph_from_adjacency_matrix(G_t)),directed=F)
# plot(Bg)
# tr1 = "T13"
# tr2 = "T7"
# B[tr1,tr2];B[tr2,tr1]
# sepsets[[tr1]][[tr2]]
# absB = abs(B)>0
# table(absB[tr1,],absB[tr2,])
# chisq.test(table(absB[tr1,],absB[tr2,]))
# shared_ivs = colnames(B)[absB[tr1,] & absB[tr2,]]
# intersect(iv_sets[[tr1]][[tr2]],shared_ivs)
# intersect(iv_sets[[tr2]][[tr1]],shared_ivs)
# plot(GWAS_effects[,tr1],GWAS_effects[,tr2])
# cor.test(GWAS_effects[,tr1],GWAS_effects[,tr2])
# table(G_vt[,tr1],G_vt[,tr2])
# 
