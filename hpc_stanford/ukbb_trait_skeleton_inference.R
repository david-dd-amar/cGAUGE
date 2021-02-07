# This script takes a matrix of traits and perform skeleton inference ignoring
# any genetic data (variants or PCs)
# Output has three objects:
# marginal association p-values
# maximal p-values for each pair (overl all tests)
# a matrix with all potential separating sets (p>1e-10) for each pair

source("~/repos/cGAUGE/R/auxil_functions.R")

# Define the input
input_data = "~/cgauge_resub/april2019_traits_data.RData"
out_path = "~/cgauge_resub/traits_skeleton/"
out_object_to_store_results = "~/cgauge_resub/Gs_skeleton.RData"

system(paste("mkdir",out_path))
FUNC = run_ci_test_logistic_linear_discrete
NJOB = 400
depth = 2
# this is an important threshold
# we ignore pairs with marginal association p > pthr
pthr = 0.01
to_include = c("sex","age")
load(input_data)

ls_objs = ls()
if(is.element("traits",ls_objs)){
  traits = union(traits,to_include)
}
if(!is.element("traits",ls_objs)){
  traits = names(code2phe_data)
}
n = length(traits)

pmax_network = matrix(NA,nrow=n,ncol=n)
rownames(pmax_network) = traits; colnames(pmax_network) = traits
sepsets = list()
for(nn in colnames(pmax_network)){
  sepsets[[nn]] = list()
}

# Initial filtering to avoid many unnecessary runs
for(i in 2:n){
  for(j in 1:(i-1)){
    xx = c(traits[i],traits[j])
    curr_inds = !apply(is.na(alldata[,xx]),1,any)
    if(sum(curr_inds) > 10000){
      pmax_network[i,j] = FUNC(traits[i],traits[j],NULL,data=alldata)
      pmax_network[j,i] = pmax_network[i,j]
    }
  }
  print(i)
}
print(paste("number of pairs remaining with p > 0.01:",sum(pmax_network > pthr,na.rm=T)/2))
save(pmax_network,file=paste(out_path,"initial_pmax_network.RData",sep=""))

load(paste(out_path,"initial_pmax_network.RData",sep=""))

############################################################################################
############################################################################################
############################################################################################
############################################################################################

exec_cmd_on_sherlock<-function(cmd,jobname,out_path,...){
  err_path = paste(out_path,jobname,".err",sep="")
  log_path = paste(out_path,jobname,".log",sep="")
  curr_cmd = paste(cmd)
  curr_sh_file = paste(out_path,jobname,".sh",sep="")
  sh_prefix = get_sh_prefix_one_node_specify_cpu_and_mem(err_path,log_path,...)
  print_sh_file(curr_sh_file, sh_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file,'&'))
}


# Now run an analysis for each pair: condition on a single other trait
for(i in 2:n){
  for(j in 1:(i-1)){
    if(is.na(pmax_network[i,j])){next}
    if(pmax_network[i,j]>pthr){next}
    tr1 = traits[i]
    tr2 = traits[j]
    curr_job_name = paste(tr1,"_",tr2,"_depth1",sep="")
    
    if(is.element(paste(curr_job_name,".RData",sep=""),set=list.files(out_path))){next}
    
    err_path = paste(out_path,curr_job_name,".err",sep="")
    log_path = paste(out_path,curr_job_name,".log",sep="")
    curr_cmd = paste(
      "Rscript ~/repos/cGAUGE/hpc_stanford/analyze_trait_pair_for_skeleton.R",
      input_data,tr1,tr2,1,pthr,paste(out_path,curr_job_name,".RData",sep="")
    )
    exec_cmd_on_sherlock(curr_cmd,curr_job_name,out_path,
                         Ncpu= 1,mem_size= 4000,time="02:00:00")
    
    job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
    num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
    while(num_active_jobs > NJOB){
      Sys.sleep(5)
      job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
      num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
    }
    
  }
  print(i)
}

# Wait for the single-level analyses to end
job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
while(num_active_jobs > 5){
  Sys.sleep(5)
  job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
  num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
}

# Now run an analysis for each pair but condition on a pair of other traits
for(i in 2:n){
  for(j in 1:(i-1)){
    if(is.na(pmax_network[i,j])){next}
    if(pmax_network[i,j]>pthr){next}
    tr1 = traits[i]
    tr2 = traits[j]
    curr_job_name = paste(tr1,"_",tr2,"_depth2",sep="")
    
    if(is.element(paste(curr_job_name,".RData",sep=""),set=list.files(out_path))){next}
    
    err_path = paste(out_path,curr_job_name,".err",sep="")
    log_path = paste(out_path,curr_job_name,".log",sep="")
    curr_cmd = paste(
      "Rscript ~/repos/cGAUGE/hpc_stanford/analyze_trait_pair_for_skeleton.R",
      input_data,tr1,tr2,depth,pthr,paste(out_path,curr_job_name,".RData",sep="")
    )
    exec_cmd_on_sherlock(curr_cmd,curr_job_name,out_path,
                         Ncpu= 1,mem_size= 4000,time="05:00:00")
    
    job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
    num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
    while(num_active_jobs > NJOB){
      Sys.sleep(5)
      job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
      num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
    }
  }
  print(i)
}

# Wait for the jobs to end
job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
while(num_active_jobs > 5){
  Sys.sleep(5)
  job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
  num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
}

# # cancel jobs
# currjobs =  system2("squeue",args = list("-u davidama"),stdout=TRUE)
# currjobs = currjobs[!grepl("bash",currjobs)]
# currjobs = sapply(currjobs,function(x)strsplit(x,split="\\s+",perl = T)[[1]][2])
# for(job in currjobs){
#   system(paste("scancel",job))
# }

############################################################################################
############################################################################################
############################################################################################
############################################################################################

# Read the results and save in a single object
load(paste(out_path,"initial_pmax_network.RData",sep=""))
marginal_pvalues = pmax_network
sepsets = list()
for(nn in colnames(pmax_network)){
  sepsets[[nn]] = list()
}

pthr2 = 1e-10 # remove all sepsets with p<pthr2 to save space

# conditioning on a single other trait
for(i in 2:n){
  print(i)
  for(j in 1:(i-1)){
    tr1 = traits[i]
    tr2 = traits[j]
    if(is.na(pmax_network[i,j]) || marginal_pvalues[i,j]>pthr){
      sepsets[[tr1]][[tr2]] = NULL
      sepsets[[tr2]][[tr1]] = NULL
      next
    }
    curr_job_name = curr_job_name = paste(tr1,"_",tr2,"_depth1",sep="")
    try({rm(sepset);rm(maxp)})
    maxp=NA;sepset=NULL
    try({(load(paste(out_path,curr_job_name,".RData",sep="")))})
    if(is.na(maxp)){
      print(paste("missing analysis for pair:",curr_job_name))
      next
    }
    if(!is.na(pmax_network[i,j])){maxp = max(maxp,pmax_network[i,j])}
    pmax_network[i,j] = maxp
    pmax_network[j,i] = maxp
    sepsets[[tr1]][[tr2]] = sepset
    sepsets[[tr2]][[tr1]] = sepset
  }
}

# conditioning on trait pairs
for(i in 2:n){
  print(i)
  for(j in 1:(i-1)){
    tr1 = traits[i]
    tr2 = traits[j]
    if(is.na(pmax_network[i,j]) || marginal_pvalues[i,j]>pthr){
      next
    }
    curr_job_name = curr_job_name = paste(tr1,"_",tr2,"_depth2",sep="")
    try({rm(sepset);rm(maxp)})
    maxp=NA;sepset=NULL
    try({(load(paste(out_path,curr_job_name,".RData",sep="")))})
    if(is.na(maxp)){
      print(paste("missing analysis for pair:",curr_job_name))
      next
    }
    if(!is.na(pmax_network[i,j])){maxp = max(maxp,pmax_network[i,j])}
    pmax_network[i,j] = maxp
    pmax_network[j,i] = maxp
    sepsets[[tr1]][[tr2]] = rbind(sepsets[[tr1]][[tr2]],sepset)
    sepsets[[tr2]][[tr1]] = sepsets[[tr1]][[tr2]]
  }
}

# clean the separating sets
# Note, the separating sets may contain duplications because:
#   We may have duplications for sets with sex, age 
#   In some subsets of the data we may have traits with zero variance
#   and these are removed.
#   By definition, in these cases, we use the first occurance of the set 
#   (because otherwise we subset by other traits first).
for(nn1 in names(sepsets)){
  for(nn2 in names(sepsets[[nn1]])){
    m = sepsets[[nn1]][[nn2]]
    if(is.null(dim(m))|| nrow(m)<2){next}
    
    g = m[,1]
    newm = tapply(m[,2],INDEX = g,FUN = function(x)x[1])
    newm = cbind(names(newm),newm)
    m = newm
    
    m_pvals = as.numeric(m[,2])
    to_rem = m_pvals <= pthr2
    m = m[!to_rem,]
    sepsets[[nn1]][[nn2]] = m
    if(is.null(dim(m))|| nrow(m)<2){next}
    # print(nrow(sepsets[[nn1]][[nn2]]) ==
    #         length(unique(sepsets[[nn1]][[nn2]][,1])))
  }
}

# Add the marginal associations to the separating sets
for(nn1 in names(sepsets)){
  for(nn2 in names(sepsets[[nn1]])){
    v = c("",marginal_pvalues[nn1,nn2])
    m = sepsets[[nn1]][[nn2]]
    if(is.null(dim(m))|| nrow(m)<2){next}
    m = rbind(v,m)
    sepsets[[nn1]][[nn2]] = m
  }
}

# QA and sanity checks
# sets are unique in each matrix
# p-values are numeric
for(nn1 in names(sepsets)){
  for(nn2 in names(sepsets[[nn1]])){
    m = sepsets[[nn1]][[nn2]]
    if(is.null(dim(m))|| nrow(m)<2){next}
    # m = unique(m)
    if(nrow(m)!=length(unique(m[,1]))){
      print(paste("Error in:",nn1,nn2))
    }
  }
}

# QA and sanity checks
# maximal observed p-value fits the sepset table
for(nn1 in names(sepsets)){
  for(nn2 in names(sepsets[[nn1]])){
    m = sepsets[[nn1]][[nn2]]
    currp = pmax_network[nn1,nn2]
    if(is.na(currp)){next}
    if(is.null(dim(m))|| nrow(m)<2){next}
    currps = as.numeric(m[,2])
    if(abs(currp - max(currps,na.rm = T))>1e-10){
      print(paste("Error in:",nn1,nn2))
    }
  }
}

save(pmax_network,marginal_pvalues,
     sepsets,file=out_object_to_store_results)

# # check file size
# system(paste("ls -lh",out_object_to_store_results))

# # compare to an older version of the pmax matrix
# all(colnames(P1)==colnames(P2))
# P1 = pmax_network
# load("~/cgauge_resub/april2019_Gs_skeleton.RData")
# P2 = pmax_network
# diffs = abs(P1-P2)
# quantile(c(diffs),na.rm = T) # there may be differences as P2 is more precise
# P1v = c(P1)
# P2v = c(P2)
# table(P1v<0.01,P2v<0.01) # diagonal matrix, perfect fit
# all((P1v<0.01) == (P2v<0.01),na.rm=T)
# all((P1v<0.001) == (P2v<0.001),na.rm=T)
# all((P1v<1e-07) == (P2v<1e-07),na.rm=T)
