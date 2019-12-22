# Get all associations between tr1 and tr2 conditioned on
# sets S of size depth (>0)

print("trait pair analysis, usage:
      <RData file with objects><tr1><tr2><depth: integer > 0>
      <pval_thr><out_file>")

source("~/repos/cGAUGE/R/auxil_functions.R")
args = commandArgs(trailingOnly=TRUE)
print(paste("input args:",args))
# Define the parameters of the current run
load(args[1])
tr1=args[2]
tr2=args[3]
depth = as.numeric(args[4])
pval_thr = as.numeric(args[5])
out_f = args[6]
FUNC = run_ci_test_logistic_linear_discrete
TO_INCLUDE = c("sex","age") # Can overlap with the actual traits
TO_INCLUDE = setdiff(TO_INCLUDE,c(tr1,tr2))

conditioned_traits = colnames(alldata)[-c(1:2)]
conditioned_traits = conditioned_traits[!grepl("^PC\\d",conditioned_traits)]
conditioned_traits = setdiff(conditioned_traits,c(tr1,tr2))

# simple helper functions for helping with indexing
is_all_adj<-function(v){
  if(length(v)==1){return(T)}
  n = length(v)
  mx = v[n]
  return(all(v ==  ((mx-n+1):mx)))
}
# is_all_adj(1:3)
# is_all_adj(3:6)
# is_all_adj(c(4,7,8))
increase_index<-function(v,base){
  if(length(v)==1){return(v+1)}
  if(is_all_adj(v) && max(v)==base){return(v)}
  n = length(v)
  if(is_all_adj(v)){
    v[n] = v[n]+1
    v[1:(n-1)] = 1:(n-1)
    return(v)
  }
  i=1
  while(i < n && is_all_adj(v[1:i])){i=i+1}
  v[i-1] = v[i-1]+1
  if(i > 2){v[1:(i-2)] = 1:(i-2)}
  return(v)
}
# v = 1:5
# counter = 1
# vs = c()
# while(any(v!=8:12)){
#   v = increase_index(v,1)
#   counter=counter+1
#   print(v)
#   vs = rbind(vs,v)
# }
# counter == choose(12,5)

inds = 1:depth
counter = 1
sepset = c()
maxp=-1
Ntuple = choose(length(conditioned_traits),depth)

while(counter <= Ntuple){
  print(inds)
  xx = c(tr1,tr2,conditioned_traits[inds])
  curr_inds = !apply(is.na(alldata[,xx]),1,any)
  if(sum(curr_inds) > 10000){
    currp = NULL
    try({currp = FUNC(tr1,tr2,conditioned_traits[inds],data=alldata)})
    # If the analysis above failed then we need to process the dataset
    if(is.null(currp)){
      condSet = conditioned_traits[inds]
      yy = alldata[curr_inds,condSet]
      # zero variance cols
      if(!is.null(dim(yy))){
        zero_sd = apply(yy,2,sd)==0
      }
      else{
        zero_sd = sd(yy) == 0
      }
      
      # TODO: add this analysis in case we are using discrete tests:
      # disc_cols = apply(yy,2,function(x)length(unique(x))<20)
      # min_class = rep(nrow(yy),ncol(yy))
      # for(j in which(disc_cols)){min_class[j] = min(table(yy[,j]))}
      
      to_rem =  zero_sd
      if(all(to_rem)){
        print(paste("cannot condition on",condSet, "skipping"))
        next
      }
      else{
        condSet = condSet[!to_rem]
        currp = NULL
        try({currp = FUNC(tr1,tr2,condSet,data=alldata)})
      }
    }
    if(!is.null(currp)){
      maxp = max(currp,maxp)
      sepset = rbind(sepset,
                     c(paste(conditioned_traits[inds],collapse=","),currp))
    }
    
    if (length(TO_INCLUDE) > 0){
      currp = NULL
      condSet = unique(c(TO_INCLUDE,conditioned_traits[inds]))
      yy = alldata[curr_inds,condSet]
      if(!is.null(dim(yy))){
        zero_sd = apply(yy,2,sd)==0
      }
      else{
        zero_sd = sd(yy) == 0
      }
      to_rem = zero_sd
      if(all(to_rem)){
        print(paste("cannot condition on",condSet, "skipping"))
        next
      }
      else{
        condSet = condSet[!to_rem]
        try({
          currp = FUNC(tr1,tr2,condSet,data=alldata)
          sepset = rbind(sepset,
                      c(paste(condSet,collapse=","),currp))
        })
      }
      
      if(!is.null(currp)){
        maxp = max(currp,maxp)
        sepset = rbind(sepset,
                       c(paste(conditioned_traits[inds],collapse=","),currp))
      }
    }
  }
  # update the loop parameters
  inds = increase_index(inds,length(conditioned_traits))
  counter = counter+1
}
save(tr1,tr2,maxp,sepset,file=out_f)
