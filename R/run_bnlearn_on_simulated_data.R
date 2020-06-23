required_libs = c("igraph","limma","bnlearn","optparse","Matrix")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
           })
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
get_distance_based_metrics<-function(arc_dists){
  arc_dists = as.numeric(arc_dists)
  fdr = sum(arc_dists<0)/length(arc_dists)
  fpr = sum(arc_dists>0)/length(arc_dists)
  ntp = sum(arc_dists>0)
  return(c(fdr,fpr,ntp))
}

get_bnlearn_results<-function(rdata,cut_breaks = 10,...){
  load(rdata)
  insts = grepl("IV",colnames(simulated_data))
  bls = c()
  for(iv in colnames(simulated_data)[insts]){
    bls = rbind(bls,cbind(colnames(simulated_data)[!insts],iv))
  }
  m1 = apply(simulated_data[,!insts],2,cut,breaks=cut_breaks)
  m2 = simulated_data[,insts]
  df = as.data.frame(cbind(m1,m2))
  
  dag = hc(df, score = "bic",blacklist = bls,...) 
  arcs = dag$arcs
  arcs_t = arcs[
    ! grepl("IV",arcs[,1]) & !grepl("IV",arcs[,2]),
    ]
  if(length(arcs_t)>0 && is.null(dim(arcs_t))){
    arcs_t = matrix(arcs_t,nrow=1)
  }
  if(length(arcs_t) != 0 && nrow(arcs_t)>0 ){
    arcs_t = add_distances(arcs_t,B_distances)
    res_with_ivs = get_distance_based_metrics(arcs_t[,3])
  }
  else{
    res_with_ivs = c(0,0,0)
  }
  
  df2 = as.data.frame(m1)
  dag2 = hc(df2, score = "bic",...) 
  arcs2 = dag2$arcs
  if(length(arcs2)>0 && is.null(dim(arcs2))){
    arcs2 = matrix(arcs2,nrow=1)
  }
  if(length(arcs2) != 0 && nrow(arcs2)>0){
    arcs2 = add_distances(arcs2,B_distances)
    res_without_ivs = get_distance_based_metrics(arcs2[,3])
  }
  else{
    res_without_ivs = c(0,0,0)
  }
  
  return(c(res_with_ivs,res_without_ivs))
}


library(parallel)

# Specify the directory with the simulation results, saved as 
# RData files using the full_causal_graph_simulations.R script
setwd("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_default/")
degs = c(1,1.5)
pleios = c(0,0.3)
p1s = c(0.001)
p2 = 0.1
param2res = list()
for(deg in degs){
  for(pleio in pleios){
    for(p1 in p1s){
      dir = paste0("deg",deg,"_pleio",pleio,"_p1",p1,"_p20.01/")
      rdata_files = list.files(dir,full.names = T)
      rdata_files = rdata_files[grepl("rdata",rdata_files,ignore.case = T)]
      
      res = mclapply(rdata_files,get_bnlearn_results,restart=100,
                     perturb=10,cut_breaks=5,mc.cores = 8)
      res = t(sapply(res,function(x)x))
      metric_names = c("FDR","TPR","N")
      colnames(res)= c(paste("bl",metric_names,sep="_"),paste("all",metric_names,sep="_"))
      print(apply(res,2,mean))
      param2res[[paste0("deg",deg,"_pleio",pleio)]] = res
      save(param2res,file="bnlearn_param2res_breaks5.RData")
    }
  }
}



# for(rdata in rdata_files){
#   #tmp = get_bnlearn_results(rdata)
# }


