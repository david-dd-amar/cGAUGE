required_libs = c("igraph","limma","bnlearn","Matrix")
# sherlock
# lib_loc = "~/R/packages3.5"
# lib_loc = c(lib_loc,.libPaths())
# scg
lib_loc = .libPaths()
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

# arc_dists: output of get_bnlearn_results: arcs ordered by their
# bootstrap proportion
get_distance_based_metrics<-function(arc_dists,topK){
  arc_dists = as.numeric(arc_dists[1:min(topK,length(arc_dists))])
  fdr = sum(arc_dists<0)/length(arc_dists)
  return(c(fdr,length(arc_dists)))
}

get_bnlearn_results<-function(rdata,cut_breaks = 5,boot.reps=100){
  load(rdata)
  insts = grepl("IV",colnames(simulated_data))
  bls = c()
  for(iv in colnames(simulated_data)[insts]){
    bls = rbind(bls,cbind(colnames(simulated_data)[!insts],iv))
  }
  # bls = rbind(bls,
  #   cbind(colnames(simulated_data)[insts],colnames(simulated_data)[insts]))
  m1 = apply(simulated_data[,!insts],2,cut,breaks=cut_breaks)
  m2 = simulated_data[,insts]
  df = as.data.frame(cbind(m1,m2))
  dags = list()
  for(j in 1:boot.reps){
    print(j)
    samp = sample(1:nrow(df),replace = T)
    dags[[j]] = hc(df[samp,], score = "bic",blacklist = bls)
  }
  t_arcs = c()
  for(j in 1:boot.reps){
    arcs = dags[[j]]$arcs
    arcs_t = arcs[! grepl("IV",arcs[,1]) & !grepl("IV",arcs[,2]),]
    arcs_string = apply(arcs_t,1,paste,collapse = "->")
    t_arcs = c(t_arcs,arcs_string)
  }
  tb = table(t_arcs) / length(dags)
  tb_n = sapply(names(tb),function(x)(strsplit(x,split="->")[[1]]))
  m = cbind(t(tb_n),tb)
  
  if(length(m) != 0 && nrow(m)>0 ){
    m = add_distances(m,B_distances)
    ord = order(as.numeric(m[,3]),decreasing = T)
    m = m[ord,]
    return(m)
  }
  else{
    return(NULL)
  }
}


library(parallel)

# Specify the directory with the simulation results, saved as 
# RData files using the full_causal_graph_simulations.R script
setwd("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_default/")
degs = c(1,1.5)
pleios = c(0,0.3)
ks = c(10,20)
param2res = list()
for(deg in degs){
  for(pleio in pleios){
      dir = paste0("deg",deg,"_pleio",pleio,"_p10.001_p20.01/")
      rdata_files = list.files(dir,full.names = T)
      rdata_files = rdata_files[grepl("rdata",rdata_files,ignore.case = T)]
      rdata_files = rdata_files[1:20]
      
      res = mclapply(rdata_files,get_bnlearn_results,cut_breaks=5,mc.cores = 20,boot.reps=200)
      
      k2scores = list()
      for(k in ks){
        k2scores[[as.character(k)]] = c()
        for(j in 1:length(res)){
          m = res[[j]]
          curr_scores = get_distance_based_metrics(as.numeric(m[,ncol(m)]),k)
          k2scores[[as.character(k)]] = rbind(
            k2scores[[as.character(k)]],curr_scores
          )
        }
      }
      param2res[[paste0("deg",deg,"_pleio",pleio)]] = list(res,k2scores)
      save(res,param2res,file="bnlearn_boot_breaks5.RData")
  }
}

# Create some plots locally
load("~/Desktop/causal_inference_projects/ms3/simulations_default/bnlearn_param2res_breaks5.RData")
load("~/Desktop/causal_inference_projects/ms3/simulations_default/bnlearn_param2res.RData")
sapply(param2res,function(x)apply(x,2,median))


t(sapply(k2scores,function(x)apply(x,2,median)))

