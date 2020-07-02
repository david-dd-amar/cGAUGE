# This was done locally in a mac machine
# library(devtools)
# devtools::install_github("jean997/cause")
library(cause)
library(parallel)
# on the hpc
# devtools::install_github("jean997/cause",lib="~/")
library(cause,lib="~/")
library(parallel)

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

run_cause_on_tr1_ivs <- function(tr1,phenos,G_it,GWAS_effects,GWAS_ses,B_distances){
  cause_res = c()
  for(tr2 in phenos){
    if(tr1==tr2){next}
    if(is.null(G_it)){
      ivs = rownames(GWAS_effects)
    }
    else{
      ivs = rownames(G_it)[G_it[,tr1]]
    }
    
    X = cbind(
      GWAS_effects[ivs,tr1],GWAS_ses[ivs,tr1],
      GWAS_effects[ivs,tr2],GWAS_ses[ivs,tr2]
    )
    colnames(X) = c("beta_hat_1","seb1", "beta_hat_2","seb2")
    X = data.frame(X)
    X$snp = ivs
    X = new_cause_data(X)
    m_p = est_cause_params(X,variants = ivs)
    m = cause(X,m_p)
    p1 = pnorm(m$elpd[2,5])
    p2 = pnorm(m$elpd[3,5]) 
    m_s = summary(m)
    v = c(tr1,tr2,m_s$p,p1,p2,B_distances[tr2,tr1])
    print(v)
    cause_res = rbind(cause_res,v)
  }
  return(cause_res)
}

#setwd("~/Desktop/causal_inference_projects/ms3/simulations_default/simulations_default/")
setwd("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_default/")
degs = c(2,1,1.5)
pleios = c(0.4,0.3,0)
p1s = c(0.001)
for(deg in degs){
  for(pleio in pleios){
    for(p1 in p1s){
      dir = paste0("deg",deg,"_pleio",pleio,"_p1",p1,"_p20.01/")
      rdata_files = list.files(dir,full.names = T)
      rdata_files = rdata_files[grepl("rdata",rdata_files,ignore.case = T)]
      rdata_files = rdata_files
      file2res = list()
      rdata_out = paste0("cause_results_pleio",pleio,"_deg",deg,"_p1",p1,"_all_instruments.RData")
      try(load(rdata_out))
      for(rdata in rdata_files){
        
        if(rdata %in% names(file2res)){
          print(paste("rdata already analyzed:",rdata,"skipping"))
          next
        }
        
        load(rdata)
        
        # Create the input for cGAUGE and MR
        df = data.frame(simulated_data)
        ivs = colnames(df)[grepl("IV",colnames(df))]
        phenos = colnames(df)[!grepl("IV",colnames(df))]
        num_ivs = length(ivs)
        p = length(phenos)
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
        
        # use G_it to get instruments and run cause
        cause_res = mclapply(phenos,run_cause_on_tr1_ivs,
              phenos=phenos,G_it=NULL,GWAS_effects=GWAS_effects,GWAS_ses=GWAS_ses,
              B_distances = B_distances,
              mc.cores = 15)
        file2res[[rdata]] = cause_res
        save(file2res,file=rdata_out)
      }
    }
  }
}


############################################################################################
# Check performance
get_distance_based_metrics<-function(arc_dists){
  arc_dists = as.numeric(arc_dists)
  fdr = sum(arc_dists<0)/length(arc_dists)
  fpr = sum(arc_dists>0)/length(arc_dists)
  ntp = sum(arc_dists>0)
  return(c(fdr,fpr,ntp))
}
get_cause_perf_tables_by_fdr<-function(file2res,thr=0.1){
  all_perf = c()
  for(l in file2res){
    m = c()
    for(tr in l){
      if(length(tr)==1){next}
      m = rbind(m,tr)
    }
    m = data.frame(m,stringsAsFactors = F)
    for(j in 3:ncol(m)){
      m[[j]] = as.numeric(as.character(m[[j]]))
    }
    m$fdr = p.adjust(m[[3]],method="BY")
    m = m[!is.na(m$fdr),]
    sub_m = m[m$fdr < thr,]
    all_perf = rbind(all_perf,
                     get_distance_based_metrics(sub_m[,"X6"]))
  }
  summ_stats = rbind(
    median=apply(all_perf,2,median,na.rm=T),
    mean=apply(all_perf,2,mean,na.rm=T),
    max=apply(all_perf,2,max,na.rm=T),
    min=apply(all_perf,2,min,na.rm=T)
  )
  colnames(summ_stats) = c("FDR","TDR","N")
  return(summ_stats)
}

get_other_methods_perf<-function(res,methodreg = "^c",pleio=0.3,p2=0.01,p1=0.001){
  for_comp = res
  x = for_comp[
    all_sim_results_fdrs$prob_pleio == pleio & all_sim_results_fdrs$deg == 1.5 &
      all_sim_results_fdrs$p2 == p2,
    ]
  p1_ind = which(names(x)=="p1")
  p1_ind = p1_ind[length(p1_ind)]
  x = x[x[[p1_ind]] == p1,]
  x = x[,!grepl("lcv|egger",colnames(x))]
  x = x[,!grepl("^p",colnames(x),perl=T)]
  x = x[,!grepl("^deg",colnames(x),perl=T)]
  x = x[,!grepl("^edge",colnames(x),perl=T)]
  x = x[,grepl(methodreg,colnames(x),perl=T)]
  
  summ_stats = rbind(
    median=apply(x,2,median,na.rm=T),
    mean=apply(x,2,mean,na.rm=T),
    max=apply(x,2,max,na.rm=T),
    min=apply(x,2,min,na.rm=T)
  )
  return(summ_stats)
}

setwd("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_default/")
result_rdata_files = list.files(".")
result_rdata_files = result_rdata_files[grepl(".RData",result_rdata_files,ignore.case = T)]

cause_comparison_summary_fdr = c()
cause_comparison_summary_n = c()

# pleio=0
load("cause_results_pleio0_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.1)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.1.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0,p1=1e-05,p2=0.001)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.1)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0,p1=1e-05,p2=0.001)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n

# pleio=0.3
load("cause_results_pleio0.3_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.1)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.1.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0.3,p1=1e-05,p2=0.001)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.1)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0.3,p1=1e-05,p2=0.001)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n

# pleio=0.4
load("cause_results_pleio0.4_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.1)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.1.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0.4,p1=1e-05,p2=0.001)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.1)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0.4,p1=1e-05,p2=0.001)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n


###########
# Same as before but with FDR = 0.01
load("cause_results_pleio0_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.01)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.01.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.01)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n

# pleio=0.3
load("cause_results_pleio0.3_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.01)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.01.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0.3)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.1)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0.3)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n

###########
# Same as before but with FDR = 0.2 for CAUSE
load("cause_results_pleio0_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.2)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.1.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.2)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n


# pleio=0.3
load("cause_results_pleio0.3_deg1.5_p10.001_all_instruments.RData")
cause_fdrs = get_cause_perf_tables_by_fdr(file2res,0.2)[,1]
load("../simulations_uniqueiv_minIV3/simulation_summ_stats_FDR0.1.RData")
other_fdrs = get_other_methods_perf(all_sim_results_fdrs,"",0.3)
all_fdrs = cbind(other_fdrs,cause_fdrs)
colnames(all_fdrs)[ncol(all_fdrs)] = "cause"
all_fdrs
cause_n = get_cause_perf_tables_by_fdr(file2res,0.2)[,3]
other_n = get_other_methods_perf(all_sim_results_preds,"",0.3)
all_n = cbind(other_n,cause_n)
colnames(all_n)[ncol(all_n)] = "cause"
all_n




