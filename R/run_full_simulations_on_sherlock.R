# This script runs the simulations using Stanford's Sherlock (a cluster)
# We have one version that runs all methods on the same simulated graphs,
# however, this is redundant for the edge sep test methods that do not use
# a p2 parameter - we therefore have two versions of the simulations in this 
# script: one for all methods and another for the edge sep tests

###################################################################################
### Helper functions to run the simulations

get_sh_prefix<-function(err="",log="",time="0:45:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=", time,sep=""),
      "#SBATCH --partition=euan,mrivas,owners,normal",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --mem=2000",
      # "#SBATCH -x sh-113-15",
      "#SBATCH --gpus-per-task=0",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load R"
    )
  )
}

print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}

exec_cmd_on_sherlock<-function(cmd,jobname,out_path){
  err_path = paste(out_path,jobname,".err",sep="")
  log_path = paste(out_path,jobname,".log",sep="")
  curr_cmd = paste("Rscript",cmd)
  curr_sh_file = paste(out_path,jobname,".sh",sep="")
  sh_prefix = get_sh_prefix(err_path,log_path)
  print_sh_file(curr_sh_file,sh_prefix,curr_cmd)
  system(paste("sbatch",curr_sh_file,'&'))
}

###################################################################################
# Set the WD and number of repeats
reps = 40
WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_uniqueiv_minIV3/"
try(system(paste("mkdir",WD)))
MAX_JOBS = 500
# Set the simulation parameters
tested_p1 = c(1e-02,1e-03,1e-04,1e-05)
tested_p2_factors = c(1,10,100)
tested_pleio_levels = c(0,0.1,0.2,0.3,0.4)
tested_degrees = c(1,1.25,1.5,1.75,2)

###################################################################################
###################################################################################
###################################################################################
# Analysis of all methods together
reps*length(tested_p1)*length(tested_p2_factors)*length(tested_pleio_levels)*length(tested_degrees)

# Run the simulations
for(p1 in tested_p1){
  print(paste("p1",p1))
  for(p2_f in tested_p2_factors){
    p2 = p1*p2_f
    print(paste("p2",p2))
    if(p2 < p1){next}
    if(p2 > 0.1){next}
    for(pleio_p in tested_pleio_levels){
      print(paste("pleio_p",pleio_p))
      for(deg in tested_degrees){
        print(paste("deg",deg))
        curr_folder = paste(WD,"deg",deg,"_pleio",pleio_p,"_p1",p1,"_p2",p2,"/",sep="")
        system(paste("mkdir",curr_folder))
        
        # check current jobs and wait if there are too much
        job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
        num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
        while(num_active_jobs > MAX_JOBS){
          Sys.sleep(5)
          job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
          num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
        }
        
        for(i in 1:reps){
          
          curr_out_file = paste(curr_folder,"sim_rep",i,".RData",sep="")
          if(file.exists(curr_out_file)){
            # print(paste("results for rep",i,"exists,skipping"))
            next
          }
          cmd = paste(
            "~/repos/cGAUGE/R/full_causal_graph_simulations.R",
            "--deg",deg,
            "--probPleio",pleio_p,
            "--cgaugeMode 1",
            "--p1",p1,
            "--p2",p2,
            "--out",curr_out_file
          )
          print(i)
          exec_cmd_on_sherlock(cmd,jobname = paste("sim_rep",i,sep=""),out_path = curr_folder)
        }
      }
    }
  }
}

# Go over the results: read and save the simulation results - all methods
FDR = 0.1
FDR_method = "BY"
is_causal<-function(dists){
  return(dists>0 )
}

all_sim_results_errs = c()
all_sim_results_preds = c()
for(p1 in tested_p1){
  print(paste("p1",p1))
  for(p2_f in tested_p2_factors){
    p2 = p1*p2_f
    print(paste("p2",p2))
    if(p2 >0.1){next}
    if(p2<p1){next}
    print(paste("p2",p2))
    for(pleio_p in tested_pleio_levels){
      print(paste("pleio_p",pleio_p))
      for(deg in tested_degrees){
        print(paste("deg",deg))
        curr_folder = paste(WD,"deg",deg,"_pleio",pleio_p,"_p1",p1,"_p2",p2,"/",sep="")
        curr_files = list.files(curr_folder)
        curr_files = curr_files[grepl("RData$",curr_files)]
        print(length(curr_files))
        if(length(curr_files)<reps){
          print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
          print(curr_folder)
          print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
          if(length(curr_files)==0){next}
        }
        
        # Measure how many of the "-1"'s are marked as TRUE - obvious errors
        method2false_discoveries = c()
        method2num_discoveries = c()
        for(curr_out_file in curr_files){
          curr_out_file = paste(curr_folder,curr_out_file,sep="")
          load(curr_out_file)
          
          num_not_causal = standard_mr_results$Egger$KnownDistance==-1
          errs = c();preds = c();edgesep_res = c()
          egger_qs = p.adjust(standard_mr_results$Egger$p,method=FDR_method)
          errs["egger"] = sum(egger_qs < FDR & !is.na(egger_qs) & 
                                !is_causal(standard_mr_results$Egger$KnownDistance),na.rm = T)
          preds["egger"] = sum(egger_qs < FDR,na.rm = T)
          
          egger_qs = p.adjust(cgauge_mr_results$Egger$p,method=FDR_method)
          errs["c-egger"] = sum(egger_qs < FDR & !is.na(egger_qs) & 
                                  !is_causal(cgauge_mr_results$Egger$KnownDistance),na.rm = T)
          preds["c-egger"] = sum(egger_qs < FDR,na.rm = T)
          
          ivw_qs = p.adjust(standard_mr_results$IVW$p,method=FDR_method)
          errs["ivw"] = sum(ivw_qs < FDR & !is.na(ivw_qs) & 
                              !is_causal(standard_mr_results$IVW$KnownDistance),na.rm = T)
          preds["ivw"] = sum(ivw_qs < FDR,na.rm = T)
          
          ivw_qs = p.adjust(cgauge_mr_results$IVW$p,method=FDR_method)
          errs["c-ivw"] = sum(ivw_qs < FDR & !is.na(ivw_qs) & 
                                !is_causal(cgauge_mr_results$IVW$KnownDistance),na.rm = T)
          preds["c-ivw"] = sum(ivw_qs < FDR,na.rm = T)
          
          presso_q = p.adjust(standard_mr_results$MRPRESSO$`P-value`,method=FDR_method)
          errs["mrpresso"] = sum(presso_q < FDR & !is.na(presso_q) &
                                   !is_causal(standard_mr_results$MRPRESSO$KnownDistance),na.rm = T)
          preds["mrpresso"] = sum(presso_q < FDR,na.rm = T)
          
          presso_q = p.adjust(cgauge_mr_results$MRPRESSO$`P-value`,method=FDR_method)
          errs["c-mrpresso"] = sum(presso_q < FDR & !is.na(presso_q) &
                                     !is_causal(cgauge_mr_results$MRPRESSO$KnownDistance),na.rm = T)
          preds["c-mrpresso"] = sum(presso_q < FDR,na.rm = T)
          
          # LCV
          lcv_q = p.adjust(standard_mr_results$LCV$P,method=FDR_method)
          errs["lcv"] = sum(lcv_q < FDR & !is.na(lcv_q) &
                              !is_causal(standard_mr_results$LCV$KnownDistance),na.rm = T)
          preds["lcv"] = sum(lcv_q < FDR,na.rm = T)
          
          # EdgeSep
          edge_sep_results = edge_sep_results[edge_sep_results$num_edgesep > 2,]
          errs["edge_sep"] = sum(!is_causal(edge_sep_results$real_distance))
          preds["edge_sep"] = nrow(edge_sep_results)
          
          # EdgeSepTest 1
          edge_sep_results_statTest1 = 
            edge_sep_results_statTest1[
              p.adjust(edge_sep_results_statTest1$`pval:trait1->trait2`,method=FDR_method)<FDR,]
          errs["edge_sep_test1"] = sum(!is_causal(edge_sep_results_statTest1$KnownDistance))
          preds["edge_sep_test1"] = nrow(edge_sep_results_statTest1)
          
          # EdgeSepTest 2
          edge_sep_results_statTest2 = 
            edge_sep_results_statTest2[
              p.adjust(edge_sep_results_statTest2$`pval:trait1->trait2`,method=FDR_method)<FDR,]
          errs["edge_sep_test2"] = sum(!is_causal(edge_sep_results_statTest2$KnownDistance))
          preds["edge_sep_test2"] = nrow(edge_sep_results_statTest2)
          
          method2num_discoveries = rbind(method2num_discoveries,preds)
          method2false_discoveries = rbind(method2false_discoveries,errs)
          
        }
        method2num_discoveries = as.data.frame(method2num_discoveries)
        method2false_discoveries = as.data.frame(method2false_discoveries)
        
        method2false_discoveries$p1 = p1
        method2false_discoveries$p2 = p2
        method2false_discoveries$prob_pleio = pleio_p
        method2false_discoveries$deg = deg
        
        method2num_discoveries$p1 = p1
        method2num_discoveries$p2 = p2
        method2num_discoveries$prob_pleio = pleio_p
        method2num_discoveries$deg = deg
        
        all_sim_results_errs = rbind(all_sim_results_errs,method2false_discoveries)
        all_sim_results_preds = rbind(all_sim_results_preds,method2num_discoveries)
      }
    }
  }
}

###############################################################
# Summarize and save
colinds = 1:10

mean_errs = aggregate(all_sim_results_errs[,colinds],
                      by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                              deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                      FUN=mean)

sd_errs = aggregate(all_sim_results_errs[,colinds],
                    by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                            deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=sd)
mean_num_discoveries = aggregate(all_sim_results_preds[,colinds],
                                 by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                                         deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                                 FUN=mean)


all_sim_results_fdrs = (all_sim_results_errs/all_sim_results_preds)[,colinds]
all_sim_results_fdrs[is.nan(as.matrix(all_sim_results_fdrs))] = 0
all_sim_results_fdrs = 
  cbind(all_sim_results_fdrs,all_sim_results_errs[,c("p1","p2","deg","prob_pleio")])
mean_fdrs = aggregate(all_sim_results_fdrs,
                      by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                              deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                      FUN=mean,na.rm=T)

sd_fdrs = aggregate(all_sim_results_fdrs,
                    by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                            deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=sd,na.rm=T)

save(
  all_sim_results_errs,all_sim_results_preds,
  all_sim_results_fdrs,
  mean_errs,sd_errs,mean_num_discoveries,
  mean_fdrs,sd_fdrs,
  file = paste(WD,"/simulation_summ_stats.RData",sep="")
)


###################################################################################
###################################################################################
###################################################################################
# Analysis of the edge sep tests
WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_ms_test_minMu1_3/"
try(system(paste("mkdir",WD)))
MAX_JOBS = 400
# Set the simulation parameters
tested_p1 = c(1e-02,1e-03,1e-04,1e-05)
tested_p2_factors = c(1,10,100)
tested_pleio_levels = c(0,0.1,0.2,0.3,0.4)
tested_degrees = c(1,1.25,1.5,1.75,2)

reps = 40

reps*length(tested_p1)*length(tested_pleio_levels)*length(tested_degrees)
# Run with the edge sep tests only (no need for p2)
for(p1 in tested_p1){
  print(paste("p1",p1))
  for(pleio_p in tested_pleio_levels){
    print(paste("pleio_p",pleio_p))
    for(deg in tested_degrees){
      print(paste("deg",deg))
      curr_folder = paste(WD,"deg",deg,"_pleio",pleio_p,"_p1",p1,"/",sep="")
      system(paste("mkdir",curr_folder))
        
      # check current jobs and wait if there are too much
      job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
      num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
      while(num_active_jobs > MAX_JOBS){
        Sys.sleep(5)
        job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
        num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
      }
        
      for(i in 1:reps){
          
        curr_out_file = paste(curr_folder,"sim_rep",i,".RData",sep="")
        if(file.exists(curr_out_file)){next}
        cmd = paste(
          "~/repos/cGAUGE/R/full_causal_graph_simulations.R",
          "--deg",deg,
          "--probPleio",pleio_p,
          "--cgaugeMode 1",
          "--edgeSepRun 1",
          "--p1",p1,
          "--p2",0.01, # a dummy input, not used for the edgesep tests
          "--out",curr_out_file
        )
        print(i)
        exec_cmd_on_sherlock(cmd,jobname = paste("sim_rep",i,sep=""),out_path = curr_folder)
      }
    }
  }
}

# Read and save the simulation results 
FDR = 0.01
FDR_method = "BY"
is_causal<-function(dists){
  return(dists>0 )
}

all_sim_results_errs = c()
all_sim_results_preds = c()
for(p1 in tested_p1){
  print(paste("p1",p1))
  for(pleio_p in tested_pleio_levels){
    print(paste("pleio_p",pleio_p))
    for(deg in tested_degrees){
      print(paste("deg",deg))
      curr_folder = paste(WD,"deg",deg,"_pleio",pleio_p,"_p1",p1,"/",sep="")
      curr_files = list.files(curr_folder)
      curr_files = curr_files[grepl("RData$",curr_files)]
      print(length(curr_files))
      if(length(curr_files)<reps){
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(curr_folder)
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        if(length(curr_files)==0){next}
      }
        
      # Measure how many of the "-1"'s are marked as TRUE - obvious errors
      method2false_discoveries = c()
      method2num_discoveries = c()
      for(curr_out_file in curr_files){
        curr_out_file = paste(curr_folder,curr_out_file,sep="")
        load(curr_out_file)
        errs = c();preds = c();edgesep_res = c()
          
        # EdgeSepTest 1
        edge_sep_results_statTest1 = 
          edge_sep_results_statTest1[
            p.adjust(edge_sep_results_statTest1$`pval:trait1->trait2`,method=FDR_method)<FDR,]
        errs["edge_sep_test1"] = sum(!is_causal(edge_sep_results_statTest1$KnownDistance))
        preds["edge_sep_test1"] = nrow(edge_sep_results_statTest1)
        
        # # EdgeSepTest 2
        # edge_sep_results_statTest2 =
        #   edge_sep_results_statTest2[
        #     p.adjust(edge_sep_results_statTest2$`pval:trait1->trait2`,method=FDR_method)<FDR,]
        # errs["edge_sep_test2"] = sum(!is_causal(edge_sep_results_statTest2$KnownDistance))
        # preds["edge_sep_test2"] = nrow(edge_sep_results_statTest2)
        
        method2num_discoveries = rbind(method2num_discoveries,preds)
        method2false_discoveries = rbind(method2false_discoveries,errs)
          
      }
      method2num_discoveries = as.data.frame(method2num_discoveries)
      method2false_discoveries = as.data.frame(method2false_discoveries)
      method2false_discoveries$p1 = p1
      method2false_discoveries$prob_pleio = pleio_p
      method2false_discoveries$deg = deg
        
      method2num_discoveries$p1 = p1
      method2num_discoveries$prob_pleio = pleio_p
      method2num_discoveries$deg = deg
        
      all_sim_results_errs = rbind(all_sim_results_errs,method2false_discoveries)
      all_sim_results_preds = rbind(all_sim_results_preds,method2num_discoveries)
    }
  }
}

###############################################################
# Summarize and save
colinds = 1:2 

mean_errs = aggregate(all_sim_results_errs[,colinds],
                    by=list(p1=all_sim_results_errs$p1,
                    deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=mean)
sd_errs = aggregate(all_sim_results_errs[,colinds],
                    by=list(p1=all_sim_results_errs$p1,
                    deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=sd)
mean_num_discoveries = aggregate(all_sim_results_preds[,colinds],
                    by=list(p1=all_sim_results_errs$p1,
                    deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=mean)

all_sim_results_fdrs = (all_sim_results_errs/all_sim_results_preds)[,colinds]
all_sim_results_fdrs[is.nan(as.matrix(all_sim_results_fdrs))] = 0
all_sim_results_fdrs = 
  cbind(all_sim_results_fdrs,all_sim_results_errs[,c("p1","deg","prob_pleio")])

mean_fdrs = aggregate(all_sim_results_fdrs,
                      by=list(p1=all_sim_results_errs$p1,
                      deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                      FUN=mean,na.rm=T)

sd_fdrs = aggregate(all_sim_results_fdrs,
                    by=list(p1=all_sim_results_errs$p1,
                    deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=sd,na.rm=T)

save(
  all_sim_results_errs,all_sim_results_preds,
  all_sim_results_fdrs,
  mean_errs,sd_errs,mean_num_discoveries,
  mean_fdrs,sd_fdrs,
  file = paste(WD,"/simulation_summ_stats_",FDR,"FDR.RData",sep="")
)

##############################################################
##############################################################
# Add the pi1 estimates for pairs
##############################################################
##############################################################
# WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_uniqueivs/"
# WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_default/"
WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_uniqueiv2/"
setwd(WD)
# Go over all RData files in WD and add their pi1_estimates
wd_folders = list.files(WD,full.names = T)
for(folder in wd_folders){
  if(!dir.exists(folder)){next}
  print(paste("***********",folder))
  folder_files = list.files(folder,full.names = T)
  folder_files = folder_files[grepl(".RData",folder_files)]
  folder_files = folder_files[grepl("sim_rep",folder_files)]
  for(rdfile in folder_files){
    pi1_estimates = NULL
    load(rdfile)
    if("pi1_estimates" %in% ls() && !is.null(pi1_estimates) && !is.null(dim(pi1_estimates))){next}
    i = strsplit(rdfile,split = "sim_rep")[[1]][2]
    i = gsub(".RData","",i)
    p1 = strsplit(rdfile,split = "_p1")[[1]][2]
    p1 = strsplit(p1,split = "_")[[1]][1]
    print(paste(i,p1))
    cmd = paste(
      "~/repos/cGAUGE/R/simulations_add_pi1_estimates.R",
      "--p1",p1,
      "--rdata_file",rdfile
    )
    exec_cmd_on_sherlock(cmd,jobname = paste("sim_rep",i,"_pi1",sep=""),out_path = folder)
  }
  
  # check current jobs and wait if there are too much
  job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
  num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
  while(num_active_jobs > MAX_JOBS){
    Sys.sleep(5)
    job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
    num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
  }
}

# The commands above create the log and err files in the WD - delete them
# (the reason is that we did not add a "/" after the folder name)
system("rm *.err")
system("rm *.log")
system("rm *.sh")

# Read and summarize the results
analyze_pi1_estimates<-function(x,thrs = seq(0.1,0.9,by=0.1),pi1col=2){
  nonedges = x[,3]==-1
  fdrs = c()
  num_discoveries = c()
  fprs = c()
  for(thr in thrs){
    s_thr = as.character(thr)
    discoveries = x[,pi1col] > thr
    fdrs[s_thr] = sum(nonedges[discoveries],na.rm = T)/sum(discoveries,na.rm = T)
    num_discoveries[s_thr] = sum(discoveries,na.rm = T)
    fprs[s_thr] = sum(nonedges[!discoveries],na.rm = T)/sum(!discoveries,na.rm = T)
  }
  return(list(fdrs=fdrs,num_discoveries = num_discoveries,fprs=fprs))
}

format_rep_results<-function(x,l){
  x = as.data.frame(x)
  for(y in names(l)){
    x[[y]] = l[[y]]
  }
  return(x)
}
  
# Set the simulation parameters
tested_p1 = c(1e-02,1e-03,1e-04,1e-05)
tested_p2_factors = c(1,10,100)
tested_pleio_levels = c(0,0.1,0.2,0.3,0.4)
tested_degrees = c(1,1.25,1.5,1.75,2)

# Read the results
all_sim_results_fdrs = c()
all_sim_results_preds  = c()
all_sim_results_fprs = c()
for(p1 in tested_p1){
  print(paste("p1",p1))
  for(p2_f in tested_p2_factors){
    p2 = p1*p2_f
    print(paste("p2",p2))
    if(p2 < p1){next}
    if(p2 > 0.1){next}
    for(pleio_p in tested_pleio_levels){
      print(paste("pleio_p",pleio_p))
      for(deg in tested_degrees){
        print(paste("deg",deg))
        curr_folder = paste(WD,"deg",deg,"_pleio",pleio_p,"_p1",p1,"_p2",p2,"/",sep="")
        currfdrs = c();currNs = c();currfprs = c()
        for(i in 1:reps){
          pi1_estimates = NULL
          curr_out_file = paste(curr_folder,"sim_rep",i,".RData",sep="")
          try({load(curr_out_file)})
          if(is.null(pi1_estimates)||length(pi1_estimates)==0){
            print(paste("missing pi1 in:",curr_out_file))
            next
          }
          pi1_res2 = analyze_pi1_estimates(pi1_estimates,pi1col=2)
          pi1_res1 = analyze_pi1_estimates(pi1_estimates,pi1col=1)
          pi1_res2 = lapply(pi1_res2,function(x){names(x)=paste("our",names(x));x})
          pi1_res1 = lapply(pi1_res1,function(x){names(x)=paste("raw",names(x));x})
          currfdrs = rbind(currfdrs,c(pi1_res1$fdrs,pi1_res2$fdrs))
          currNs = rbind(currNs,c(pi1_res1$num_discoveries,pi1_res2$num_discoveries))
          currfprs = rbind(currfprs,c(pi1_res1$fprs,pi1_res2$fprs))
        }
        curr_params = list(p1=p1,p2=p2,deg=deg,pleio_p=pleio_p)
        currfdrs = format_rep_results(currfdrs,curr_params)
        currNs = format_rep_results(currNs,curr_params)
        currfprs = format_rep_results(currfprs,curr_params)
        all_sim_results_fdrs = rbind(all_sim_results_fdrs,currfdrs)
        all_sim_results_fprs = rbind(all_sim_results_fprs,currfprs)
        all_sim_results_preds = rbind(all_sim_results_preds,currNs)
      }
    }
  }
}

mean_fdrs = aggregate(all_sim_results_fdrs,
     by=list(p1=all_sim_results_fdrs$p1,p2 = all_sim_results_fdrs$p2,
     deg = all_sim_results_fdrs$deg,prob_pleio = all_sim_results_fdrs$pleio_p),
     FUN=mean,na.rm=T)

sd_fdrs = aggregate(all_sim_results_fdrs,
      by=list(p1=all_sim_results_fdrs$p1,p2 = all_sim_results_fdrs$p2,
      deg = all_sim_results_fdrs$deg,prob_pleio = all_sim_results_fdrs$pleio_p),
      FUN=sd,na.rm=T)

mean_fprs = aggregate(all_sim_results_fprs,
      by=list(p1=all_sim_results_fprs$p1,p2 = all_sim_results_fprs$p2,
      deg = all_sim_results_fprs$deg,prob_pleio = all_sim_results_fprs$pleio_p),
      FUN=mean,na.rm=T)

sd_fprs = aggregate(all_sim_results_fprs,
       by=list(p1=all_sim_results_fprs$p1,p2 = all_sim_results_fprs$p2,
       deg = all_sim_results_fprs$deg,prob_pleio = all_sim_results_fprs$pleio_p),
       FUN=sd,na.rm=T)

mean_num_discoveries = aggregate(all_sim_results_preds,
       by=list(p1=all_sim_results_preds$p1,p2 = all_sim_results_preds$p2,
       deg = all_sim_results_preds$deg,prob_pleio = all_sim_results_preds$pleio_p),
       FUN=mean,na.rm=T)

sd_num_discoveries = aggregate(all_sim_results_preds,
       by=list(p1=all_sim_results_preds$p1,p2 = all_sim_results_preds$p2,
       deg = all_sim_results_preds$deg,prob_pleio = all_sim_results_preds$pleio_p),
       FUN=sd,na.rm=T)

save(
  all_sim_results_fdrs,all_sim_results_fprs,all_sim_results_preds,
  mean_num_discoveries,sd_num_discoveries,
  mean_fdrs,sd_fdrs,
  mean_fprs,sd_fprs,
  file = paste(WD,"/simulation_pi1_summ_stats.RData",sep="")
)

##############################################################
##############################################################
# Write down the supplementary tables
##############################################################
##############################################################

setwd("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/")
library("xlsx",lib.loc = "~/R/packages3.5")

supp_readme = c()

load("simulations_default/simulation_summ_stats.RData")
mean_fdrs_thm21 = mean_fdrs
sd_fdrs_thm21 = sd_fdrs
mean_num_discoveries_thm21 = mean_num_discoveries
load("simulations_uniqueiv2/simulation_summ_stats.RData")
mean_fdrs_thm22 = mean_fdrs
sd_fdrs_thm22 = sd_fdrs
mean_num_discoveries_thm22 = mean_num_discoveries

all(mean_fdrs_thm21[,1:4]==mean_fdrs_thm22[,1:4])
mr_fdr_simulations_supp_table = mean_fdrs_thm21[,1:4]
for(j in c(5,7,9)){
  x1 = format(mean_fdrs_thm21[,j],digits=3)
  x1 = paste(x1," (",format(sd_fdrs_thm21[,j],digits=3),")",sep="")
  x21 = format(mean_fdrs_thm21[,j+1],digits=3)
  x21 = paste(x21," (",format(sd_fdrs_thm21[,j+1],digits=3),")",sep="")
  x22 = format(mean_fdrs_thm22[,j+1],digits=3)
  x22 = paste(x22," (",format(sd_fdrs_thm22[,j+1],digits=3),")",sep="")
  m = cbind(x1,x21,x22)
  currname = colnames(mean_fdrs_thm21)[j]
  colnames(m) = c(
    currname,
    paste(currname,"ImpIV",sep=","),
    paste(currname,"UniqueIV",sep=",")
  )
  mr_fdr_simulations_supp_table = cbind(mr_fdr_simulations_supp_table,m)
}

# Add the naive edge sep
j=12
exsep = format(mean_fdrs_thm21[,j],digits=3)
exsep = paste(exsep," (",format(sd_fdrs_thm21[,j],digits=3),")",sep="")
mr_fdr_simulations_supp_table = cbind(mr_fdr_simulations_supp_table,exsep)
colnames(mr_fdr_simulations_supp_table)[ncol(mr_fdr_simulations_supp_table)] = "ExSep, Naive"
write.xlsx2(mr_fdr_simulations_supp_table,file="./supp_tables/Supplementary_Tables.xlsx",
            sheetName = "ST1",row.names=F)
supp_readme = c(supp_readme,
                "ST1: Empirical FDRs of the MR and the naive ExSep methods")


mr_N_simulations_supp_table = mean_num_discoveries_thm21[,1:4]
for(j in c(5,7,9)){
  x1 = format(mean_num_discoveries_thm21[,j],digits=3)
  x21 = format(mean_num_discoveries_thm21[,j+1],digits=3)
  x22 = format(mean_num_discoveries_thm22[,j+1],digits=3)
  m = cbind(x1,x21,x22)
  currname = colnames(mean_num_discoveries_thm21)[j]
  colnames(m) = c(
    currname,
    paste(currname,"ImpIV",sep=","),
    paste(currname,"UniqueIV",sep=",")
  )
  mr_N_simulations_supp_table = cbind(mr_N_simulations_supp_table,m)
}

# Add the naive edge sep
j=12
exsep = format(mean_num_discoveries_thm21[,j],digits=3)
mr_N_simulations_supp_table = cbind(mr_N_simulations_supp_table,exsep)
colnames(mr_N_simulations_supp_table)[ncol(mr_N_simulations_supp_table)] = "ExSep, Naive"
write.xlsx2(mr_N_simulations_supp_table,file="./supp_tables/Supplementary_Tables.xlsx",
            sheetName = "ST2",row.names=F,append=T)
supp_readme = c(supp_readme,
                "ST2: Number of discoveries of the MR and the naive ExSep methods")

load("./simulations_edgesep/simulation_summ_stats.RData")
edesep_fdr_simulations_supp_table = mean_fdrs[,1:3]
for(j in c(4,5)){
  x1 = format(mean_fdrs[,j],digits=3)
  x1 = paste(x1," (",format(sd_fdrs[,j],digits=3),")",sep="")
  edesep_fdr_simulations_supp_table = cbind(edesep_fdr_simulations_supp_table,x1)
}
colnames(edesep_fdr_simulations_supp_table)[4:5] = c("MS test","TDR test")

edgesep_N_simulations_supp_table = mean_num_discoveries[,1:3]
for(j in c(4,5)){
  x1 = format(mean_num_discoveries[,j],digits=3)
  edgesep_N_simulations_supp_table = cbind(edgesep_N_simulations_supp_table,x1)
}
colnames(edgesep_N_simulations_supp_table)[4:5] = c("MS test","TDR test")

write.xlsx2(edesep_fdr_simulations_supp_table,
            file="./supp_tables/Supplementary_Tables.xlsx",
            sheetName = "ST3",row.names=F,append=T)
supp_readme = c(supp_readme,
                "ST3: Empirical FDRs of the ExSep tests")

write.xlsx2(edgesep_N_simulations_supp_table,
            file="./supp_tables/Supplementary_Tables.xlsx",
            sheetName = "ST4",row.names=F,append=T)
supp_readme = c(supp_readme,
                "ST4: Number of discoveries of the ExSep tests")

# Pi1 results
load("simulations_default/simulation_pi1_summ_stats.RData")
mean_fdrs_thm21 = mean_fdrs
sd_fdrs_thm21 = sd_fdrs
mean_fprs_thm21 = mean_fprs
sd_fprs_thm21 = sd_fprs
mean_num_discoveries_thm21 = mean_num_discoveries
load("simulations_uniqueiv2/simulation_pi1_summ_stats.RData")
mean_fdrs_thm22 = mean_fdrs
sd_fdrs_thm22 = sd_fdrs
mean_fprs_thm22 = mean_fprs
sd_fprs_thm22 = sd_fprs
mean_num_discoveries_thm22 = mean_num_discoveries
colnames(mean_fdrs_thm21) = gsub(
  "raw ",replacement = "FDR,raw,pi1>",colnames(mean_fdrs_thm21)
)
colnames(mean_fdrs_thm21) = gsub(
  "our ",replacement = "FDR,ImpIV,pi1>",colnames(mean_fdrs_thm21)
)
colnames(mean_fprs_thm21) = gsub(
  "raw ",replacement = "FPR,raw,pi1<",colnames(mean_fprs_thm21)
)
colnames(mean_fprs_thm21) = gsub(
  "our ",replacement = "FPR,ImpIV,pi1<",colnames(mean_fprs_thm21)
)
colnames(mean_num_discoveries_thm21) = gsub(
  "raw ",replacement = "N,raw,pi1>",colnames(mean_num_discoveries_thm21)
)
colnames(mean_num_discoveries_thm21) = gsub(
  "our ",replacement = "N,ImpIV,pi1<",colnames(mean_num_discoveries_thm21)
)
colnames(mean_fdrs_thm22) = gsub(
  "our ",replacement = "FDR,ImpIV,pi1>",colnames(mean_fdrs_thm22)
)
colnames(mean_fprs_thm22) = gsub(
  "our ",replacement = "FPR,ImpIV,pi1<",colnames(mean_fprs_thm22)
)
colnames(mean_num_discoveries_thm22) = gsub(
  "our ",replacement = "N,UniqueIV,pi1>",colnames(mean_num_discoveries_thm22)
)

pi1_fdr_fpr_supp_table = mean_fdrs_thm21[,1:4]
for(j in c(5:22)){
  x1 = format(mean_fdrs_thm21[,j],digits=3)
  x1 = paste(x1," (",format(sd_fdrs_thm21[,j],digits=3),")",sep="")
  x2 = format(mean_fprs_thm21[,j],digits=3)
  x2 = paste(x2," (",format(sd_fprs_thm21[,j],digits=3),")",sep="")
  m = cbind(x1,x2)
  colnames(m) = c(
    colnames(mean_fdrs_thm21)[j],
    colnames(mean_fprs_thm21)[j]
  )
  pi1_fdr_fpr_supp_table = cbind(pi1_fdr_fpr_supp_table,m)
}
for(j in c(14:22)){
  x1 = format(mean_fdrs_thm22[,j],digits=3)
  x1 = paste(x1," (",format(sd_fdrs_thm22[,j],digits=3),")",sep="")
  x2 = format(mean_fprs_thm22[,j],digits=3)
  x2 = paste(x2," (",format(sd_fprs_thm22[,j],digits=3),")",sep="")
  m = cbind(x1,x2)
  colnames(m) = c(
    colnames(mean_fdrs_thm22)[j],
    colnames(mean_fprs_thm22)[j]
  )
  pi1_fdr_fpr_supp_table = cbind(pi1_fdr_fpr_supp_table,m)
}

# Add the naive edge sep
write.xlsx2(pi1_fdr_fpr_supp_table,file="./supp_tables/Supplementary_Tables.xlsx",
            sheetName = "ST5",row.names=F,append=T)
supp_readme = c(supp_readme,
                "ST5: Empirical FDRs and FPRs using the pi1 estimates")

pi1_N_simulations_supp_table = mean_num_discoveries_thm21[,1:4]
for(j in c(5:22)){
  x1 = format(mean_num_discoveries_thm21[,j],digits=3)
  pi1_N_simulations_supp_table = cbind(pi1_N_simulations_supp_table,x1)
  colnames(pi1_N_simulations_supp_table)[ncol(pi1_N_simulations_supp_table)] = 
    colnames(mean_num_discoveries_thm21)[j]
}
for(j in c(14:22)){
  x1 = format(mean_num_discoveries_thm22[,j],digits=3)
  pi1_N_simulations_supp_table = cbind(pi1_N_simulations_supp_table,x1)
  colnames(pi1_N_simulations_supp_table)[ncol(pi1_N_simulations_supp_table)] = 
    colnames(mean_num_discoveries_thm22)[j]
}

write.xlsx2(pi1_N_simulations_supp_table,
            file="./supp_tables/Supplementary_Tables.xlsx",
            sheetName = "ST6",row.names=F,append=T)
supp_readme = c(supp_readme,
                "ST6: Number of discoveries using the pi1 estimates (pi1>thr)")

write(supp_readme,file = "./supp_tables/supp_readme.txt")

##############################################################
##############################################################
# Locally - figures
##############################################################
##############################################################
# install.packages("fmsb")
setwd("~/Desktop/causal_inference_projects/ms3/")
library(fmsb)

deg = 1.5
p1 = 1e-05
p2 = 0.001

#######
# Simulations using Thm 22 - better WC
#######

method2col = c(hcl.colors(5)[1:2],heat.colors(5)[1:2],
               rainbow(5)[1:2],"black","gray","blue")
names(method2col) = c(
  "mrpresso","c-mrpresso","ivw","c-ivw","egger","c-egger",
  "edge_sep","edge_sep_test1","edge_sep_test2"
)

load("simulations_uniqueiv/simulation_summ_stats.RData")
par(mfrow=c(1,2),mar=c(0,1,4,1),xpd=TRUE)
# (A) IVW, MR-PRESSO: Num discoveries
resultsdf = mean_num_discoveries
ndigits = 0
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","mrpresso","c-mrpresso","ivw","c-ivw")]
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,pcol=cols,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,rownames(df1)[-c(1:2)],fill = cols,ncol = 2,border = F)
text(0,0,"N",cex = 1.5,font = 2)

# (B) IVW, MR-PRESSO: FDRs
resultsdf = mean_fdrs
ndigits = 2
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","mrpresso","c-mrpresso","ivw","c-ivw")]
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,pcol=cols,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,rownames(df1)[-c(1:2)],fill = cols,ncol = 2,border = F)
text(0,0,"FDR",cex = 1.2,font = 2)

# plot the legend separately
dev.off()
par(mar = c(5,5,5,5))
plot(c(2,2))
legend(x="top",c("MRPRESSO","UnIV+MRPRESSO","IVW","UnIV+IVW"),
       fill = method2col[1:4],ncol = 1,border = F,cex=2)

#######
# Simulations using Thm 21 - lower WC performance
#######

method2col = c(hcl.colors(5)[c(1,4)],heat.colors(5)[c(1,4)],
               rainbow(5)[1:2],"black","gray","blue")
names(method2col) = c(
  "mrpresso","c-mrpresso","ivw","c-ivw","egger","c-egger",
  "edge_sep","edge_sep_test1","edge_sep_test2"
)

load("simulations_def/simulation_summ_stats.RData")
par(mfrow=c(1,2),mar=c(0,1,4,1),xpd=TRUE)
# (A) IVW, MR-PRESSO: Num discoveries
resultsdf = mean_num_discoveries
ndigits = 0
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","mrpresso","c-mrpresso","ivw","c-ivw")]
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,pcol=cols,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,rownames(df1)[-c(1:2)],fill = cols,ncol = 2,border = F)
text(0,0,"N",cex = 1.5,font = 2)

# (B) IVW, MR-PRESSO: FDRs
resultsdf = mean_fdrs
ndigits = 2
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","mrpresso","c-mrpresso","ivw","c-ivw")]
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,pcol=cols,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,rownames(df1)[-c(1:2)],fill = cols,ncol = 2,border = F)
text(0,0,"FDR",cex = 1.2,font=2)

# plot the legend separately
dev.off()
par(mar = c(5,5,5,5))
plot(c(2,2))
legend(x="top",c("MRPRESSO","ImpIV+MRPRESSO","IVW","ImpIV+IVW"),
       fill = method2col[1:4],ncol = 1,border = F,cex=2)

#######
# Edge sep results: 
# Take the naive analysis from the results above
# Take the edge sep tests from their simulations
#######
load("simulations_def/simulation_summ_stats.RData")
naive_mean_fdrs = mean_fdrs
naive_mean_num_discoveries = mean_num_discoveries
load("simulations_edgesep/simulation_summ_stats.RData")
par(mfrow=c(1,2),mar=c(0,1,4,1),xpd=TRUE)

# Num discoveries
ndigits = 0
resultsdf = naive_mean_num_discoveries
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","edge_sep")]
resultsdf = mean_num_discoveries
inds = resultsdf$deg==deg & resultsdf$p1 == p1 
df2 = resultsdf[inds,c("prob_pleio","edge_sep_test1","edge_sep_test2")]
if(all(df1[,1]==df2[,1])){
  df1 = cbind(df1,df2[,-1])
}
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,pcol=cols,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,c("Naive count","EM test","TDR test"),fill = cols,ncol = 2,border=F)
text(0,0,"N",cex = 1.5,font = 2)

# (D) IVW, MR-PRESSO: FDRs
ndigits = 2
resultsdf = naive_mean_fdrs
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","edge_sep")]
resultsdf = mean_fdrs
inds = resultsdf$deg==deg & resultsdf$p1 == p1 
df2 = resultsdf[inds,c("prob_pleio","edge_sep_test1","edge_sep_test2")]
if(all(df1[,1]==df2[,1])){
  df1 = cbind(df1,df2[,-1])
}
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,pcol=cols,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,c("Naive count","EM test","TDR test"),fill = cols,ncol = 2,border=F
text(0,0,"FDR",cex = 1.2,font = 2)

# plot the legends separately
dev.off()
par(mar = c(5,5,5,5))
plot(c(2,2))
legend(x="top",c("Naive count","MS test","TDR test"),fill = cols,ncol = 1,border=F,cex=2)

#######
# Look at the Pi1 results
#######
load("./simulations_strict/simulation_pi1_summ_stats.RData")

# Num discoveries
ndigits = 0
resultsdf = mean_num_discoveries
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","our 0.9","raw 0.9")]
resultsdf = mean_num_discoveries
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,c("Naive count","EM test","TDR test"),fill = cols,ncol = 2,border=F)
text(0,0,"N",cex = 1.5,font = 2)

# (D) IVW, MR-PRESSO: FDRs
ndigits = 2
resultsdf = mean_fprs
inds = resultsdf$deg==deg & resultsdf$p1 == p1 & resultsdf$p2==p2
df1 = resultsdf[inds,c("prob_pleio","our 0.1","raw 0.1")]
resultsdf = mean_num_discoveries
rownames(df1) = df1[,1]
df1 = df1[,-1]
df1 = t(df1)
df1 = as.data.frame(df1)
df1 = rbind(min(0.1,min(df1)),df1)
df1 = rbind(max(df1),df1)
axslabs = round(seq(min(df1),max(df1),length.out = 6),digits = ndigits)
cols = method2col[rownames(df1)[-c(1:2)]]
radarchart(df1,axistype=1,seg=5,plwd=2,caxislabels=axslabs,plty=5,vlcex = 1.5,calcex=1.2)
# legend(x=-1,y=1.8,c("Naive count","EM test","TDR test"),fill = cols,ncol = 2,border=F
text(0,0,"FDR",cex = 1.2,font = 2)



# library(reshape2);library(ggplot2)
# deg = 1.5
# p1 = 1e-03
# p2 = 0.001
# inds = mean_fdrs$deg==deg & mean_fdrs$p1 == p1 & mean_fdrs$p2==p2
# df1 = mean_fdrs[inds,c("prob_pleio","p1","p2","egger","c-egger","ivw","c-ivw","mrpresso","c-mrpresso")]
# df2 = sd_fdrs[inds,names(df1)]
# df1 = melt(df1,id.vars = c("prob_pleio","p1","p2"))
# df2 = melt(df2,id.vars = c("prob_pleio","p1","p2"))
# names(df1) = c("prob_pleio","p1","p2","Method","EmpiricalFDR")
# names(df2) = c("prob_pleio","p1","p2","Method","SD")
# df1$SD = df2$SD
# print(
#   ggplot(df1, aes(x=as.character(prob_pleio), y=EmpiricalFDR, fill=Method)) +
#     geom_bar(position=position_dodge(), stat="identity", colour='black') +
#     geom_errorbar(aes(ymin=EmpiricalFDR-SD, ymax=EmpiricalFDR+SD),na.rm=T, 
#                   width=.2,position=position_dodge(.9))
# )
# 
# inds = mean_fdrs$deg==deg & mean_fdrs$p1 == p1 & mean_fdrs$p2==p2
# df1 = mean_fdrs[inds,c("prob_pleio","p1","p2","edge_sep","edge_sep_test1","edge_sep_test2")]
# df2 = sd_fdrs[inds,names(df1)]
# df1 = melt(df1,id.vars = c("prob_pleio","p1","p2"))
# df2 = melt(df2,id.vars = c("prob_pleio","p1","p2"))
# names(df1) = c("prob_pleio","p1","p2","Method","EmpiricalFDR")
# names(df2) = c("prob_pleio","p1","p2","Method","SD")
# df1$SD = df2$SD
# print(
#   ggplot(df1, aes(x=as.character(prob_pleio), y=EmpiricalFDR, fill=Method)) +
#     geom_bar(position=position_dodge(), stat="identity", colour='black') +
#     geom_errorbar(aes(ymin=EmpiricalFDR-SD, ymax=EmpiricalFDR+SD),na.rm=T, 
#                   width=.2,position=position_dodge(.9))
# )
# 
# # Try boxplots
# df1 = melt(all_sim_results_fdrs,id.vars = c("prob_pleio","p1","p2","deg"))
# df2 = melt(all_sim_results_preds,id.vars = c("prob_pleio","p1","p2","deg"))
# 
# par(mfrow=c(3,2),mar=c(3,3,3,3))
# deg = 1.5
# p1 = 1e-02
# p2 = 0.01
# pleio = 0
# inds = df$deg==deg & df$p1 == p1 & df$p2==p2 & df$prob_pleio == pleio
# boxplot(value~variable,data=df1[inds,]);abline(h=0.1,col="red",lty=3,lwd=4)
# boxplot(value~variable,data=df2[inds,])
# pleio = 0.2
# inds = df$deg==deg & df$p1 == p1 & df$p2==p2 & df$prob_pleio == pleio
# boxplot(value~variable,data=df1[inds,]);abline(h=0.1,col="red",lty=3,lwd=4)
# boxplot(value~variable,data=df2[inds,])
# pleio = 0.4
# inds = df$deg==deg & df$p1 == p1 & df$p2==p2 & df$prob_pleio == pleio
# boxplot(value~variable,data=df1[inds,]);abline(h=0.1,col="red",lty=3,lwd=4)
# boxplot(value~variable,data=df2[inds,])





