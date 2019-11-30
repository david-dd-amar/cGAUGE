# This script runs the simulations using Stanford's Sherlock (a cluster)

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
reps = 20
WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_2/"
MAX_JOBS = 300

tested_p1 = c(1e-02,1e-03,1e-04,1e-05)
tested_p2_factors = c(1,10,100)
tested_pleio_levels = c(0,0.1,0.2,0.3,0.4)
tested_degrees = c(1,1.25,1.5,1.75,2)

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
        # job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
        job_state = system2("sacct",args = list("| grep PEND | wc"),stdout=TRUE)
        num_curr_waiting_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
        while(num_curr_waiting_jobs > MAX_JOBS){
          Sys.sleep(5)
          job_state = system2("sacct",args = list("| grep PEND | wc"),stdout=TRUE)
          num_curr_waiting_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
        }
        
        for(i in 1:reps){
          
          curr_out_file = paste(curr_folder,"sim_rep",i,".RData",sep="")
          if(file.exists(curr_out_file)){next}
          cmd = paste(
            "~/repos/cGAUGE/R/full_causal_graph_simulations.R",
            "--deg",deg,
            "--probPleio",pleio_p,
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

##############################################################################################
# Go over the results

# Check the errors
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
        curr_files = curr_files[grepl("err$",curr_files)]
        for(f in curr_files){
          l = readLines(paste(curr_folder,f,sep=""))
          err_lines = l[grepl("error",l,ignore.case = T)]
          if(length(err_lines)>0){
            print(l)
          }
        }
      }
    }
  }
}

# Read the simulation results
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
          
          # EdgeSepTest
          edge_sep_results_statTest = 
            edge_sep_results_statTest[
              p.adjust(edge_sep_results_statTest$`pval:trait1->trait2`,method=FDR_method)<FDR,]
          errs["edge_sep_test"] = sum(!is_causal(edge_sep_results_statTest$KnownDistance))
          preds["edge_sep_test"] = nrow(edge_sep_results_statTest)
          
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
mean_errs = aggregate(all_sim_results_errs[,1:9],
          by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
          deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
          FUN=mean)

sd_errs = aggregate(all_sim_results_errs[,1:9],
                      by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                              deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                      FUN=sd)


all_sim_results_fdrs = (all_sim_results_errs/all_sim_results_preds)[,1:9]
all_sim_results_fdrs[is.nan(all_sim_results_fdrs)] = 0
all_sim_results_fdrs = cbind(all_sim_results_fdrs,all_sim_results_errs[,c("p1","p2","deg","prob_pleio")])
mean_fdrs = aggregate(all_sim_results_fdrs,
                      by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                              deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                      FUN=median,na.rm=T)

sd_fdrs = aggregate(all_sim_results_fdrs,
                    by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                            deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
                    FUN=sd,na.rm=T)
