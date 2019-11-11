# This script runs the simulations using Stanford's Sherlock (a cluster)

###################################################################################
### Helper functions to run the simulations

get_sh_prefix<-function(err="",log="",time="0:20:00"){
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
reps = 50
WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations_v2/"
MAX_JOBS = 500

tested_p1 = c(1e-02,1e-03,1e-04,1e-05)
tested_p2 = c(0.001,0.01,0.1)
tested_pleio_levels = c(0,0.1,0.2,0.3,0.4,0.5)
tested_degrees = c(1,1.5,2,2.5,3)

for(p1 in tested_p1){
  print(paste("p1",p1))
  for(p2 in tested_p2){
    print(paste("p2",p2))
    if(p2<p1){next}
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
FDR = 0.1

all_sim_results_errs = c()
all_sim_results_preds = c()
for(p1 in tested_p1){
  print(paste("p1",p1))
  for(p2 in tested_p2){
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
          egger_qs = p.adjust(standard_mr_results$Egger$p)
          errs["egger"] = sum(egger_qs < FDR & 
                   standard_mr_results$Egger$KnownDistance==-1,na.rm = T)
          preds["egger"] = sum(egger_qs < FDR,na.rm = T)
          ivw_qs = p.adjust(standard_mr_results$IVW$p)
          errs["ivw"] = sum(ivw_qs < FDR & standard_mr_results$IVW$KnownDistance==-1,na.rm = T)
          preds["ivw"] = sum(ivw_qs < FDR,na.rm = T)
          presso_q = p.adjust(standard_mr_results$MRPRESSO$`P-value`)
          errs["mrpresso"] = sum(presso_q < FDR & standard_mr_results$MRPRESSO$KnownDistance==-1,na.rm = T)
          preds["mrpresso"] = sum(presso_q < FDR,na.rm = T)
          lcv_q = p.adjust(standard_mr_results$LCV$P)
          errs["lcv"] = sum(lcv_q < FDR & standard_mr_results$LCV$KnownDistance==-1,na.rm = T)
          preds["lcv"] = sum(lcv_q < FDR,na.rm = T)
          
          egger_qs = p.adjust(cgauge_mr_results$Egger$p)
          errs["c-egger"] = sum(egger_qs < FDR & cgauge_mr_results$Egger$KnownDistance==-1
                                & cgauge_mr_results$PleioProperty,na.rm = T)
          preds["c-egger"] = sum(egger_qs < FDR,na.rm = T)
          ivw_qs = p.adjust(cgauge_mr_results$IVW$p)
          errs["c-ivw"] = sum(ivw_qs < FDR & cgauge_mr_results$IVW$KnownDistance==-1,na.rm = T)
          preds["c-ivw"] = sum(ivw_qs < FDR,na.rm = T)
          presso_q = p.adjust(cgauge_mr_results$MRPRESSO$`P-value`)
          errs["c-mrpresso"] = sum(presso_q < FDR & cgauge_mr_results$MRPRESSO$KnownDistance==-1,na.rm = T)
          preds["c-mrpresso"] = sum(presso_q < FDR,na.rm = T)
          method2false_discoveries = rbind(method2false_discoveries,errs)
          method2num_discoveries = rbind(method2num_discoveries,preds)
          
          
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
aggregate(all_sim_results_errs[,1:7],
          by=list(p1=all_sim_results_errs$p1,p2=all_sim_results_errs$p2,
                  deg = all_sim_results_errs$deg,prob_pleio = all_sim_results_errs$prob_pleio),
          FUN=mean)




