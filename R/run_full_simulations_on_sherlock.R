# This script runs the simulations using Stanford's Sherlock (a cluster)

###################################################################################
### Helper functions to run the simulations

get_sh_prefix<-function(err="",log="",time="1:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=", time,sep=""),
      "#SBATCH --partition=euan,mrivas,owners,normal",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --mem=2000",
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
  system(paste("sbatch -x sh-113-15 ",curr_sh_file,'&'))
}

###################################################################################
reps = 1
WD = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/simulations/"

tested_p1 = c(1e-02,1e-03,1e-04,1e-05)
tested_p2 = c(0.001,0.01,0.1)
tested_pleio_levels = c(0,0.1,0.2,0.3,0.4,0.5)
tested_degrees = c(0.5,1,1.5,2)

for(p1 in tested_p1){
  for(p2 in tested_p2){
    for(pleio_p in tested_pleio_levels){
      for(deg in tested_degrees){
        curr_folder = paste(WD,"deg",deg,"_pleio",pleio_p,"_p1",p1,"_p2",p2,"/",sep="")
        system(paste(mkdir,curr_folder))
        for(i in 1:reps){
          cmd = paste(
            "~/repos/cGAUGE/R/full_causal_graph_simulations.R",
            "--deg",deg,
            "--probPleio",pleio_p,
            "--p1",p1,
            "--p2",p2,
            "--out",paste(curr_folder,"sim_rep",i,".RData",sep="")
          )
          exec_cmd_on_sherlock(cmd,jobname = paste("sim_rep",i,sep=""),out_path = curr_folder)
        }
        
      }
    }
    break
  }
  break
}







