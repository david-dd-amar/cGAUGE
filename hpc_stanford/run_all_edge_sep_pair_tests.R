##############################################################
# Helper functions for running the analysis
get_sh_prefix<-function(err="",log="",time="5:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=", time,sep=""),
      "#SBATCH --partition=euan,mrivas,owners,normal",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --mem=4000",
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

MAX_JOBS = 300
##############################################################
# Input data - files
genetic_ci_tests_plink_path = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_causal_analysis_flow_results/genetic_CI_tests_results.RData"
skeleton_file = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/Gs_skeleton.RData"
geno_data_path = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_genotypes/all_genotypes"
gwas_res_data = "/oak/stanford/groups/mrivas/users/davidama/april2019_causal_analysis_flow_input.RData"
gwas_res_path = "/oak/stanford/groups/mrivas/users/davidama/gwas_res/"
out_path = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/ukbb_res/em_edge_sep_jobs"
system(paste("mkdir",out_path))
# Input data- load into session
load(genetic_ci_tests_plink_path)
load(gwas_res_data)
load(skeleton_file)
skeleton_pmax = pmax_network
p1 = 1e-05

# Define the GWAS results
# Select SNPs by their MAF
pruned_snp_list = cl_unified_list
pruned_snp_lists = code2clumped_list
MAF = 0.01
maf_file = paste(geno_data_path,".frq",sep="")
mafs = read.table(maf_file,stringsAsFactors = F,header=T)
bim = read.table(paste(geno_data_path,".bim",sep=""),stringsAsFactors = F)
mhc_snps = bim[bim[,1]==6 & bim[,4]>23000000 & bim[,4]<35000000,2]
our_snps = mafs$SNP[mafs$MAF >= MAF]
pruned_snp_list = intersect(pruned_snp_list,our_snps)
pruned_snp_list = intersect(pruned_snp_list,rownames(snp_matrix))
pruned_snp_list = setdiff(pruned_snp_list,mhc_snps)
for(nn in names(pruned_snp_lists)){
  pruned_snp_lists[[nn]] = intersect(pruned_snp_lists[[nn]],pruned_snp_list)
}
iv2trait_p = snp_matrix[pruned_snp_list,]
maf_as_weights = mafs$MAF;names(maf_as_weights) = mafs$SNP
GWAS_Ps = iv2trait_p[pruned_snp_list,]

# define the skeleton using p1
G_t = skeleton_pmax < p1
diag(G_t) = F;mode(G_t)="numeric"

# Avoid NAs, make the data over the same snp names
NonNA_GWAS_Ps = GWAS_Ps
NonNA_GWAS_Ps[is.na(NonNA_GWAS_Ps)] = 0.5
for(n1 in names(trait_pair_pvals)){
  for(n2 in names(trait_pair_pvals[[n1]])){
    m = trait_pair_pvals[[n1]][[n2]]
    m = m[rownames(GWAS_Ps),]
    m[is.na(m)] = 0.5
    trait_pair_pvals[[n1]][[n2]] = m
  }
  gc()
}

# Run all pairs
for(tr1 in colnames(GWAS_Ps)){
  for(tr2 in colnames(GWAS_Ps)){
    if(tr1==tr2){next}
    if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]==0){next}
    ps1 = GWAS_Ps[,tr2]
    if(is.null(rownames(GWAS_Ps))){
      ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][,text_col_name]
    }
    else{
      ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][rownames(GWAS_Ps),text_col_name]
    }
    
    ps = cbind(ps1,ps_with_tr2_cond_tr1)
    rdata_name = paste(tr1,"_",tr2,"_input.RData",sep="")
    out_name = paste(tr1,"_",tr2,"_edgesep_em_output.RData",sep="")
    job_name = paste(tr1,"_vs_",tr2,sep="")
    save(ps,file=rdata_name)
    cmd = paste(
      "~/repos/cGAUGE/hpc_stanford/run_edge_sep_test_for_pair.R",
      "--file",rdata_name,
      "--out",out_name
    )
    exec_cmd_on_sherlock(cmd,jobname = job_name,out_path = out_path)
    
    job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
    num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
    while(num_active_jobs > MAX_JOBS){
      Sys.sleep(5)
      job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
      num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
    }
  }
}


