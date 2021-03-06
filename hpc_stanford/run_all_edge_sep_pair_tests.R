##############################################################
# Helper functions for running the analysis
get_sh_prefix<-function(err="",log="",time="1:00:00"){
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

MAX_JOBS = 400
##############################################################
# Input data - files
genetic_ci_tests_plink_path = "~/cgauge_resub/genetic_CI_tests_results.RData"
skeleton_file = "~/cgauge_resub/Gs_skeleton.RData"
geno_data_path = "~/cgauge_resub/april2019_traits_genotypes/"
gwas_res_data = "~/cgauge_resub/april2019_causal_analysis_flow_input.RData"
gwas_res_path = "~/cgauge_resub/gwas_res/"
out_path = "~/cgauge_resub/ukbb_res/em_edge_sep_jobs/"
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
maf_file = paste(geno_data_path,"all_genotypes.frq",sep="")
mafs = read.table(maf_file,stringsAsFactors = F,header=T)
bim = read.table(paste(geno_data_path,"all_genotypes.bim",sep=""),stringsAsFactors = F)
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

# Get the overall pruned snp list
prlist = read.delim(paste(geno_data_path,"/plink.prune.in",sep=""),stringsAsFactors = F)[,1]
prlist = intersect(rownames(GWAS_Ps),prlist)
GWAS_Ps = GWAS_Ps[prlist,]

# define the skeleton using p1
G_t = skeleton_pmax < p1
diag(G_t) = F;mode(G_t)="numeric"

# Run all pairs
for(tr1 in colnames(GWAS_Ps)){
  for(tr2 in colnames(GWAS_Ps)){
    if(tr1==tr2){next}
    if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]==0){next}

    job_name = paste(tr1,"_vs_",tr2,sep="")
    rdata_name = paste(tr1,"_",tr2,"_input.RData",sep="")
    out_name = paste(tr1,"_",tr2,"_edgesep_em_output.RData",sep="")
    
    # Skip if results file exists
    if(out_name %in% list.files(out_path)){next}
    # add the full path to the file names
    rdata_name = paste(out_path,rdata_name,sep="")
    out_name = paste(out_path,out_name,sep="")
    
    # save the data of the current pair
    ps1 = GWAS_Ps[,tr2]
    ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][rownames(GWAS_Ps),"test3"]
    ps = cbind(ps1,ps_with_tr2_cond_tr1)
    save(ps,file=rdata_name)
    cmd = paste(
      "~/repos/cGAUGE/hpc_stanford/run_edge_sep_test_for_pair.R",
      "--file",rdata_name,
      "--testName grid",
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

# Read the results
edge_sep_em_res = c()
for(tr1 in colnames(GWAS_Ps)){
  for(tr2 in colnames(GWAS_Ps)){
    if(tr1==tr2){next}
    if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]==0){next}
    out_name = paste(tr1,"_",tr2,"_edgesep_em_output.RData",sep="")
    out_name = paste(out_path,out_name,sep="")
    currp = NA
    try({currp = get(load(out_name))})
    edge_sep_em_res = rbind(edge_sep_em_res,c(tr1,tr2,currp))
  }
}
edge_sep_em_res = data.frame(edge_sep_em_res,stringsAsFactors = F)
edge_sep_em_res[[3]] = as.numeric(as.character(edge_sep_em_res[[3]]))
save(edge_sep_em_res,file=paste(out_path,"edge_sep_em_res.RData",sep=""))
quantile(p.adjust(edge_sep_em_res[[3]]),na.rm=T)

##############################################################
# A few local tests
setwd("~/Desktop/causal_inference_projects/ms3/edge_sep_em/")
setwd("~/cgauge_resub/ukbb_res/em_edge_sep_jobs/")

load("./edge_sep_em_res.RData")

pval = get(load("Glycated_haemoglobin_HbA1c_HC221_edgesep_em_output.RData"))
ps = get(load("Glycated_haemoglobin_HbA1c_HC221_input.RData"))
p1 = ps[,1]
p2 = ps[,2]
univar_mixtools_em(p1,p2,reps=10)

tr1 = "Alanine_aminotransferase"
tr2 = "statins"
p1 = GWAS_Ps[,tr2]
p2 = trait_pair_pvals[[tr2]][[tr1]][rownames(GWAS_Ps),"test3"]
univar_mixtools_em(p1,p2,reps=3)
simple_lfdr_test(p1,p2)


