###################################################################################
# Set the session
required_libs = c("igraph","bnlearn","MRPRESSO",
                  "optparse","limma","MendelianRandomization",
                  "mixtools","locfdr")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
           })
}
# Add the cGAUGE functions and auxiliary functions for MR
# From GitHub
try({
  source("https://raw.githubusercontent.com/david-dd-amar/cGAUGE/master/R/cGAUGE.R")
  source("https://raw.githubusercontent.com/david-dd-amar/cGAUGE/master/R/twogroups_em_tests.R")
})
# From local clone (GitHub server sometimes has issues)
try({
  source("~/repos/cGAUGE/R/cGAUGE.R")
  source("~/repos/cGAUGE/R/twogroups_em_tests.R")
})
print("Completed loading libraries and code")
###################################################################################
###################################################################################
###################################################################################
# Set input data (see explanation of each file)

#' This is a table with phenotype code, description, and sample size
#' We also use this table to create a mapping from phenotype code to its name
#' (the object name is called icd2name but is not limited to icds only)
rivaslab_pheno_codes_file = "/oak/stanford/groups/mrivas/users/magu/repos/rivas-lab/wiki/ukbb/icdinfo/icdinfo.txt"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=F)
rownames(rivaslab_codes) = rivaslab_codes[,1]
icd2name = rivaslab_codes[,3];names(icd2name) = rownames(rivaslab_codes)

#### Analysis of 96 traits for the paper ####

#' This file contains an object with the p-values of all conditional independence tests
#'  for each variant G vs. a trait X. 
#' We represent this object using a named list of lists in which element [[tr1]][[tr2]]
#'is a matrix with the conditional independence results (p-values) for trait 1 conditioned on
#' trait 2 (rows are variants).
#' One of the columns in this matrix has the p-values, and these are used within cGAUGE's filters
genetic_ci_tests_plink_path = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_causal_analysis_flow_results/genetic_CI_tests_results.RData"
#' This file contains three objects:
#' marginal association p-values
#' maximal p-values for each pair (overl all tests)
#' a matrix with all potential separating sets (p>1e-10) for each pair
skeleton_file = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/Gs_skeleton.RData"
#' This is a prefix of the frq and bim files that contain the MAFs and SNP ids and locations
#' These are used to filter out low MAF SNPs and the MHC region
geno_data_path = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_genotypes/all_genotypes"
#' This RData file contains the standard (marginal) GWAS results,
#' adjusted for sex, age, and PCs
gwas_res_data = "/oak/stanford/groups/mrivas/users/davidama/april2019_causal_analysis_flow_input.RData"
#' Set a path for the output files
out_path = "/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/ukbb_res/"

###################################################################################
###################################################################################
###################################################################################

# Load the data and set output paths
system(paste("mkdir",out_path))
# load the plink analysis results
# this loads the list of pairwise CI tests
load(genetic_ci_tests_plink_path)
# load the skeleton: get a matrix of the maximal association p-values
# and an object of the separating sets for non-edges
load(skeleton_file)
skeleton_pmax = pmax_network
all_sepsets = sepsets

# this loads the standard GWAS results
load(gwas_res_data)

# set the phenotype names
pheno_names = icd2name[colnames(skeleton_pmax)]
pheno_names[colnames(skeleton_pmax)=="sex"] = "sex"
pheno_names[colnames(skeleton_pmax)=="age"] = "age"
names(pheno_names)[colnames(skeleton_pmax)=="sex"] = "sex"
names(pheno_names)[colnames(skeleton_pmax)=="age"] = "age"
# shorten some names
pheno_names[pheno_names=="Age_when_periods_started_(menarche)"] = "Menarche"
pheno_names[pheno_names=="Number_of_children_fathered"] = "Num children (f)"
pheno_names[pheno_names=="Red_blood_cell_(erythrocyte)_count"] = "RBC"
pheno_names[pheno_names=="Mean_corpuscular_haemoglobin"] = "Hemoglobin"
pheno_names[pheno_names=="Red_blood_cell_(erythrocyte)_distribution_width"] = "RBC_distr_width"
pheno_names[pheno_names=="White_blood_cell_(leukocyte)_count"] = "WBC"
pheno_names[pheno_names=="Platelet_distribution_width"] = "Platelet_distr_width"
pheno_names[pheno_names=="Mean_platelet_(thrombocyte)_volume"] = "Platelet_volume"
pheno_names[pheno_names=="Forced_vital_capacity_(FVC)"] = "FVC"
pheno_names[pheno_names=="Number_of_days/week_of_vigorous_physical_activity_10+_minutes"] = "Vigorous_PA"
pheno_names[pheno_names=="Number_of_days/week_of_moderate_physical_activity_10+_minutes"] = "Moderate_PA"
pheno_names[pheno_names=="Time_spend_outdoors_in_summer"] = "Summer_outdoor_time"
pheno_names[pheno_names=="Time_spent_outdoors_in_winter"] = "Winter_outdoor_time"
pheno_names[pheno_names=="Standing_height"] = "Height"
pheno_names[pheno_names=="Pulse_rate,_automated_reading"] = "Pulse rate"
pheno_names[pheno_names=="Average_weekly_beer_plus_cider_intake"] = "Beer intake"
pheno_names[pheno_names=="Systolic_blood_pressure,_automated_reading"] = "SBP"
pheno_names[pheno_names=="Diastolic_blood_pressure,_automated_reading" ] = "DBP"
pheno_names[pheno_names=="Fluid_intelligence_score" ] = "Intelligence"
pheno_names[pheno_names=="heart_attack/myocardial_infarction"] = "Myocardial infarction"
pheno_names[pheno_names=="heart_failure/pulmonary_odema"] = "Heart failure"
pheno_names[pheno_names=="large_bowel_cancer/colorectal_cancer" ] = "Colorectal cancer"
pheno_names[pheno_names=="inflammatory_bowel_disease"] = "IBD"
pheno_names[pheno_names=="gastro-oesophageal_reflux_(gord)_/_gastric_reflux" ] = "Gastric_reflux"
pheno_names[pheno_names=="hayfever/allergic_rhinitis" ] = "Allergic_rhinitis"
pheno_names[pheno_names=="urinary_tract_infection/kidney_infection" ] = "Kidney infection"
pheno_names[pheno_names=="hypothyroidism/myxoedema"] = "Hypothyroidism"
pheno_names = gsub(pheno_names,pattern="_",replacement = " ")
pheno_names = sapply(pheno_names,function(x)paste(toupper(substr(x,1,1)),substr(x,2,nchar(x)),sep=""))

# add missing names as is
missing_names = setdiff(colnames(snp_matrix),names(pheno_names))
pheno_names[missing_names] = missing_names
missing_names2 = setdiff(colnames(snp_matrix),names(icd2name))
icd2name[missing_names2] = missing_names2
save(pheno_names,icd2name,rivaslab_codes,file=
       paste(out_path,"all_pheno_name_metadata.RData",sep=""))

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

###################################################################################
###################################################################################
###################################################################################

P1s = c(1e-5,1e-6,1e-7,1e-8)
P2s = c(0.1,0.01,0.001)

#### Perform the G_t skeleton learning for the different p1 values #####
p12G_t = list()
for(p1 in P1s){
  sepsets = all_sepsets
  G_t = skeleton_pmax < p1
  diag(G_t) = F;mode(G_t)="numeric"
  # Print G_T to a text file
  G_t_edges = c()
  for(i in 2:nrow(G_t)){
    for(j in 1:(i-1)){
      if(!is.na(G_t[i,j]) && G_t[i,j]>0){
        G_t_edges = rbind(G_t_edges,c(pheno_names[colnames(G_t)[i]],pheno_names[colnames(G_t)[j]]))
      }
    }
  }
  write.table(G_t_edges,file=paste(out_path,"G_t_edges","p1_",p1,".txt",sep=""),sep="\t",
              row.names = F,col.names = F,quote = F)
  
  # Clean the separating sets
  for(nn1 in names(sepsets)){
    for(nn2 in names(sepsets[[nn1]])){
      m = sepsets[[nn1]][[nn2]]
      if(is.null(dim(m))|| nrow(m)<2){next}
      m_pvals = as.numeric(m[,2])
      to_rem = m_pvals <= p1
      m = m[!to_rem,]
      # print(nrow(m))
      sepsets[[nn1]][[nn2]] = m
    }
  }
  # transform to lists and remove the non minimal separating sets
  # (this may take some time)
  p1_sepsets = list()
  for(nn1 in names(sepsets)){
    print(nn1)
    for(nn2 in names(sepsets[[nn1]])){
      m = unique(sepsets[[nn1]][[nn2]])
      if(is.null(dim(m))|| nrow(m)<2){next}
      l = lapply(m[,1], function(x)strsplit(x,split=",")[[1]])
      l1 = remove_non_minimal_sepsets(l)
      p1_sepsets[[nn1]][[nn2]] = l1
    }
  }
  
  p12G_t[[as.character((p1))]] = list(
    G_t=G_t,mseps = p1_sepsets
  )
}

save(p12G_t,file = paste(out_path,"p12G_t.RData",sep=""))

#### Run the EdgeSep tests for the maximum p1 (contains the results for lower p1 values) #####
# The code below is commented out because it is rather slow for real data
# A parallel version of this code was run on Stanford's HPC, see the scripts in hpc_stanford:
# run_edge_sep_test_for_pair.R - wrapper for a single job per pair
# run_all_edge_sep_pair_tests.R - a script that loads the data and runs the EM test for all pairs
# Below the commented code we load the results of these scripts and use them subsequently.
# We leave this code for the readers interetsed in replication.
# p1 = max(P1s)
# G_t = p12G_t[[as.character((p1))]]$G_t
# 
# # Perform the EdgeSep statistical tests
# NonNA_GWAS_Ps = GWAS_Ps
# NonNA_GWAS_Ps[is.na(NonNA_GWAS_Ps)] = 0.5
# for(n1 in names(trait_pair_pvals)){
#   for(n2 in names(trait_pair_pvals[[n1]])){
#     m = trait_pair_pvals[[n1]][[n2]]
#     m = m[rownames(GWAS_Ps),]
#     m[is.na(m)] = 0.5
#     trait_pair_pvals[[n1]][[n2]] = m
#   }
#   gc()
# }
# 
# edge_sep_results_statTest2 = EdgeSepTest(NonNA_GWAS_Ps,G_t,trait_pair_pvals,
#                                          text_col_name="test3",test = simple_lfdr_test)
# edge_sep_results_statTest1 = EdgeSepTest(NonNA_GWAS_Ps,G_t,trait_pair_pvals,
#                                          text_col_name="test3",test = univar_mixtools_em)

# Load the EdgeSep results, create output files and networks
load(paste(out_path,"p12G_t.RData",sep=""))
load("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/ukbb_res/em_edge_sep_jobs/edge_sep_em_res.RData")

# For MR: check different combinations of the input parameters p1 and p2
for(p1 in P1s){
  G_t = p12G_t[[as.character((p1))]]$G_t
  mseps = p12G_t[[as.character((p1))]]$mseps
  for (p2 in P2s){
    G_vt = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,P1=p1,
                                 P2=p2,test_columns = c("test2","test3"))[[1]]
    uniquely_mapped_ivs = rownames(G_vt)[rowSums(G_vt)==1]
    
    # Get the iv sets for each MR analysis
    iv_sets_thm21 = list()
    iv_sets_thm22 = list()
    for(tr1 in colnames(GWAS_Ps)){
      iv_sets_thm21[[tr1]] = list()
      iv_sets_thm22[[tr1]] = list()
      for(tr2 in colnames(GWAS_Ps)){
        if(tr1==tr2){next}
        iv_sets_thm21[[tr1]][[tr2]] = rownames(GWAS_Ps)[!is.na(GWAS_Ps[,tr1]) & GWAS_Ps[,tr1]<p1]
        currseps = unique(unlist(mseps[[tr1]][[tr2]]))
        currseps = setdiff(currseps,c("sex","age"))
        # remove IVs into separating variables
        for(sep in currseps){
          curr_sep_ivs = rownames(G_vt)[G_vt[,sep]]
          iv_sets_thm21[[tr1]][[tr2]] = setdiff(iv_sets_thm21[[tr1]][[tr2]],curr_sep_ivs)
        }
        
        iv_sets_thm22[[tr1]][[tr2]] = rownames(GWAS_Ps)[!is.na(GWAS_Ps[,tr1]) & GWAS_Ps[,tr1]<p1]
        iv_sets_thm22[[tr1]][[tr2]] = intersect(iv_sets_thm22[[tr1]][[tr2]],uniquely_mapped_ivs)
      }
    }
    
    # Analysis 5.1: simple meta-analysis on all pairs
    meta_anal_res_thm21 = run_pairwise_pval_combination_analysis_from_iv_sets(iv_sets_thm21,GWAS_Ps)
    meta_anal_res_thm22 = run_pairwise_pval_combination_analysis_from_iv_sets(iv_sets_thm22,GWAS_Ps)
    ivw_res_thm21 = run_pairwise_mr_analyses_with_iv_sets(
      sum_stat_matrix,sum_stat_se_matrix,iv_sets_thm21,func=mr_ivw,robust=T)
    ivw_res_thm22 = run_pairwise_mr_analyses_with_iv_sets(
      sum_stat_matrix,sum_stat_se_matrix,iv_sets_thm22,func=mr_ivw,robust=T)
    print("Done updating the MR results")
    thm21_res = combine_mm_mr_analyses(meta_anal_res_thm21,ivw_res_thm21,
                                        p_h_thr = -1,minIVs = 3)
    thm22_res = combine_mm_mr_analyses(meta_anal_res_thm22,ivw_res_thm22,
                                       p_h_thr = -1,minIVs = 3)
    thm21_res = add_is_non_edge_column(thm21_res,G_t)
    thm22_res = add_is_non_edge_column(thm22_res,G_t)
    thm21_res = add_edgesep_res_column(thm21_res,edge_sep_em_res,G_t)
    thm22_res = add_edgesep_res_column(thm22_res,edge_sep_em_res,G_t)
    
    # Represent the selected edges nicely
    thm21_res[,1] = pheno_names[thm21_res[,1]]
    thm21_res[,2] = pheno_names[thm21_res[,2]]
    thm22_res[,1] = pheno_names[thm22_res[,1]]
    thm22_res[,2] = pheno_names[thm22_res[,2]]
    
    write.table(thm21_res,file=paste(out_path,
                "mr_thm21_res_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    thm21_res2 = thm21_res[as.numeric(thm21_res[,"numIVs"])>9,]
    write.table(thm21_res2,file=paste(out_path,
                "mr_thm21_res_atleast_10_ivs_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    
    write.table(thm22_res,file=paste(out_path,
                "mr_thm22_res_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    thm22_res2 = thm22_res[as.numeric(thm22_res[,"numIVs"])>9,]
    write.table(thm22_res2,file=paste(out_path,
                  "mr_thm22_res_atleast_10_ivs_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    
    write.table(thm22_res[grepl("cancer",thm22_res[,2]),c(1:4,9)],quote=F,sep="\t")
    write.table(thm22_res[grepl("LDL",thm22_res[,1]),c(1:4,9:10)],quote=F,sep="\t")
    write.table(thm22_res[grepl("HDL",thm22_res[,1]),c(1:4,9:10)],quote=F,sep="\t")
    
    ####################################################################################################
    # save the results of the analysis for further examination
    save(
      GWAS_Ps, # The original associations of the GWAS without conditional independence filtering
      meta_anal_res_thm21,meta_anal_res_thm22,# meta-analysis results
      G_t,G_vt, # the skeletons
      iv_sets_thm21,iv_sets_thm22, # selected instruments per pair
      ivw_res_thm21,ivw_res_thm22, # raw MR results 
      thm21_res,thm22_res, # MR+Meta-analysis results after filters
      file = paste(out_path,"cgauge_res_",p1,"_",p2,".RData",sep="")
    )
  }
}

################################################################
################################################################
################################################################
################################################################
# Add MR-PRESSO estimates
################################################################
################################################################
################################################################
################################################################
presso_runs_path = paste(out_path,"mrpresso_runs/",sep="")
system(paste("mkdir",presso_runs_path))
# Subset of the data - save some time
P1s = c(1e-05,1e-6,1e-7,1e-8)
P2s = c(0.01,0.001)
# Load the EdgeSep results, create output files and networks
load(paste(out_path,"p12G_t.RData",sep=""))
load("/oak/stanford/groups/mrivas/users/davidama/cgauge_resub/ukbb_res/em_edge_sep_jobs/edge_sep_em_res.RData")

#####
# Helper functions for running in Stanford's HPC
get_sh_prefix<-function(err="",log="",time="3:00:00"){
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
#######

# First, run all the jobs
for(p1 in P1s){
  print(paste("p1",p1))
  for (p2 in P2s){
    print(paste("p2",p2))
    load(paste(out_path,"cgauge_res_",p1,"_",p2,".RData",sep=""))
    if("mrpresso_thm22_res" %in% ls()){next}
    for(tr1 in names(iv_sets_thm22)){
      for(tr2 in names(iv_sets_thm22)){
        if(tr1==tr2){next}
        ivs = iv_sets_thm22[[tr1]][[tr2]]
        if(length(ivs)<5){next}
        curr_job_name = paste(tr1,"_",tr2,"_",p1,"_",p2,sep="")
        curr_out_file = paste(presso_runs_path,curr_job_name,".RData",sep="")
        curr_ivs_file = paste(presso_runs_path,curr_job_name,"_ivs.RData",sep="")
        if(file.exists(curr_out_file)){next}
        save(ivs,file=curr_ivs_file)
        cmd = paste(
          "~/repos/cGAUGE/hpc_stanford/run_mrpresso_on_pair.R",
          "--tr1",tr1,
          "--tr2",tr2,
          "--ivs_rdata",curr_ivs_file,
          "--out",curr_out_file
        )
        print(curr_job_name)
        exec_cmd_on_sherlock(cmd,jobname = curr_job_name,out_path = presso_runs_path)
        
        job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
        num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
        while(num_active_jobs > 500){
          Sys.sleep(5)
          job_state = system2("squeue",args = list("-u davidama | wc"),stdout=TRUE)
          num_active_jobs = as.numeric(strsplit(job_state,split="\\s+",perl = T)[[1]][2])
        }
      }
    }
  }
}

# Next, read and parse the results, save the new cgauge output objects
load(paste(out_path,"all_pheno_name_metadata.RData",sep=""))
for(p1 in P1s){
  for (p2 in P2s){
    load(paste(out_path,"cgauge_res_",p1,"_",p2,".RData",sep=""))
    # Get the results in a similar format to the other MR analyses
    if("mrpresso_thm22_res" %in% ls()){next}
    cgauge_mrpresso_thm22 = c()
    for(tr1 in names(iv_sets_thm22)){
      print(tr1)
      for(tr2 in names(iv_sets_thm22)){
        if(tr1==tr2){next}
        ivs = iv_sets_thm22[[tr1]][[tr2]]
        if(length(ivs)<5){next}
        curr_job_name = paste(tr1,"_",tr2,"_",p1,"_",p2,sep="")
        curr_out_file = paste(presso_runs_path,curr_job_name,".RData",sep="")
        if(!file.exists(curr_out_file)){next}
        mrpresso_res = NULL
        load(curr_out_file)
        if(is.null(mrpresso_res)){
          print(paste("NULL mrpresso results",curr_job_name))
          next
        }
        causalp = mrpresso_res[[1]][2,"P-value"]
        causalest = mrpresso_res[[1]][2,"Causal Estimate"]
        if(is.na(causalp)){
          causalp = mrpresso_res[[1]][1,"P-value"]
          causalest = mrpresso_res[[1]][1,"Causal Estimate"]
        }
        p_het = mrpresso_res[["MR-PRESSO results"]][["Global Test"]]$Pvalue
        v = c(mrpresso_res$tr1,mrpresso_res$tr2,causalp,p_het,causalest,NA,
            length(iv_sets_thm22[[mrpresso_res$tr1]][[mrpresso_res$tr2]]))
        names(v) = c("Exposure","Outcome","p","p_het","est","Q","NumIVs")
        cgauge_mrpresso_thm22 = rbind(cgauge_mrpresso_thm22,v)
      }
    }
    
    cgauge_mrpresso_thm22 = as.data.frame(cgauge_mrpresso_thm22)
    for(j in 3:ncol(cgauge_mrpresso_thm22)){
      cgauge_mrpresso_thm22[[j]] = as.numeric(as.character(cgauge_mrpresso_thm22[[j]]))
    }
    # combine with the other analyses and save the results
    meta_anal_res_thm21 = run_pairwise_pval_combination_analysis_from_iv_sets(iv_sets_thm21,GWAS_Ps)
    meta_anal_res_thm22 = run_pairwise_pval_combination_analysis_from_iv_sets(iv_sets_thm22,GWAS_Ps)
    presso_thm22_res = combine_mm_mr_analyses(meta_anal_res_thm22,cgauge_mrpresso_thm22,
                                       p_h_thr = -1,minIVs = 3)
    presso_thm22_res = add_is_non_edge_column(presso_thm22_res,G_t)
    presso_thm22_res = add_edgesep_res_column(presso_thm22_res,edge_sep_em_res,G_t)
    # Represent the selected edges nicely
    presso_thm22_res[,1] = pheno_names[presso_thm22_res[,1]]
    presso_thm22_res[,2] = pheno_names[presso_thm22_res[,2]]
    
    write.table(presso_thm22_res,file=paste(out_path,
                                     "mr_thm22_res_mrpresso_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    presso_thm22_res2 = presso_thm22_res[as.numeric(presso_thm22_res[,"numIVs"])>9,]
    write.table(presso_thm22_res2,file=paste(out_path,
                                      "mr_thm22_res_mrpresso_atleast_10_ivs_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    
    write.table(presso_thm22_res[grepl("cancer",presso_thm22_res[,2]),c(1:4,9)],quote=F,sep="\t")
    write.table(presso_thm22_res[grepl("LDL",presso_thm22_res[,1]),c(1:4,9:10)],quote=F,sep="\t")
    write.table(presso_thm22_res[grepl("HDL",presso_thm22_res[,1]),c(1:4,9:10)],quote=F,sep="\t")
    
    ####################################################################################################
    # save the results of the analysis for further examination
    mrpresso_thm22_res_raw = cgauge_mrpresso_thm22
    mrpresso_thm22_res = presso_thm22_res
    save(
      GWAS_Ps, # The original associations of the GWAS without conditional independence filtering
      G_t,G_vt, # the skeletons
      meta_anal_res_thm21,meta_anal_res_thm22,
      iv_sets_thm21,iv_sets_thm22, # selected instruments per pair
      ivw_res_thm21,ivw_res_thm22, # raw MR results 
      thm21_res,thm22_res, # MR+Meta-analysis results after filters
      mrpresso_thm22_res_raw, # MRPRESSO raw results
      mrpresso_thm22_res, # MRPRESSO filtered results
      file = paste(out_path,"cgauge_res_",p1,"_",p2,".RData",sep="")
    )
    
    # remove iv, err and log files
    setwd(presso_runs_path)
    before_del = length(list.files("."))
    system(paste("rm ","*_",p1,"_",p2,"*.log",sep=""))
    system(paste("rm ","*_",p1,"_",p2,"*.err",sep=""))
    system(paste("rm ","*_",p1,"_",p2,"*_ivs.RData",sep=""))
    after_del = length(list.files("."))
    print(paste("deleted log, err, ivs files, number of deleted:",before_del-after_del))
    
    # rm objects
    rm(mrpresso_thm22_res)
    rm(mrpresso_thm22_res_raw)
  }
}


if(!is.null(dim(cgauge_mrpresso_res))){
  cgauge_mrpresso_res = as.data.frame(cgauge_mrpresso_res)
  for(j in 3:ncol(cgauge_mrpresso_res)){
    cgauge_mrpresso_res[[j]] = as.numeric(as.character(cgauge_mrpresso_res[[j]]))
  }
}
  
