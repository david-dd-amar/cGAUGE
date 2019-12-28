###################################################################################
# Set the session
required_libs = c("igraph","bnlearn","MRPRESSO",
                  "optparse","limma","MendelianRandomization",
                  "mixtools")
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

geno_data_path = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_genotypes/all_genotypes"
gwas_res_data = "/oak/stanford/groups/mrivas/users/davidama/april2019_causal_analysis_flow_input.RData"
gwas_res_path = "/oak/stanford/groups/mrivas/users/davidama/gwas_res/"

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
# this loads the standard GWAS results
load(gwas_res_data)
# load the skeleton: get a matrix of the maximal association p-values
# and an object of the separating sets for non-edges
load(skeleton_file)
skeleton_pmax = pmax_network
all_sepsets = sepsets

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

# Check different combinations of the input parameters p1 and p2
P1s = c(1e-5,1e-6,1e-7,1e-8)
P2s = c(0.1,0.01,0.001)
for(p1 in P1s){
  
  #######################################################################################
  # Infer the G_T skeleton (use p1 as the threshold for significance)
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
  skeleton_pthr[is.na(skeleton_pthr)]=0
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
    for(nn2 in names(sepsets[[nn1]])){
      m = unique(sepsets[[nn1]][[nn2]])
      if(is.null(dim(m))|| nrow(m)<2){next}
      l = lapply(m[,1], function(x)strsplit(x,split=",")[[1]])
      l = remove_non_minimal_sepsets(l)
      p1_sepsets[[nn1]][[nn2]] = l
    }
  }
  
  for (p2 in P2s){
    
    G_vt = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,P1=p1,
                                 P2=p2,test_columns = c("test2","test3"))[[1]]
    
    # Perform the EdgeSep statistical tests
    edge_sep_results_statTest1 = EdgeSepTest(GWAS_Ps,G_t,trait_pair_pvals,text_col_name="test3",
                                             test = univar_mixtools_em)
    edge_sep_results_statTest2 = EdgeSepTest(GWAS_Ps,G_t,trait_pair_pvals,text_col_name="test3",
                                             test = simple_lfdr_test)
        
    # Print the resulting scored network to files
    edge_orientation_res = c()
    for(nn in names(detected_cis_per_edge)){
      arr = strsplit(nn,split=";")[[1]][c(4,6)]
      score1 = length(detected_cis_per_edge[[nn]]$variants)
      score2 = score1/detected_cis_per_edge[[nn]][[1]]
      score3 = score1/detected_cis_per_edge[[nn]][[2]]
      new_edge = c(arr,score1,score2,score3)
      print(new_edge)
      edge_orientation_res = rbind(edge_orientation_res,new_edge)
    }
    rownames(edge_orientation_res) = NULL
    colnames(edge_orientation_res) = c("X->","Y","NumSepIVs","Percentage_vs_tr1","Percentage_vs_tr1_tr2")
    write.table(edge_orientation_res,file=paste(out_path,"edge_orientation_res_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    edge_orientation_res2 = edge_orientation_res[as.numeric(edge_orientation_res[,3])>4,]
    write.table(edge_orientation_res2,file=paste(out_path,"edge_orientation_res__atleast_5_ivs_",
                                                 p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    
    
    # Run MR analysis
    # Get the iv sets for each MR analysis
    # Analysis 5.1: simple meta-analysis on all pairs
    meta_anal_res = run_pairwise_pval_combination_analyses(G_it,GWAS_Ps,
                        pruned_lists=code2clumped_list,weights=maf_as_weights,maxp=0.001)
    iv2_res = run_pairwise_mr_analyses_with_iv_sets(sum_stat_matrix,sum_stat_se_matrix,iv_sets,
                                                    func=mr_ivw,robust=T)
    print("Done updating the MR results")
    
    cleaned_Egger_res = combine_mm_mr_analyses(meta_anal_res,mr_anal_res[["Egger"]],
                                               pi1_thr=2,p_h_thr = -1,minIVs = 3)
    colnames(cleaned_Egger_res) = c("tr1->","tr2","p_MR","Est","pi1","numIVs")
    cleaned_Egger_res_non_pleio = clean_non_pleio_pairs(cleaned_Egger_res,G_t)
    is_non_pleio = is.element(rownames(cleaned_Egger_res),set=rownames(cleaned_Egger_res_non_pleio))
    cleaned_Egger_res = cbind(cleaned_Egger_res,is_non_pleio)
    effect_direction = rep("Up",nrow(cleaned_Egger_res))
    effect_direction[as.numeric(cleaned_Egger_res[,"Est"])<0] = "Down"
    cleaned_Egger_res = cbind(cleaned_Egger_res,effect_direction)
    
    # Represent the selected edges nicely
    mr_edges = cleaned_Egger_res
    rownames(mr_edges) = NULL
    mr_edges[,1] = pheno_names[mr_edges[,1]]
    mr_edges[,2] = pheno_names[mr_edges[,2]]
    mr_edges[,3] = -log(as.numeric(as.character(mr_edges[,3]))+1e-100,10)
    write.table(mr_edges,file=paste(out_path,"mr_Egger_non_pleio_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    mr_edges2 = mr_edges[as.numeric(mr_edges[,6])>9,]
    write.table(mr_edges2,file=paste(out_path,"mr_Egger_non_pleio_atleast_10_ivs_",p1,"_",p2,".txt",sep=""),
                quote=F,row.names = F,col.names = T,sep="\t")
    
    ####################################################################################################
    # save the results of the analysis for further examination
    save(
      G_it_0, # The original associations of the GWAS without conditional independence filtering
      G_t,G_it, # the skeleton of the traits
      G_t_edges, # nice representation of skeleton edges
      detected_cis_per_edge, # Edges that cause separation
      edge_orientation_res, # Summary of the disappearing assoc analysis
      updated_skel, # loci-trait skeleton
      newly_formed_sigs, # Emerging associations
      partial_ev_res, # summary of the emerging assoc analysis
      cleaned_Egger_res, # MR+Meta-analysis results after filters
      mr_edges, # MR+MM discovered edges - same as cleaned_Egger_res but formatted nicely
      mr_anal_res, # All MR results
      meta_anal_res, # All meta-analysis results with Stouffers method and maxp=0.001
      file = paste(out_path,"three_rule_analysis_",p1,"_",p2,".RData",sep="")
    )
  }
}



