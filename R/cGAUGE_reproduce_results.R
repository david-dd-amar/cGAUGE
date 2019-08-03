
# We assume that summary statistic data were downloaded from the following links
# and put into the working dir.
# The usage of each of these files is given in the code below. 
#
#
# 1.The results of the single GWAS (one per trait) - single_gwas_res.RData
#       https://drive.google.com/file/d/1XNZSYlDnepnPdLgG5qBrtTHrlo2Yq7IG/view?usp=sharing
#       This single GWAS result objects contain the data of variants that had p<0.001 in at least one 
#       of the traits analyzed in the study.
#       The useful objects in this RData are: 
#       snp_matrix - GWAS p-values (rows are variants, columns are traits)
#       sum_stat_matrix - GWAS effect sizes (rows are variants, columns are traits)
#       sum_stat_se_matrix - GWAS effect size standard error (rows are variants, columns are traits)
#       cl_unified_list - a list of clumped variant ids per trait (rows are variants, columns are traits)
#       code2clumped_list - lists of clumped variant sets per trait
# 2.icdinfo.txt - a (large) table with information about the traits covered in the Rivaslab GBE. We use this file
#   here mainly for names (in the loaded objects from the RData files the trait names are the codes in the GBE, which
#   is not always human readable)
#     link: https://drive.google.com/file/d/1TveaMn38xAu-r7KKq4v2NtMl5u5LwWuI/view?usp=sharing
# 3.genotypes bim and mafs to filter variants by MHC and/or by minor allele frequency, links:
#       frq - https://drive.google.com/file/d/1uieq63XKxuCAKGRxaknm1bVWNvUiGLAY/view?usp=sharing
#       bim - https://drive.google.com/file/d/1I9WqATOQ2SjEQDNSlvXvTi9wbXyZyTZl/view?usp=sharing      
# 4. Gs_skeleton.RData - the skeleton information from the trait-trait association analysis
#       https://drive.google.com/file/d/1CGav4eGQLi-G1zCdqyrSXbGL8b_aseGM/view?usp=sharing
# 5. genetic_CI_tests_results.RData - this file is large (11Gb) and contains the list of ALL conditional
#   independence tests among 96 traits form the UK Biobank. We use this below by loading the whole list into
#   the R session. This requires a fair amount of memory and should preferably used in a cluster.
#     https://drive.google.com/file/d/10nJEydJ_FpcRYzzZYq1xEWk8qtSlQl1X/view?usp=sharing
# 
#
# We also provide a smaller toy dataset with the data of 5 traits only.
# This is a single RData file with the objects explained above. See:
#
# Thus, the following lines that load the data are irrelevant for this
# simple dataset.

# Set the working dir to where the data files are and load the data
setwd("~/Desktop/causal_inference_projects/ms2/data")
# Source the code with the functions
source("https://raw.githubusercontent.com/rivas-lab/cGAUGE/master/R/cGAUGE.R?token=ADT2FBSPNAMHJZPBVDS2MKK5J5JJM")

# Load the real data from the paper

# The are the trait annotation
# Their raw codes are used in the Rivaslab GBE and are not human
# readable. Here we create the map from codes to their description
rivaslab_pheno_codes_file = "icdinfo.txt"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=F)
rownames(rivaslab_codes) = rivaslab_codes[,1]
icd2name = rivaslab_codes[,3];names(icd2name) = rownames(rivaslab_codes)

# Define the paths to the data files
conditional_indep_tests = "genetic_CI_tests_results.RData"
skeleton_file = "Gs_skeleton.RData"
maf_file = "genotypes.frq"
bim_file = "genotypes.bim"
gwas_res_data = "single_gwas_res.RData"

load(conditional_indep_tests)
load(gwas_res_data)
load(skeleton_file)
# Define the trait skeleton
skeleton_pmax = pmax_network
# Use the clumped/prune variant lists (per trait)
pruned_snp_list = cl_unified_list
# Define the MAF for the downstream analysis
MAF = 0.05
# Read the MAF and location info
mafs = read.table(maf_file,stringsAsFactors = F,header=T)
bim = read.table(bim_file,stringsAsFactors = F)
# Filter by MHC
mhc_snps = bim[bim[,1]==6 & bim[,4]>20000000 & bim[,4]<35000000,2]
our_snps = mafs$SNP[mafs$MAF >= MAF]
pruned_snp_list = intersect(pruned_snp_list,our_snps)
pruned_snp_list = intersect(pruned_snp_list,rownames(snp_matrix))
pruned_snp_list = setdiff(pruned_snp_list,mhc_snps)
# restrict the analysis to use the filters defined above
GWAS_Ps = snp_matrix[pruned_snp_list,]
# Add missing names to the mapping vector
missing_names = setdiff(colnames(GWAS_Ps),names(icd2name))
icd2name[missing_names] = missing_names
icd2name = icd2name[colnames(GWAS_Ps)]

# # Create the toy dataset - BMI, SBP, LDL, MI, angina
# toy_dataset_traits = c("INI21001","INI4080","LDL_direct","HC326","HC132")
# is.element(toy_dataset_traits,set=names(trait_pair_pvals))
# for(tr in toy_dataset_traits){
#   trait_pair_pvals[[tr]] = trait_pair_pvals[[tr]][toy_dataset_traits]
# }
# trait_pair_pvals = trait_pair_pvals[toy_dataset_traits]
# skeleton_pmax = skeleton_pmax[toy_dataset_traits,toy_dataset_traits]
# GWAS_Ps = GWAS_Ps[,toy_dataset_traits]
# code2pruned_list = code2pruned_list[toy_dataset_traits]
# code2clumped_list = code2clumped_list[toy_dataset_traits]
# icd2name = icd2name[toy_dataset_traits]
# icd2name["LDL_direct"] = "LDL_direct"
# sum_stat_matrix = sum_stat_matrix[,toy_dataset_traits]
# sum_stat_se_matrix = sum_stat_se_matrix[,toy_dataset_traits]
# save(trait_pair_pvals,GWAS_Ps,skeleton_pmax,
#      pruned_snp_list,code2clumped_list,code2pruned_list,icd2name,
#      sum_stat_se_matrix, sum_stat_matrix,
#      file = "toy_dataset.RData")
# # Add missing names to the mapping vector
# missing_names = setdiff(colnames(GWAS_Ps),names(icd2name))
# icd2name[missing_names] = missing_names
# icd2name = icd2name[colnames(GWAS_Ps)]
# gc()

# Load the toy dataset
load("./toy_dataset.RData")

# Make sure these are installed:
library(MendelianRandomization)
library(limma)

p1 = 1e-07
p2 = 0.001

# Define the skeletons
all_skel_ps = skeleton_pmax[lower.tri(skeleton_pmax)]
skeleton_pthr = max(all_skel_ps[all_skel_ps<p1],na.rm = T)
G_t = skeleton_pmax < skeleton_pthr
diag(G_t) = F;mode(G_t)="numeric"

# Run cGAUGE's constraint-based methods
cGAUGE_G_VT = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,p1,p2)
cGAUGE_G_VT = cGAUGE_G_VT[[1]]
cGAUGE_DepEmerge = DepEmerge(GWAS_Ps,trait_pair_pvals,p1,p2)
cGAUGE_EdgeSep = EdgeSep(GWAS_Ps,G_t,trait_pair_pvals,p1,p2)

# Run cGAUGE's MR analysis
meta_anal_res = run_pairwise_pval_combination_analyses(cGAUGE_G_VT,
    GWAS_Ps,pruned_lists=code2clumped_list,maxp=0.001)
# Analysis 5.2: Various MR methods
mr_anal_res = list(
  "Egger" = run_pairwise_mr_analyses(cGAUGE_G_VT,sum_stat_matrix,sum_stat_se_matrix,
                                     pleio_size=1,pruned_lists=code2clumped_list,func=mr_egger,robust=T),
  "IVW" = run_pairwise_mr_analyses(cGAUGE_G_VT,sum_stat_matrix,sum_stat_se_matrix,
                                   pleio_size=1,pruned_lists=code2clumped_list,func=mr_ivw,robust=T)
)
cleaned_Egger_res = combine_mm_mr_analyses(meta_anal_res,mr_anal_res[["Egger"]],
                                           pi1_thr=2,p_h_thr = -1,minIVs = 3)
colnames(cleaned_Egger_res) = c("tr1->","tr2","p_MR","Est","pi1","numIVs")
cleaned_Egger_res_non_pleio = clean_non_pleio_pairs(cleaned_Egger_res,G_t)
is_non_pleio = is.element(rownames(cleaned_Egger_res),set=rownames(cleaned_Egger_res_non_pleio))
cleaned_Egger_res = cbind(cleaned_Egger_res,is_non_pleio)
effect_direction = rep("Up",nrow(cleaned_Egger_res))
effect_direction[as.numeric(cleaned_Egger_res[,"Est"])<0] = "Down"
cleaned_Egger_res = cbind(cleaned_Egger_res,effect_direction)
# map to readable names
cleaned_Egger_res[,1] = icd2name[cleaned_Egger_res[,1]]
cleaned_Egger_res[,2] = icd2name[cleaned_Egger_res[,2]]
# cleaned_Egger_res - the final table with the MR-Egger results

# Use if cGAUGE was run on the real data (i.e., NOT the toy example)
# Here we load the data file with the results for p1=0.001 and p2 = 1e-07
# that were used as a part of the interpretation in the paper
# We check the agreement between the results.
# 
# load( paste(out_path,"three_rule_analysis_",p1,"_",p2,".RData",sep=""))
# all(table(G_it == cGAUGE_G_VT))
# length(cGAUGE_EdgeSep) == length(detected_cis_per_edge)
# all(newly_formed_sigs == cGAUGE_DepEmerge[[2]])
# for(i in 1:length(detected_cis_per_edge)){
#   if(is.element("list",set=sapply(detected_cis_per_edge[[i]],class))){
#     print(all(unlist(cGAUGE_EdgeSep[[i]]) == unlist(detected_cis_per_edge[[i]][-1])))
#   }
#   else{
#     print(all(unlist(cGAUGE_EdgeSep[[i]]) == unlist(detected_cis_per_edge[[i]])))
#   }
# }



