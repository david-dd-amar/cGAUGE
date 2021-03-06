
# We assume that summary statistic data were downloaded from the links provided in the Readme
# and put into the working dir.
# The usage of each of these files is given in the code below. 
#
# 1.The results of the single GWAS (one per trait) - single_gwas_res.RData
#       This single GWAS result objects contain the data of variants that had p<0.001 in at least one 
#       of the traits analyzed in the study.
#       The useful objects in this RData are: 
#       snp_P_matrix - GWAS p-values (rows are variants, columns are traits)
#       sum_stat_matrix - GWAS effect sizes (rows are variants, columns are traits)
#       sum_stat_se_matrix - GWAS effect size standard error (rows are variants, columns are traits)
#       cl_unified_list - a list of clumped variant ids per trait (rows are variants, columns are traits)
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

setwd("~/Desktop/causal_inference_projects/data")

rivaslab_pheno_codes_file = "icdinfo.txt"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=F)
rownames(rivaslab_codes) = rivaslab_codes[,1]
icd2name = rivaslab_codes[,3];names(icd2name) = rownames(rivaslab_codes)

# Analysis of 96 traits
conditional_indep_tests = "genetic_CI_tests_results.RData"
skeleton_file = "Gs_skeleton.RData"
maf_file = "genotypes.frq"
bim_file = "genotypes.bim"
gwas_res_data = "single_gwas_res.RData"

# define the output path
out_path = "./res/"
# Load the data
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

# Make sure these are installed:
library(MendelianRandomization)
library(limma)

p1 = 1e-07
p2 = 0.001

p1_seps = get_minimal_separating_sets(sepsets,p1=1e-07)
G_t = pmax_network < p1
diag(G_t) = F;mode(G_t)="numeric"

# Run cGAUGE's methods and compare
G_vt = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,P1=p1,
                             P2=p2,test_columns = c("test2","test3"))[[1]]
ivs = cGAUGE_instrument_filters(G_t,G_vt,p1_seps,GWAS_Ps,p1,code2clumped_list)
uniqueivs = ivs[[1]]
impivs = ivs[[2]]
cGAUGE_EdgeSep = EdgeSep(GWAS_Ps,G_t,trait_pair_pvals,p1,p2)


