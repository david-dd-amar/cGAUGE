rivaslab_pheno_codes_file = "/oak/stanford/groups/mrivas//users/magu/repos/rivas-lab/wiki/ukbb/icdinfo/icdinfo.txt"
rivaslab_codes = read.delim(rivaslab_pheno_codes_file,stringsAsFactors = F,header=F)
rownames(rivaslab_codes) = rivaslab_codes[,1]
icd2name = rivaslab_codes[,3];names(icd2name) = rownames(rivaslab_codes)

# Nov 2018: Analysis of 70 traits
out_plink_path = "/oak/stanford/groups/mrivas/users/davidama/nov2018_traits_causal_analysis_flow_results/"
skeleton_file = "/oak/stanford/groups/mrivas/users/davidama/nov2018_Gs_skeleton.RData"
geno_data_path = "/oak/stanford/groups/mrivas/users/davidama/nov2018_traits_genotypes/all_genotypes"
gwas_res_data = "/oak/stanford/groups/mrivas/users/davidama/nov2018_causal_analysis_flow_input.RData"
gwas_res_path = "/oak/stanford/groups/mrivas/users/davidama/gwas_res/"
# Load and define the data
out_path = paste(out_plink_path,"ccd_res/",sep="")
load(paste(out_plink_path,"genetic_CI_tests_results.RData",sep=""))
load(gwas_res_data)
load(skeleton_file)
skeleton_pmax = pmax_network
pruned_snp_list = cl_unified_list
MAF = 0.05
maf_file = paste(geno_data_path,".frq",sep="")
mafs = read.table(maf_file,stringsAsFactors = F,header=T)
bim = read.table(paste(geno_data_path,".bim",sep=""),stringsAsFactors = F)
mhc_snps = bim[bim[,1]==6 & bim[,4]>20000000 & bim[,4]<35000000,2]
our_snps = mafs$SNP[mafs$MAF >= MAF]
pruned_snp_list = intersect(pruned_snp_list,our_snps)
pruned_snp_list = intersect(pruned_snp_list,rownames(snp_matrix))
pruned_snp_list = setdiff(pruned_snp_list,mhc_snps)
GWAS_Ps = snp_matrix[pruned_snp_list,]
all_skel_ps = skeleton_pmax[lower.tri(skeleton_pmax)]
skeleton_pthr = max(all_skel_ps[all_skel_ps<p1],na.rm = T)
G_t = skeleton_pmax < skeleton_pthr
diag(G_t) = F;mode(G_t)="numeric"

library(MendelianRandomization,lib="/home/users/davidama/R/packages")
library(limma,lib="/home/users/davidama/R/packages")

p1 = 1e-07
p2 = 0.001
load(paste(out_path,"three_rule_analysis_",p1,"_",p2,".RData",sep=""))

# Run cGAUGE's methods and compare
cGAUGE_G_VT = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,p1,p2)
cGAUGE_G_VT = cGAUGE_G_VT[[1]]
all(table(G_it == cGAUGE_G_VT))

cGAUGE_DepEmerge = DepEmerge(GWAS_Ps,trait_pair_pvals,p1,p2)
all(newly_formed_sigs == cGAUGE_DepEmerge[[2]])

cGAUGE_EdgeSep = EdgeSep(GWAS_Ps,G_t,trait_pair_pvals,p1,p2)
length(cGAUGE_EdgeSep) == length(detected_cis_per_edge)
for(i in 1:length(detected_cis_per_edge)){
  if(is.element("list",set=sapply(detected_cis_per_edge[[i]],class))){
    print(all(unlist(cGAUGE_EdgeSep[[i]]) == unlist(detected_cis_per_edge[[i]][-1])))
  }
  else{
    print(all(unlist(cGAUGE_EdgeSep[[i]]) == unlist(detected_cis_per_edge[[i]])))
  }
}
