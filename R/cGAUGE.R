# This is an implementation of the cGAUGE procedures.
# It generally takes as input several parameters and preprocessed data containing
# both GWAS results, a skeleton of the traits (G_T), and pairwise conditional independence
# results.
# To run cGAUGE the data must be formatted as discussed in the input of each of the functions
# below.
# skeleton_pmax

# Try loading required packages
try({library(limma)})
try({library(MendelianRandomization)})

#' Use conditional independencies to infer the variants-traits skeleton.
#' 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param P1 A number. A threshold used to define significant association.
#' @param P2 A number. A threshold used to define independent association (when p >P2).
#' @return A list wit two objects: a binary matrix with the skeleton and the list of the separating sets that rendered associations independent.
extract_skeleton_G_VT<-function(GWAS_Ps,trait_pair_pvals,P1,P2){
  G_VT = GWAS_Ps <= P1
  num_ivs = nrow(GWAS_Ps)
  num_traits = ncol(GWAS_Ps)
  iv_trait_sepsets = list()
  for(i in 1:num_traits){
    tr1 = colnames(GWAS_Ps)[i]
    print(paste("#####",tr1,icd2name[tr1]))
    for(j in 1:num_traits){
      if(i==j){next}
      tr2 = colnames(GWAS_Ps)[j]
      curr_p_mat = trait_pair_pvals[[tr1]][[tr2]][rownames(GWAS_Ps),]
      for (testname in colnames(curr_p_mat)){
        # Do we have significant SNPs that become non-sig?
        curr_ci_test_results = curr_p_mat[,testname] >= P2
        curr_ci_test_results[is.na(curr_ci_test_results)]=F
        curr_sep_snps = rownames(curr_p_mat)[curr_ci_test_results]
        curr_new_sep_snps = curr_sep_snps[G_VT[curr_sep_snps,tr1]]
        G_VT[curr_sep_snps,tr1] = F
        num_new_seps = length(curr_new_sep_snps)
        if(num_new_seps>0){print(paste(num_new_seps, "new sepsets were found for tr2:",tr2,icd2name[tr2]))}
        curr_sepset = tr2
        if(num_new_seps > 0){
          if(is.null(iv_trait_sepsets[[tr1]])){iv_trait_sepsets[[tr1]] = list()}
          for(snp in curr_new_sep_snps){
            iv_trait_sepsets[[tr1]][[snp]] = union(iv_trait_sepsets[[tr1]][[snp]],curr_sepset)
          }
        }
      }
    }
  }
  return(list(G_VT,iv_trait_sepsets))
}

#' Search for newly-formed associations when conditioned on phenotypes.
#' 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param trait_pair_pvals. A named list. Each element is a list. Element [[tr1]][[tr2]] in the list is the conditional independence results for trait 1 conditioned on trait 2.
#' @param P1 A number. A threshold used to define significant association.
#' @param text_col_name a string. The column name to take for the pairwise p-value (i.e., from each element of trait_pair_pvals)
#' @return A list with two matrices. In each matrix, each row has four elements: the snp id, trait1, trait2, p-value of trait 1, p-value of trait 2. The first matrix contains all emerging associations. The second matrix contains only the variants that have significant GWAS association with tr2 (at p1).
DepEmerge<-function(GWAS_Ps,trait_pair_pvals,P1,P2,text_col_name="test3"){
  num_traits = ncol(GWAS_Ps)
  iv_trait_pairs_that_become_dep = c()
  for(i in 1:num_traits){
    tr1 = colnames(GWAS_Ps)[i]
    ps_tr1 = GWAS_Ps[,tr1]
    indep_tr1_snps = rownames(GWAS_Ps)[ps_tr1 >= P2]
    for(j in 1:num_traits){
      if(i==j){next}
      tr2 = colnames(GWAS_Ps)[j]
      curr_p_mat = trait_pair_pvals[[tr1]][[tr2]]
      if(!is.element(text_col_name,set=colnames(curr_p_mat))){next}
      curr_sig_snps = rownames(curr_p_mat)[curr_p_mat[,text_col_name] <= P1]
      new_sig_prev_nonsig_snps = intersect(curr_sig_snps,indep_tr1_snps) # having snps here may point out to t2->t1
      for(new_sig_s in new_sig_prev_nonsig_snps){
        iv_trait_pairs_that_become_dep = rbind(
          iv_trait_pairs_that_become_dep,c(new_sig_s,tr1,tr2,GWAS_Ps[new_sig_s,tr1],GWAS_Ps[new_sig_s,tr2]))
      }
    }
  }
  ps = as.numeric(iv_trait_pairs_that_become_dep[,5])
  newly_formed_sigs = iv_trait_pairs_that_become_dep[ps<p1,]
  return(list("all"=iv_trait_pairs_that_become_dep,"only_tr2_genetic_variants"=newly_formed_sigs))
}

#' Get all cases of disappearing correlations (based on p1,p2)
#'
#' @details For each analyzed (x,y) G_t skeleton we compute the number of variants that are associated with x and y but not with y given x. We also compute the fraction of these cases out of all of x's variants. 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param G_t A binary matrix. TRUE values represent trait skeleton edges.
#' @param trait_pair_pvals. A named list. Each element is a list. Element [[tr1]][[tr2]] in the list is the conditional independence results for trait 1 conditioned on trait 2.
#' @param P1 A number. A threshold used to define significant association.
#' @param text_col_name a string. The column name to take for the pairwise p-value (i.e., from each element of trait_pair_pvals)
#' @param pheno_names A named character vector. 
#' @return A list with two matrices. In each matrix, each row has four elements: the snp id, trait1, trait2, p-value of trait 1, p-value of trait 2. The first matrix contains all emerging associations. The second matrix contains only the variants that have significant GWAS association with tr2 (at p1).
EdgeSep<-function(GWAS_Ps,G_t,trait_pair_pvals,p1,p2,text_col_name="test3",pheno_names=NULL){
  detected_cis_per_edge = list()
  for(tr1 in colnames(GWAS_Ps)){
    for(tr2 in colnames(GWAS_Ps)){
      if(tr1==tr2){next}
      # Go over skeleton edges only
      if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]==0){next}
      # Check which variants associated with both tr1 and tr2 lose the association with tr2
      ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][pruned_snp_list,text_col_name]
      currN = sum(GWAS_Ps[,tr1] < p1,na.rm = T)
      curr_test_inds = GWAS_Ps[,tr2] < p1 & GWAS_Ps[,tr1] < p1 & ps_with_tr2_cond_tr1 > p2
      if(sum(curr_test_inds,na.rm = T)==0){next}
      if(!is.null(pheno_names)){
        currname = paste(tr1,"cause_of",tr2,pheno_names[tr1],"cause_of",pheno_names[tr2],sep = ";")
      }
      else{
        currname = paste(tr1,"cause_of",tr2,sep = ";")
      }
      curr_test_inds = which(curr_test_inds)
      detected_cis_per_edge[[currname]] = list(num_tests = currN,variants=rownames(GWAS_Ps)[curr_test_inds])
    }
  }
  # Filter out intersection between reverse edges
  rev_edge_intersection_bias_sign=list()
  for(nn in names(detected_cis_per_edge)){
    arr = strsplit(nn,split=";")[[1]]
    arr2 = arr[c(3:1,6:4)]
    e1 = paste(arr,collapse = ";")
    e2 = paste(arr2,collapse = ";")
    if(!is.element(e2,set=names(detected_cis_per_edge))){next}
    g1 = detected_cis_per_edge[[e1]]$variants
    g2 = detected_cis_per_edge[[e2]]$variants
    print(paste(e1,length(intersect(g1,g2))))
    if(length(intersect(g1,g2))>0){
      rev_edge_intersection_bias_sign[[paste(arr[1],arr[3],sep=";")]]
      detected_cis_per_edge[[e1]]$variants = setdiff(detected_cis_per_edge[[e1]],intersect(g1,g2))
      detected_cis_per_edge[[e2]]$variants = setdiff(detected_cis_per_edge[[e2]],intersect(g1,g2))
    }
  }
  return(detected_cis_per_edge)
}

#' Go over all trait pairs and run MR
#' 
#' @param G_VT A binary matrix. The instruments-trait skeleton from which genetic variants are chosen. 
#' @param sum_stats A numeric matrix. Rows are instruments, columns are phenotypes. Values are effect sizes.
#' @param sum_stats_se A numeric matrix. Rows are instruments, columns are phenotypes. Values are effect size standard errors.
#' @param pleio_size A number. The maximal number of phenotypes added per variant 
#' @param pruned_lists A named list. Contains a set of pruned or clumped variants per phenotype. The variant names should fit the rownames in the matrices above.
#' @param ... Additional parameters to run_single_mr_analysis.
#' @return A matrix. A row for each analyzed pair. First elements are phenotype 1 (cause), phenotype 2. Then the MR results (depend on the MR method used). Last element is the number of variants of phenotype 1 that were used in the analysis.
run_pairwise_mr_analyses<-function(G_VT,sum_stats,sum_stats_se,
                                   pleio_size=1, pruned_lists=NULL,...){
  trait_pairs_analysis = c()
  traits = colnames(G_VT)
  num_tests = 0
  iv2num_traits = rowSums(G_VT)
  for(tr1 in traits){
    iv2num_traits = rowSums(G_VT,na.rm=T)
    ivs = G_VT[,tr1]==1 & iv2num_traits <= pleio_size
    ivs[is.na(ivs)]=F
    ivs = rownames(G_VT)[ivs]
    if(!is.null(pruned_lists)){ivs = intersect(ivs,pruned_lists[[tr1]])}
    if(length(ivs)<3){next}
    for(tr2 in traits){
      if(tr1==tr2){next}
      try({ # required as some MR methods may fail
        curr_mr_res = run_single_mr_analysis(ivs,tr1,tr2,sum_stats,sum_stats_se,...);
        trait_pairs_analysis = rbind(trait_pairs_analysis,c(tr1,tr2,curr_mr_res,length(ivs)))        
      })
    }
  }
  return(trait_pairs_analysis)
}

#' Go over all trait pairs compute the proportion of non-null p-values.
#'  
#' @param G_VT A binary matrix. The instruments-trait skeleton from which genetic variants are chosen. 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param pleio_size A number. The maximal number of phenotypes added per variant 
#' @param pruned_lists A named list. Contains a set of pruned or clumped variants per phenotype. The variant names should fit the rownames in the matrices above.
#' @return A matrix. A row for each analyzed pair. First elements are phenotype 1 (cause), phenotype 2. Then the estimated proportion and the number of variants of phenotype 1 that were used in the analysis.
run_pairwise_pval_combination_analyses<-function(G_VT,GWAS_Ps,pleio_size=1,pruned_lists=NULL,maxp=0.001){
  trait_pairs_analysis = c()
  traits = colnames(G_VT)
  num_tests = 0
  iv2num_traits = rowSums(G_VT)
  for(tr1 in traits){
    ivs = G_VT[,tr1]==1 & iv2num_traits <= pleio_size
    ivs = rownames(G_VT)[ivs]
    if(!is.null(pruned_lists)){ivs = intersect(ivs,pruned_lists[[tr1]])}
    if(length(ivs)<1){next}
    for(tr2 in traits){
      if(tr1==tr2){next}
      curr_ps = GWAS_Ps[ivs,tr2]
      curr_ps = pmax(curr_ps,maxp)
      curr_ps = curr_ps[!is.na(curr_ps)]
      if(length(curr_ps)==0){next}
      curr_prop = 1-propTrueNull(curr_ps)
      trait_pairs_analysis = rbind(trait_pairs_analysis,c(tr1,tr2,curr_prop,length(ivs)))
    }
  }
  return(trait_pairs_analysis)
}

# For NonPleioMR we need two types of filters: one using analysis thresholds and one using skeleton Edges
combine_mm_mr_analyses<-function(mm,mr,p_thr=0.1,p_h_thr=0,pi1_thr=0.5,minIVs=5){
  names1 = paste(mm[,1],mm[,2],sep="->")
  names2 = paste(mr[,1],mr[,2],sep="->")
  rownames(mm) = names1;rownames(mr)=names2
  inds =intersect(names1,names2)
  mm=mm[inds,];mr=mr[inds,]
  # filters of the results
  corrected_ps = p.adjust(as.numeric(as.character(mr[,3])),method="fdr")
  filter1 = corrected_ps <= p_thr
  pi1s = as.numeric(as.character(mm[,4]))
  filter2 = pi1s >= pi1_thr
  p_hs = p.adjust(as.numeric(as.character(mr[,4])))
  filter3 = p_hs >= p_h_thr
  filter4 = as.numeric(as.character(mr[,ncol(mr)])) >= minIVs
  res_inds = (filter2 & filter4) | (filter1 & filter4 & filter3)
  res_inds[is.na(res_inds)]=F
  res = cbind(mr[res_inds,c(1:3,5)],mm[res_inds,3:4])
  colnames(res) = c("Cause","Effect","P","Est","Pi1","NumIVs")
  return(res)
}

clean_non_pleio_pairs<-function(res,G_t){
  to_keep = rep(T,nrow(res))
  for(i in 1:nrow(res)){
    tr1 = res[i,1];tr2=res[i,2]
    if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]>0){to_keep[i]=F}
  }
  non_pleio_res = res[to_keep,]
  return(non_pleio_res)
}

# Helper functions for running MR using MendelianRandomization
run_single_mr_analysis<-function(snpset,tr1,tr2,X,Y,func=mr_egger,...){
  mr_in = mr_input(X[snpset,tr1],Y[snpset,tr1],X[snpset,tr2],Y[snpset,tr2])
  xx = func(mr_in,...)
  p = 1
  if(is.element("Pvalue.Est",set=slotNames(xx))){p=xx@Pvalue.Est}
  if(is.element("Pvalue",set=slotNames(xx))){p=xx@Pvalue}
  p_het = 1;Q=0;I2=100
  if(is.element("Heter.Stat",set=slotNames(xx))){
    p_het = xx@Heter.Stat[2]
    Q = xx@Heter.Stat[1]
  }
  est = 0
  if(is.element("Estimate",set=slotNames(xx))){est = xx@Estimate}
  return(c(p,p_het,est,Q))
}

