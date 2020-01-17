# This is an implementation of the cGAUGE procedures.
# It generally takes as input several parameters and preprocessed data containing
# both GWAS results, a skeleton of the traits (G_T), and pairwise conditional independence
# results.
# To run cGAUGE the data must be formatted as discussed in the input of each of the functions
# below.
# skeleton_pmax

# Try loading required packages
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

#' Use conditional independencies to infer the variants-traits skeleton.
#' 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param trait_pair_pvals. A named list. Each element is a list. Element [[tr1]][[tr2]] in the list is the conditional independence results for trait 1 conditioned on trait 2.
#' @param P1 A number. A threshold used to define significant association.
#' @param P2 A number. A threshold used to define independent association (when p >P2).
#' @param test_columns A vector. Contains either the index or name of the column with the P-values. May have more than one index/name. If NULL, we assume that all columns are p-values of different tests.
#' @return A list wit two objects: a binary matrix with the skeleton and the list of the separating sets that rendered associations independent.
extract_skeleton_G_VT<-function(GWAS_Ps,trait_pair_pvals,P1,P2,test_columns = NULL){
  G_VT = GWAS_Ps <= P1
  num_ivs = nrow(GWAS_Ps)
  num_traits = ncol(GWAS_Ps)
  iv_trait_sepsets = list()
  for(i in 1:num_traits){
    tr1 = colnames(GWAS_Ps)[i]
    print(paste("#####",tr1))
    for(j in 1:num_traits){
      if(i==j){next}
      tr2 = colnames(GWAS_Ps)[j]
      curr_p_mat = trait_pair_pvals[[tr1]][[tr2]][rownames(GWAS_Ps),]
      curr_columns = test_columns
      if(is.null(curr_columns)){
        curr_columns = colnames(curr_p_mat)
      }
      if(is.character(curr_columns)){
        curr_columns = intersect(curr_columns,colnames(curr_p_mat))
      }
      for (testname in curr_columns){
        # Do we have significant SNPs that become non-sig?
        curr_ci_test_results = curr_p_mat[,testname] >= P2
        curr_ci_test_results[is.na(curr_ci_test_results)]=F
        curr_sep_snps = rownames(curr_p_mat)[curr_ci_test_results]
        curr_new_sep_snps = curr_sep_snps[G_VT[curr_sep_snps,tr1]]
        G_VT[curr_sep_snps,tr1] = F
        num_new_seps = length(curr_new_sep_snps)
        if(num_new_seps>0){print(paste(num_new_seps,
               "new sepsets were found for tr2:",tr2))}
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
EdgeSep<-function(GWAS_Ps,G_t,trait_pair_pvals,p1,p2,pruned_snp_lists = NULL,
                  text_col_name="test3",pheno_names=NULL){
  detected_cis_per_edge = list()
  for(tr1 in colnames(GWAS_Ps)){
    for(tr2 in colnames(GWAS_Ps)){
      if(tr1==tr2){next}
      tr1_tr2_ivs = rownames(GWAS_Ps)[GWAS_Ps[,tr1]<p1 & GWAS_Ps[,tr2]<p1]
      if(!is.null(pruned_snp_lists)){
        tr1_tr2_ivs = intersect(pruned_snp_lists[[tr1]],tr1_tr2_ivs)
      }
      # Go over skeleton edges only
      if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]==0){next}
      # Check which variants associated with both tr1 and tr2 lose the association with tr2
      ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][tr1_tr2_ivs,text_col_name]
      currN = sum(GWAS_Ps[,tr1] < p1,na.rm = T)
      curr_test_inds = ps_with_tr2_cond_tr1 > p2
      if(sum(curr_test_inds,na.rm = T)==0){next}
      if(!is.null(pheno_names)){
        currname = paste(tr1,"cause_of",tr2,pheno_names[tr1],
                         "cause_of",pheno_names[tr2],sep = ";")
      }
      else{
        currname = paste(tr1,"cause_of",tr2,sep = ";")
      }
      curr_test_inds = which(curr_test_inds)
      detected_cis_per_edge[[currname]] = list(num_tests = currN,
                          variants=tr1_tr2_ivs[curr_test_inds])
    }
  }
  
  # Filter out intersection between reverse edges
  rev_edge_intersection_bias_sign=list()
  pairs_as_str = gsub(";cause_of;","",names(detected_cis_per_edge))
  for(nn in names(detected_cis_per_edge)){
    arr = strsplit(nn,split=";")[[1]]
    curr_expo = arr[1];curr_out = arr[3]
    rev_edge_regex = paste(curr_out,curr_expo,"$",sep="")
    if(!is.null(pheno_names)){
      rev_edge_regex = paste(curr_out,curr_expo,";",sep="")
    }
    rev_edge_ind = which(grepl(rev_edge_regex,pairs_as_str))
    if(length(rev_edge_ind)<1){next}
    e1 = nn
    e2 = names(detected_cis_per_edge)[rev_edge_ind]
    g1 = detected_cis_per_edge[[e1]]$variants
    g2 = detected_cis_per_edge[[e2]]$variants
    if(length(intersect(g1,g2))>0){
      rev_edge_intersection_bias_sign[[paste(arr[1],arr[3],sep=";")]]
      detected_cis_per_edge[[e1]]$variants = setdiff(detected_cis_per_edge[[e1]],intersect(g1,g2))
      detected_cis_per_edge[[e2]]$variants = setdiff(detected_cis_per_edge[[e2]],intersect(g1,g2))
    }
  }
  return(detected_cis_per_edge)
}

#' Statistical test analysis for separation along skeleton edges
#'
#' @details 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param G_t A binary matrix. TRUE values represent trait skeleton edges.
#' @param trait_pair_pvals. A named list. Each element is a list. Element [[tr1]][[tr2]] in the list is the conditional independence results for trait 1 conditioned on trait 2.
#' @param text_col_name a string. The column name to take for the pairwise p-value (i.e., from each element of trait_pair_pvals)
#' @param test_func a function. Takes two (paired) p-value vectors and returns a statistic (typically a p-value) testing if there is evidence for "disappearing associations"
#' @return A matrix with three columns.
EdgeSepTest<-function(GWAS_Ps,G_t,trait_pair_pvals,text_col_name="test3",
                      test_func = simple_lfdr_test,...){
  edge_sep_tests = c()
  for(tr1 in colnames(GWAS_Ps)){
    for(tr2 in colnames(GWAS_Ps)){
      if(tr1==tr2){next}
      if(is.na(G_t[tr1,tr2]) || G_t[tr1,tr2]==0){next}
      p1 = GWAS_Ps[,tr2]
      if(is.null(rownames(GWAS_Ps))){
        ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][,text_col_name]
      }
      else{
        ps_with_tr2_cond_tr1 = trait_pair_pvals[[tr2]][[tr1]][rownames(GWAS_Ps),text_col_name]
      }
      test1 = test_func(p1,ps_with_tr2_cond_tr1,...)
      edge_sep_tests = rbind(edge_sep_tests,c(tr1,tr2,test1))
    }
  }
  colnames(edge_sep_tests) = c("trait1","trait2","pval:trait1->trait2")
  edge_sep_tests = as.data.frame(edge_sep_tests)
  edge_sep_tests[[3]] = as.numeric(as.character(edge_sep_tests[[3]]))
  return(edge_sep_tests)
}

#' Remove non-minimal separating sets (naive implementation, for QC).
#'
#' @param sepsets a list of sets
#' @return a subset of the input lits
#' @details remove all sets such that there exists another set that is contained in them; useful for cleaning instrument sets
remove_non_minimal_sepsets_naive<-function(sepsets){
  if(length(sepsets)==0){return(list())}
  if(length(sepsets)==1){return(l)}
  to_rem = rep(F,length(sepsets))
  # go over all pairs
  for(i in 2:length(sepsets)){
    for(j in 1:(i-1)){
      set1 = sepsets[[i]]
      set2 = sepsets[[j]]
      if(all(set1 %in% set2)){
        to_rem[j] = T
      }
      if(all(set2 %in% set1)){
        to_rem[i] = T
      }
    }
  }
  return(sepsets[!to_rem])
}

# a faster version of the above
#' Remove non-minimal separating sets.
#' 
#' @param sepsets a list of sets
#' @return a subset of the input lits
#' @details remove all sets such that there exists another set that is contained in them; useful for cleaning instrument sets
remove_non_minimal_sepsets<-function(l){
  if(length(l)==0){return(list())}
  if(length(l)==1){return(l)}
  
  sizes = sapply(l,length)
  l = l[order(sizes)]
  sizes = sapply(l,length)
  i = 1
  while(i < length(l)){
    set1 = l[[i]]
    to_rem = rep(F,length(l))
    for(j in (i+1):length(l)){
      if(sizes[i]==sizes[j]){next}
      set2 = l[[j]]
      if(all(set1 %in% set2)){
        to_rem[j] = T
      }
    }
    i = i+1
    l = l[!to_rem]
    sizes = sizes[!to_rem]
  }
  
  return(l)
}

# tr1 = "T3"
# tr2 = "T13"
# p1 = GWAS_Ps[,tr2]
# p2 = trait_pair_pvals[[tr2]][[tr1]][,1]
# z1 = -qnorm(p1)
# z2 = -qnorm(p2)
# # z2[z2>z1+2] = z2[z2>z1+2] -1
# # ind = which(z1>4 & z2 > 4)[1:2]
# # z2[ind] = 0
# # par(mar = c(9,9,9,9))
# plot(x=z1,y=z2,pch=20, cex = 1.2,
#      xlab = "Assoc with Y (z-scores)",
#      ylab = "Assoc with Y given X",
#      cex.lab = 1.4)
# abline(0,1,xpd=F,lwd=2,lty=2,col="red")

#' Go over all trait pairs and run MR from a variants x traits graph
#'
#' @param G_VT A binary matrix. The instruments-trait skeleton from which genetic variants are chosen.
#' @param sum_stats A numeric matrix. Rows are instruments, columns are phenotypes. Values are effect sizes.
#' @param sum_stats_se A numeric matrix. Rows are instruments, columns are phenotypes. Values are effect size standard errors.
#' @param pleio_size A number. The maximal number of phenotypes added per variant
#' @param pruned_lists A named list. Contains a set of pruned or clumped variants per phenotype. The variant names should fit the rownames in the matrices above.
#' @param ... Additional parameters to run_single_mr_analysis.
#' @return A matrix. A row for each analyzed pair. First elements are phenotype 1 (cause), phenotype 2. Then the MR results (depend on the MR method used). Last element is the number of variants of phenotype 1 that were used in the analysis.
run_pairwise_mr_analyses<-function(G_VT,sum_stats,sum_stats_se,
                                   pleio_size=1, minIVs = 5,pruned_lists=NULL,...){
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
    if(length(ivs)<minIVs){next}
    for(tr2 in traits){
      if(tr1==tr2){next}
      try({ # required as some MR methods may fail
        curr_mr_res = run_single_mr_analysis(ivs,tr1,tr2,sum_stats,sum_stats_se,...);
        trait_pairs_analysis = rbind(trait_pairs_analysis,c(tr1,tr2,curr_mr_res,length(ivs)))
      })
    }
  }
  if(!is.null(dim(trait_pairs_analysis))){
    colnames(trait_pairs_analysis) = c("Exposure","Outcome","p","p_het","est","Q","NumIVs")
    trait_pairs_analysis = as.data.frame(trait_pairs_analysis)
    for(j in 3:ncol(trait_pairs_analysis)){
      trait_pairs_analysis[[j]] = as.numeric(as.character(trait_pairs_analysis[[j]]))
    }
  }
  return(trait_pairs_analysis)
}

#' Go over all trait pairs and run MR using cGAUGE's instrument sets
#' 
#' @param sum_stats A numeric matrix. Rows are instruments, columns are phenotypes. Values are effect sizes.
#' @param sum_stats_se A numeric matrix. Rows are instruments, columns are phenotypes. Values are effect size standard errors.
#' @param iv_sets a list of lists, entry [[a]][[b]] is the snp ids or indices for the MR of a->b 
#' @param minIVs a number, ignore iv sets with less than minIVs snps
#' @param ... Additional parameters to run_single_mr_analysis.
#' @return A matrix. A row for each analyzed pair. First elements are phenotype 1 (cause), phenotype 2. Then the MR results (depend on the MR method used). Last element is the number of variants of phenotype 1 that were used in the analysis.
run_pairwise_mr_analyses_with_iv_sets<-function(sum_stats,sum_stats_se,iv_sets,
                                                minIVs = 5,...){
  trait_pairs_analysis = c()
  traits = colnames(sum_stats)
  num_tests = 0
  for(tr1 in traits){
    for(tr2 in traits){
      if(tr1==tr2){next}
      ivs = iv_sets[[tr1]][[tr2]]
      if(length(ivs)<minIVs){next}
      try({ # required as some MR methods may fail
        curr_mr_res = run_single_mr_analysis(ivs,tr1,tr2,sum_stats,sum_stats_se,...);
        trait_pairs_analysis = rbind(trait_pairs_analysis,c(tr1,tr2,curr_mr_res,length(ivs)))        
      })
    }
  }
  if(!is.null(dim(trait_pairs_analysis))){
    colnames(trait_pairs_analysis) = c("Exposure","Outcome","p","p_het","est","Q","NumIVs")
    trait_pairs_analysis = as.data.frame(trait_pairs_analysis)
    for(j in 3:ncol(trait_pairs_analysis)){
      trait_pairs_analysis[[j]] = as.numeric(as.character(trait_pairs_analysis[[j]]))
    }
  }
  return(trait_pairs_analysis)
}

#' Go over all trait pairs compute the proportion of non-null p-values.
#'  
#' @param iv_sets a list of lists, entry [[a]][[b]] is the snp ids or indices for the MR of a->b 
#' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' @param pruned_lists A named list. Contains a set of pruned or clumped variants per phenotype. The variant names should fit the rownames in the matrices above.
#' @return A matrix. A row for each analyzed trait pair with the estimates of the outcome non-null exposure instruments. We use limma for estimation using the theoretical null.
run_pairwise_pval_combination_analysis_from_iv_sets<-function(iv_sets,GWAS_Ps){
  trait_pairs_analysis = c()
  traits = colnames(GWAS_Ps)
  num_tests = 0
  for(tr1 in traits){
    for(tr2 in traits){
      if(tr1==tr2){next}
      ivs = iv_sets[[tr1]][[tr2]]
      if(length(ivs)<2){next}
      curr_ps = GWAS_Ps[ivs,tr2]
      curr_ps = curr_ps[!is.na(curr_ps)]
      if(length(curr_ps)==0){next}
      limma_lfdr_prop = 1-propTrueNull(curr_ps)
      limma_hist_prop = 1-propTrueNull(curr_ps,method = "hist")
      trait_pairs_analysis = 
        rbind(trait_pairs_analysis,c(tr1,tr2,limma_lfdr_prop,limma_hist_prop,length(ivs)))
    }
  }
  colnames(trait_pairs_analysis) = c("tr1->","tr2","lfdr_pi_1","convest_pi_1","numIVs")
  return(trait_pairs_analysis)
}

#' Combine and filter the pairwise MR and meta-analysis results
#' @param mm a matrix with the meta-analysis results, each row represents a pair
#' @param mr a matrix with the MR results, each row represents a pair
#' @param p_thr numeric, a threshold for the MR adjusted p-values
#' @param adj_method character, the name of the p-value adjustment method, default is BY, see p.adjust for details
#' @param p_h_thr numeric, trait pairs with heterogeneity p-value < p_h_thr are removed, default is 1 - avoids using this filter
#' @param minIVs numeric, trait pairs whose analysis is based on less than this number are excluded
#' @return a matrix with the combined meta-analysis and Mendelian randomization results
combine_mm_mr_analyses<-function(mm,mr,p_thr=0.01,adj_method="BY",
                                 p_h_thr=1,minIVs=5){
  names1 = paste(mm[,1],mm[,2],sep="->")
  names2 = paste(mr[,1],mr[,2],sep="->")
  rownames(mm) = names1;rownames(mr)=names2
  inds =intersect(names1,names2)
  mm=mm[inds,];mr=mr[inds,]
  # filters of the results
  corrected_ps = p.adjust(as.numeric(as.character(mr[,3])),method=adj_method)
  filter1 = corrected_ps <= p_thr
  p_hs = p.adjust(as.numeric(as.character(mr[,4])))
  filter2 = (p_hs >= p_h_thr) | is.na(p_hs)
  filter3 = as.numeric(as.character(mr[,ncol(mr)])) >= minIVs
  res_inds = filter1 & filter2 & filter3
  res_inds[is.na(res_inds)]=F
  res = cbind(mr[res_inds,c(1:3,5)],mm[res_inds,c("lfdr_pi_1","numIVs")])
  colnames(res) = c("tr1->","tr2","p_MR","Est","pi1","numIVs")
  # add a binary column for edge direction
  effect_direction = rep("Up",nrow(res))
  effect_direction[as.numeric(res[,"Est"])<0] = "Down"
  res = cbind(res,effect_direction)
  log10p = -log10(as.numeric(res[,"p_MR"]))
  res = cbind(res,log10p)
  colnames(res) = c("tr1->","tr2","p_MR","Est","pi1","numIVs","EdgeDirection","-log10p_MR")
  for(j in c(1,2)){
    res[[j]] = as.character(res[[j]])
  }
  for(j in c(3,4,5,6,8)){
    res[[j]] = as.numeric(as.character(res[[j]]))
  }
  return(res)
}

#' Add a column that indicates for each pair if it is a skeleton edge or not
#' 
#' @param res a matrix, the first two columns are the exposure (tr1->) and the outcome
#' @param G_t a binary matrix that represents the skeleton of the traits
#' @return res with a new column, is_skeleton_edge
add_is_non_edge_column<-function(res,G_t){
  is_skeleton_edge = rep(F,nrow(res))
  for(i in 1:nrow(res)){
    tr1 = res[i,1];tr2=res[i,2]
    if(!is.na(G_t[tr1,tr2]) && G_t[tr1,tr2]>0){
      is_skeleton_edge[i]=T
    }
  }
  res = cbind(res,is_skeleton_edge)
  return(res)
}

add_edgesep_res_column<-function(mr_res,edgesep_res,G_t){
  # filter out rows that are not edges in G_t
  to_rem = rep(F,nrow(edgesep_res))
  for(i in 1:nrow(edgesep_res)){
    curre = G_t[edgesep_res[i,1],edgesep_res[i,2]]
    if(is.na(curre) || curre == 0){
      to_rem[i] = T
    }
  }
  edgesep_res = edgesep_res[!to_rem,]
  
  names1 = paste(edgesep_res[,1],edgesep_res[,2],sep="->")
  names2 = paste(mr_res[,1],mr_res[,2],sep="->")
  rownames(edgesep_res) = names1
  rownames(mr_res)=names2
  inds =intersect(names1,names2)
  edge_sep_col = rep(NA,length(names2))
  names(edge_sep_col) = names2
  edge_sep_col[inds] = edgesep_res[inds,3]
  mr_res = cbind(mr_res,edge_sep_col)
  mr_res[[ncol(mr_res)]] = as.numeric(as.character(mr_res[[ncol(mr_res)]]))
  return(mr_res)
}

##############################################################
# Helper functions for running MR using MendelianRandomization
##############################################################
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


#' #' Search for newly-formed associations when conditioned on phenotypes.
#' #' 
#' #' @param GWAS_Ps A matrix. Rows are variants and columns are phenotypes. Cells are P-values.
#' #' @param trait_pair_pvals. A named list. Each element is a list. Element [[tr1]][[tr2]] in the list is the conditional independence results for trait 1 conditioned on trait 2.
#' #' @param P1 A number. A threshold used to define significant association.
#' #' @param text_col_name a string. The column name to take for the pairwise p-value (i.e., from each element of trait_pair_pvals)
#' #' @return A list with two matrices. In each matrix, each row has four elements: the snp id, trait1, trait2, p-value of trait 1, p-value of trait 2. The first matrix contains all emerging associations. The second matrix contains only the variants that have significant GWAS association with tr2 (at p1).
#' DepEmerge<-function(GWAS_Ps,trait_pair_pvals,P1,P2,text_col_name="test3"){
#'   num_traits = ncol(GWAS_Ps)
#'   iv_trait_pairs_that_become_dep = c()
#'   for(i in 1:num_traits){
#'     tr1 = colnames(GWAS_Ps)[i]
#'     ps_tr1 = GWAS_Ps[,tr1]
#'     indep_tr1_snps = rownames(GWAS_Ps)[ps_tr1 >= P2]
#'     for(j in 1:num_traits){
#'       if(i==j){next}
#'       tr2 = colnames(GWAS_Ps)[j]
#'       curr_p_mat = trait_pair_pvals[[tr1]][[tr2]]
#'       if(!is.element(text_col_name,set=colnames(curr_p_mat))){next}
#'       curr_sig_snps = rownames(curr_p_mat)[curr_p_mat[,text_col_name] <= P1]
#'       new_sig_prev_nonsig_snps = intersect(curr_sig_snps,indep_tr1_snps) # having snps here may point out to t2->t1
#'       for(new_sig_s in new_sig_prev_nonsig_snps){
#'         iv_trait_pairs_that_become_dep = rbind(
#'           iv_trait_pairs_that_become_dep,c(new_sig_s,tr1,tr2,GWAS_Ps[new_sig_s,tr1],GWAS_Ps[new_sig_s,tr2]))
#'       }
#'     }
#'   }
#'   if(!is.null(dim(iv_trait_pairs_that_become_dep))){
#'     colnames(iv_trait_pairs_that_become_dep) = c("Variant","Exposure","Outcome","PvalExposure","PvalOutcome")
#'     ps = as.numeric(iv_trait_pairs_that_become_dep[,5])
#'     newly_formed_sigs = iv_trait_pairs_that_become_dep[ps<p1,]
#'     return(list("all"=iv_trait_pairs_that_become_dep,"only_tr2_genetic_variants"=newly_formed_sigs))  
#'   }
#'   return(iv_trait_pairs_that_become_dep)
#' }

