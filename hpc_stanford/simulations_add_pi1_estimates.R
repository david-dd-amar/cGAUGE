required_libs = c("igraph","limma","locfdr","optparse")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
           })
}

# Helper functions
run_lm<-function(x,y,z,df){
  if(is.null(z)){
    df = data.frame(x=df[,x],y=df[,y])
  }
  else{
    df = data.frame(x=df[,x],y=df[,y],df[,z])
  }
  model = lm(x~.,data=df)
  coefs = summary(model)$coefficients
  return(coefs[2,])
}
# Helper for running the pi1 analysis
get_pi1_estimates<-function(curr_out_file,p1){
  load(curr_out_file)
  
  phenos = colnames(B)[grepl("^T",colnames(B))]
  ivs = colnames(B)[grepl("^IV",colnames(B))]
  df = data.frame(simulated_data)
  
  # Get all IV-phenotype associations
  GWAS_Ps = matrix(1,length(ivs),length(phenos),
                   dimnames = list(ivs,phenos))
  for(pheno in phenos){
    gwas_res = sapply(ivs,run_lm,x=pheno,z=NULL,df = df)
    GWAS_Ps[,pheno] = gwas_res[4,]
  }
  # get the pi1 estimates
  pi1_estimates = c()
  for(tr1 in colnames(B_distances)){
    for (tr2 in colnames(B_distances)){
      if(tr1==tr2){next}
      currdist = B_distances[tr2,tr1] # distance from tr1 to tr2
      # if(currdist!=-1 && currdist!=1){next}
      raw_ivs = rownames(GWAS_Ps)[GWAS_Ps[,tr1]<p1]
      our_ivs = iv_sets[[tr1]][[tr2]]
      raw_ivs_pi1 = 1-propTrueNull(GWAS_Ps[raw_ivs,tr2])
      our_ivs_pi1 = 1-propTrueNull(GWAS_Ps[our_ivs,tr2])
      pi1_estimates = rbind(pi1_estimates,
                            c(raw_ivs_pi1,our_ivs_pi1,currdist))
    }
  }
  colnames(pi1_estimates) = c("raw_ivs","our_ivs","real_dist")
  pi1_estimates = as.data.frame(pi1_estimates)
  for(j in names(pi1_estimates)){
    pi1_estimates[[j]] = as.numeric(as.character(pi1_estimates[[j]]))
  }
  return(pi1_estimates)
}

option_list <- list( 
  make_option(c("--rdata_file"), action="store",default="",type="character",
              help="RData file with the output of the full_causal_graph_simulations.R script, 
              this script will add a new object called pi1_estimates to this RData"),
  make_option(c("--p1"), action="store", default=0.001,type="double",
              help="p1 value for the naive way of taking instruments")
)

opt <- parse_args(OptionParser(option_list=option_list))
rdata_file = opt$rdata_file
p1 = opt$p1
pi1_estimates = get_pi1_estimates(rdata_file,p1)
load(rdata_file)
success = F
try({
  save(
    opt, # input parameters
    B,Bg,simulated_data,B_distances, # simulated data
    cgauge_mr_results,standard_mr_results, # MR results
    edge_sep_results, # EdgeSep
    edge_sep_results_statTest1, # EdgeSepStatTest using EM
    edge_sep_results_statTest2, # EdgeSepStatTest using lfdrs
    G_it,G_vt,G_t, iv_sets, # Skeletons
    pi1_estimates, # pi1_estimates for each pair
    file = rdata_file
  )
  success = T
})




