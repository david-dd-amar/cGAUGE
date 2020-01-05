# Set the session
required_libs = c("optparse","MRPRESSO","MendelianRandomization")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
           })
}

option_list <- list( 
  make_option(c("--gwas_rdata_file"), action="store", 
              default="/oak/stanford/groups/mrivas/users/davidama/april2019_causal_analysis_flow_input.RData",type="character",
              help="RData file with two matrices: sum_stat_matrix and sum_stat_se_matrix"),
  make_option(c( "--ivs_rdata"), action="store", default="",type="character",
              help="RData file with a vector called ivs that contain the instrument ids that should also fit the rownames of the matrices from the first input"),
  make_option(c( "--tr1"), action="store", default="",type="character",
              help="id of exposure"),
  make_option(c( "--tr2"), action="store", default="",type="character",
              help="id of outcome"),
  make_option(c( "--out"), action="store", default="",type="character",
              help="name of output file (RData)")
)

opt <- parse_args(OptionParser(option_list=option_list))
gwas_rdata_file = opt$gwas_rdata_file
ivs_rdata = opt$ivs_rdata
out_file = opt$out
tr1 = opt$tr1
tr2 = opt$tr2

#' A helper function for running MRPRESSO
#' tr1 - exposure
#' tr2 - outcome
#' iv_sets - a list with tr1's instruments for every tr2
mrpresso_wrapper <-function(ivs,GWAS_effects,GWAS_ses,tr1,tr2){
  X = data.frame(E1b=GWAS_effects[ivs,tr1],O1b=GWAS_effects[ivs,tr2],
                 E1sd=GWAS_ses[ivs,tr1],O1sd=GWAS_ses[ivs,tr2])
  try({
    res = mr_presso(BetaOutcome = "O1b", BetaExposure = "E1b", 
                    SdOutcome = "O1sd", SdExposure = "E1sd",data=X,
                    OUTLIERtest=T,
                    DISTORTIONtest = T,
                    NbDistribution = 1000,SignifThreshold = 0.1)
    res$tr1 = tr1
    res$tr2 = tr2
    return(res)
  })
  return(NULL)
}

load(gwas_rdata_file)
load(ivs_rdata)
mrpresso_res = mrpresso_wrapper(ivs,sum_stat_matrix,sum_stat_se_matrix,tr1,tr2)
save(mrpresso_res,file = out_file)



