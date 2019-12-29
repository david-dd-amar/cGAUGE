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


option_list <- list( 
  make_option(c("--file"), action="store", default="",type="character",
              help="RData file with a matrix with two columns: first for associations of genetic variables with Y, second is the same as the first but with X added to the conditioned set"),
  make_option(c( "--testName"), action="store", default="em",type="character",
              help="either lfdr or em"),
  make_option(c( "--out"), action="store", default="",type="character",
              help="name of output file")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
rdata_path = opt$file
testName = opt$testName
out_file = opt$out

ps = get(load(rdata_path))
p1 = ps[,1]
p2 = ps[,2]
res = NULL
if(testName == "em"){
  res = univar_mixtools_em(p1,p2)
}
if(testName == "lfdr"){
  res = simple_lfdr_test(p1,p2)
}
save(res,out_file)


