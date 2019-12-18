# This script takes a matrix of traits and perform skeleton inference ignoring
# any genetic data (variants or PCs)

auxil_functions_R_code = "/home/users/davidama/repos/CCDfromD/scripts/auxil_functions.R"
ccd_functions_R_code = "/home/users/davidama/repos/CCDfromD/scripts/ccdfromd_functions.R"
library("bnlearn",lib.loc = "~/R/packages/")
library("speedglm",lib.loc = "~/R/packages/")
source(auxil_functions_R_code)
source(ccd_functions_R_code)

# Define the input
input_data = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_data.RData"
out_path = "/oak/stanford/groups/mrivas/users/davidama/april2019_traits_skeleton/"
out_object_to_store_results = "/oak/stanford/groups/mrivas/users/davidama/april2019_Gs_skeleton.RData"

system(paste("mkdir",out_path))
FUNC = run_ci_test_logistic_linear_discrete
NJOB = 2000
depth = 2
pthr = 0.01
to_include = c("sex","age")
load(input_data)

ls_objs = ls()
if(is.element("traits",ls_objs)){
  traits = union(traits,to_include)
}
if(!is.element("traits",ls_objs)){
  traits = names(code2phe_data)
}
n = length(traits)

pmax_network = matrix(NA,nrow=n,ncol=n)
rownames(pmax_network) = traits; colnames(pmax_network) = traits
sepsets = list()
for(nn in colnames(pmax_network)){
  sepsets[[nn]] = list()
}

# Initial filtering to avoid some unecessary runs
for(i in 2:n){
  for(j in 1:(i-1)){
    xx = c(traits[i],traits[j])
    curr_inds = !apply(is.na(alldata[,xx]),1,any)
    if(sum(curr_inds) > 10000){
      pmax_network[i,j] = FUNC(traits[i],traits[j],NULL,data=alldata)
      pmax_network[j,i] = pmax_network[i,j]
    }
  }
  print(i)
}
table(pmax_network < pthr)/2
save(pmax_network,file=paste(out_path,"initial_pmax_network.RData",sep=""))

load(paste(out_path,"initial_pmax_network.RData",sep=""))
# Now run an analysis for each pair
for(i in 2:n){
  for(j in 1:(i-1)){
    if(is.na(pmax_network[i,j])){next}
    if(pmax_network[i,j]>pthr){next}
    tr1 = traits[i]
    tr2 = traits[j]
    curr_job_name = paste(tr1,"_",tr2,sep="")
    
    if(is.element(paste(curr_job_name,".RData",sep=""),set=list.files(out_path))){next}
    
    err_path = paste(out_path,curr_job_name,".err",sep="")
    log_path = paste(out_path,curr_job_name,".log",sep="")
    curr_cmd = paste(
      "module load R\n",
      "Rscript ~/repos/CCDfromD/R_batch_tools/analyze_trait_pair_for_skeleton.R",
      input_data,tr1,tr2,depth,pthr,paste(out_path,curr_job_name,".RData",sep="")
    )
    curr_sh_file = paste(out_path,curr_job_name,".sh",sep="")
    print_sh_file(curr_sh_file,
                  get_sh_prefix_one_node_specify_cpu_and_mem(
                        err_path,log_path,"plink/2.0a1",1,4000,time="04:00:00")
                  ,curr_cmd)
    system(paste("sbatch -x sh-113-15 ",curr_sh_file,'&'))
    while(nrow(get_my_jobs())>=NJOB){Sys.sleep(60)}
  }
  print(i)
}

# for(rr in rownames(get_my_jobs())){system(paste("scancel",rr))}

# Read the results
load(paste(out_path,"initial_pmax_network.RData",sep=""))

for(i in 2:n){
  print(i)
  for(j in 1:(i-1)){
    tr1 = traits[i]
    tr2 = traits[j]
    if(is.na(pmax_network[i,j]) || pmax_network[i,j]>pthr){
      sepsets[[tr1]][[tr2]] = NULL
      sepsets[[tr2]][[tr1]] = NULL
      next
    }
    curr_job_name = paste(tr1,"_",tr2,sep="")
    maxp=NA;sepset=NULL
    try({(load(paste(out_path,curr_job_name,".RData",sep="")))})
    if(is.na(maxp)){
      print(paste("missing analysis for pair:",curr_job_name))
      next
    }
    if(!is.na(pmax_network[i,j])){maxp = max(maxp,pmax_network[i,j])}
    pmax_network[i,j] = maxp
    pmax_network[j,i] = maxp
    sepsets[[tr1]][[tr2]] = sepset
    sepsets[[tr2]][[tr1]] = sepset
  }
}

table(is.na(pmax_network[lower.tri(pmax_network)]))
table(pmax_network[lower.tri(pmax_network)] < pthr)

save(pmax_network,sepsets,file=out_object_to_store_results)

