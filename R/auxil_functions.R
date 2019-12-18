required_libs = c("speedglm","bnlearn","pcalg")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
           })
}

get_sh_default_prefix<-function(err="",log="",time="10:00:00"){
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=", time,sep=""),
      "#SBATCH --partition=euan,mrivas,normal,owners",
      "#SBATCH --nodes=1",
      "#SBATCH --cpus-per-task=2",
      "#SBATCH --mem=16000",
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      "module load plink/1.90b5.3"
    )
  )
}

# plink2: plink/2.0a1
get_sh_prefix_one_node_specify_cpu_and_mem<-function(
  err="",log="",plink_pkg = "plink/1.90b5.3",Ncpu,mem_size,time="6:00:00"){
  partition_line = "#SBATCH --partition=euan,mrivas,normal,owners"
  if(mem_size>=256000){
    partition_line = "#SBATCH --partition=bigmem,euan,mrivas"
  }
  return(
    c(
      "#!/bin/bash",
      "#",
      paste("#SBATCH --time=", time,sep=""),
      partition_line,
      "#SBATCH --nodes=1",
      paste("#SBATCH -c",Ncpu),
      paste("#SBATCH --mem=",mem_size,sep=""),
      paste("#SBATCH --error",err),
      paste("#SBATCH --out",log),
      "",
      "module load biology",
      paste("module load",plink_pkg)
    )
  )
}

print_sh_file<-function(path,prefix,cmd){
  cmd = c(prefix,cmd)
  write.table(file=path,t(t(cmd)),row.names = F,col.names = F,quote = F)
}


disc_data_using_cut<-function(x,cuts=5,min_bin_size=10,useQuantiles=T){
  y = NULL
  if(!is.numeric(x)){y = factor(x)}
  if(length(unique(x))<=cuts){y = factor(x)}
  if(is.null(y)&&!useQuantiles){y=factor(cut(x,breaks=cuts, ordered_result=T))}
  if(is.null(y)&&useQuantiles){y=factor(cut(x,breaks=quantile(x,probs=seq(from=0,to=1,length.out = cuts+1)), ordered_result=T))}
  table_y = table(y)
  while(any(table_y<min_bin_size) && !all(y==y[1],na.rm=T)){
    curr_levels = levels(y)
    j = which(table_y==min(table_y))[1]
    ll = names(j)
    j2 = j-1
    if(j==1){j2=2}
    ll2 = names(table_y)[j2]
    new_levels=curr_levels
    new_levels[j]=paste(ll2,ll,sep=",")
    new_levels[j2]=paste(ll2,ll,sep=",")
    levels(y) = new_levels
    table_y = table(y)
  }
  return(y)
}
transform_to_column_name<-function(x,data){
  if(is.numeric(x) && length(x) == 1){
    x = names(data)[x];return(x)
  }
  z = x
  if(is.numeric(z) && length(z)>0){
    zz = c()
    for(ii in z){zz = c(zz,names(data)[ii])}
    z = zz
  }
  return(z)
}
run_discrete_ci_test<-function(x,y,z,data,test="x2",...){
  x = transform_to_column_name(x,data)
  y = transform_to_column_name(y,data)
  z = transform_to_column_name(z,data)
  if(length(z)==0){
    inds = !apply(is.na(data[,c(x,y)]),1,any)
    testp = ci.test(x,y,data=data[inds,c(x,y)],test=test,...)$p.value
    gc()
    return(testp)
  }
  inds = !apply(is.na(data[,c(x,y,z)]),1,any)
  testp = ci.test(x,y,z,data[inds,c(x,y,z)],test=test,...)$p.value
  gc()
  return(testp)
}
run_discrete_ci_test_fast<-function(x,y,z,data,test="mc-mi",...){
  return(run_discrete_ci_test(x,y,z,data,test="x2-adf"))
}
run_ci_test_one_is_numeric<-function(x,y,z,data,...){
  x = transform_to_column_name(x,data)
  y = transform_to_column_name(y,data)
  z = transform_to_column_name(z,data)
  yv = data[,y];xv = data[,x];zzv = NULL
  if(length(z)>0){zzv = data[,z]}
  if(is.numeric(yv) && is.numeric(xv)){
    if(length(z)==0){
      inds = !is.na(xv) & !is.na(yv)
      data1 = data[inds,c(x,y)]
      data1[[1]] = as.numeric(data1[[1]])
      data1[[2]] = as.numeric(data1[[2]])
      yv=NULL;xv=NULL;gc()
      return(ci.test(x,y,data=data1,test="cor")$p.value)
    }
    else{
      d1 = data.frame(x=xv,zzv);d2=data.frame(y=yv,zzv)
      lm1  = lm(x~.,data=d1)$residuals
      lm2  = lm(y~.,data=d2)$residuals
      inds = intersect(names(lm1),names(lm2))
      yv=NULL;xv=NULL;zzv=NULL;d1=NULL;d2=NULL;gc()
      return(cor.test(lm1[inds],lm2[inds])$p.value)
    }
  }
  if(is.numeric(yv)&&length(z)>0){
    summ = summary(lm(y~.,data=data.frame(y=yv,x=xv,zzv)))
    yv=NULL;xv=NULL;gc()
    lmres = summ$coefficients[2,4]
    summ=NULL;gc()
    return(lmres)
  }
  if(is.numeric(xv)&&length(z)>0){
    summ = summary(lm(x~.,data=data.frame(x=xv,y=yv,zzv)))
    yv=NULL;xv=NULL;zzv=NULL;gc()
    lmres = summ$coefficients[2,4]
    summ=NULL;gc()
    return(lmres)
  }
  if(is.numeric(yv)&&length(z)==0){
    summ = summary(lm(y~.,data=data.frame(y=yv,x=xv)))
    yv=NULL;xv=NULL;gc()
    lmres = summ$coefficients[2,4]
    summ=NULL;gc()
    return(lmres)
  }
  if(is.numeric(xv)&&length(z)==0){
    summ = summary(lm(x~.,data=data.frame(x=xv,y=yv)))
    yv=NULL;xv=NULL;zzv=NULL;gc()
    lmres = summ$coefficients[2,4]
    summ=NULL;gc()
    return(lmres)
  }
  yv=NULL;xv=NULL;gc()
  return(run_discrete_ci_test(x,y,z,data,...))
}

# assumes that x is binary and y is binary or linear
run_ci_logistic_test<-function(x,y,z,data,usespeedglm=T,...){
  #print("in logistic")
  x = transform_to_column_name(x,data)
  y = transform_to_column_name(y,data)
  z = transform_to_column_name(z,data)
  yv = data[,y];xv = data[,x];zzv = NULL
  glm_func = glm
  if(usespeedglm){
    glm_func=speedglm
  }
  if(length(z)>0){
    zzv = data[,z]
    summ = NULL
    try({summ = summary(glm_func(factor(x)~.,family=binomial(link="logit"),data=data.frame(x=xv,y=yv,zzv)))})
    if(is.null(summ)){try({summ = summary(glm(factor(x)~.,
                                              family=binomial(link="logit"),data=data.frame(x=xv,y=yv,zzv)))})}
    if(is.null(summ)){return(NULL)}
    return(as.numeric(as.character(summ$coefficients[2,4])))
  }
  summ = NULL
  try({summ = summary(glm_func(factor(x)~.,family=binomial(link="logit"),data=data.frame(x=xv,y=yv)))})
  if(is.null(summ)){
    try({summ = summary(glm(factor(x)~.,family=binomial(link="logit"),data=data.frame(x=xv,y=yv)))})}
  if(is.null(summ)){return(NULL)}
  yv=NULL;xv=NULL;gc()
  return(as.numeric(as.character(summ$coefficients[2,4])))
}

# order in the name implies order of try
# if x or y are binary and others are binary or linear use logistic
# then, if x or y are linear use linear regression
# finally, both are discrete and at least one is not binary: use discrete tests
run_ci_test_logistic_linear_discrete<-function(x,y,z,data,...){
  x = transform_to_column_name(x,data)
  y = transform_to_column_name(y,data)
  z = transform_to_column_name(z,data)
  yv = data[,y];xv = data[,x]
  N_x = length(unique(xv)); N_y = length(unique(yv))
  #print(c(N_x,N_y))
  # check logistic
  if(N_x==2){
    return(run_ci_logistic_test(x,y,z,data,...))
  }
  if(N_y==2){
    return(run_ci_logistic_test(y,x,z,data,...))
  }
  gc()
  return(run_ci_test_one_is_numeric(x,y,z,data,...))
}

#run_ci_test_one_is_numeric(2,3,NULL,data=DATA)
# data is a list with DATA and discDATA
run_ci_test<-function(x,y,z,data,test="mi-adf",...){
  yv = data$DATA[,y];xv = data$DATA[,x]
  if(is.numeric(yv)||is.numeric(xv)){
    return(run_ci_test_one_is_numeric(x,y,z,data$DATA))
  }
  return(run_discrete_ci_test(x,y,z,data$discDATA,test=test))
}
get_CI_pairwise_network<-function(Xs,Z,data,func=run_discrete_ci_test,...){
  n = length(Xs)
  m = matrix(1,nrow=n,ncol=n)
  colnames(m) = Xs; rownames(m)=Xs
  for(i in 2:n){
    for(j in 1:(i-1)){
      m[i,j] = func(Xs[i],Xs[j],Z,data=data,...)
      m[j,i] = m[i,j]
    }
  }
  return(m)
}