# Set the session
required_libs = c("igraph","bnlearn","MRPRESSO","mixtools","mclust",
                  "optparse","limma","MendelianRandomization")
lib_loc = "~/R/packages3.5"
lib_loc = c(lib_loc,.libPaths())
for (lib_name in required_libs){
  tryCatch({library(lib_name,character.only = T,lib.loc = lib_loc)},
           error = function(e) {
             print(paste("Cannot load",lib_name,", please install"))
           })
}
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
source("~/Desktop/repos/cGAUGE/R/cGAUGE.R")
N = 2000
minBeta = 0.1
maxBeta = 0.9
minMAF = 0.05
maxMAF = 0.4
p=3
p1 = 1e-4
p2 = 0.001

all_mr_res = list()
for(rep in 1:100){
  print(rep)
  # We simulate from a mixture of real instruments and confounder instruments
  # traits graph: U->X and U->Y
  # instruments of X: G->X
  # instruments of U: G->U
  B = matrix(0,3,3)
  phenos = c("X","Y","U")
  colnames(B) = phenos
  rownames(B) = colnames(B)
  # edges are cloumns to rows
  B["X","U"] = (-1^rbinom(1,1,0.5))*runif(1,minBeta,maxBeta)
  B["Y","U"] = (-1^rbinom(1,1,0.5))*runif(1,minBeta,maxBeta)
  Bg = igraph::graph_from_adjacency_matrix(t(abs(B)>0))

  # Add instruments
  num_ivs1 = 10
  num_ivs2 = 20
  # instruments into X
  for(j in 1:num_ivs1){
    v = rep(0,nrow(B))
    v[1] = (-1^rbinom(1,1,0.5))*runif(1,minBeta,maxBeta)
    B = cbind(B,v)
    B = rbind(B,rep(0,ncol(B)))
  }
  # instruments into U
  for(j in 1:num_ivs2){
    v = rep(0,nrow(B))
    v[3] = (-1^rbinom(1,1,0.5))*runif(1,minBeta,maxBeta)
    B = cbind(B,v)
    B = rbind(B,rep(0,ncol(B)))
  }
  num_ivs = num_ivs1+num_ivs2
  ivs = paste("IV",1:num_ivs,sep="")
  colnames(B)[4:ncol(B)] = ivs
  rownames(B) = colnames(B)
  Bg = igraph::graph_from_adjacency_matrix(t(abs(B)>0))
  
  mafs = runif(num_ivs,minMAF,maxMAF)
  # parameters of the current SEM
  I = diag(ncol(B))
  A = solve(I-B)
  theoretical_cov = A%*%t(A)
  # Use the simulated B to generate data, use Fisher's fixed point method
  simulated_data = c()
  for(s in 1:N){
    v_t = rnorm(p)
    v_iv = rep(0,num_ivs)
    for(j in 1:num_ivs){v_iv[j] = rbinom(1,1,mafs[j])}
    e = c(v_t,v_iv)
    v = (A%*%e)[,1]
    simulated_data = rbind(simulated_data,v)
  }
  colnames(simulated_data) = colnames(B)
  
  df = data.frame(simulated_data)
  GWAS_Ps = matrix(1,num_ivs,p,dimnames = list(ivs,phenos))
  GWAS_effects = matrix(0,num_ivs,p,dimnames = list(ivs,phenos))
  GWAS_ses = matrix(0,num_ivs,p,dimnames = list(ivs,phenos))
  for(pheno in phenos){
    gwas_res = sapply(ivs,run_lm,x=pheno,z=NULL,df = df)
    GWAS_Ps[,pheno] = gwas_res[4,]
    GWAS_effects[,pheno] = gwas_res[1,]
    GWAS_ses[,pheno] = gwas_res[2,]
  }
  
  G_it = GWAS_Ps < p1
  # Run MR
  # pleio size is set to 100 - no filtering of variants (a standard MR analysis)
  mr_anal_res = list(
    "Egger" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
                                       pleio_size=100,pruned_lists=NULL,func=mr_egger,robust=T),
    "IVW" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
                                     pleio_size=100,pruned_lists=NULL,func=mr_ivw,robust=T)
  )
  # Add MRPRESSO
  mrpresso_res = c()
  cgauge_mrpresso_res = c()
  try({
    for(tr1 in phenos){
      currivs = rownames(GWAS_Ps)[G_it[,tr1]]
      for(tr2 in phenos){
        if(tr1==tr2){next}
        X = data.frame(E1b=GWAS_effects[currivs,tr1],O1b=GWAS_effects[currivs,tr2],
                       E1sd=GWAS_ses[currivs,tr1],O1sd=GWAS_ses[currivs,tr2])
        try({
          res = mr_presso(BetaOutcome = "O1b", BetaExposure = "E1b", 
                          SdOutcome = "O1sd", SdExposure = "E1sd",data=X,
                          OUTLIERtest=T,
                          DISTORTIONtest = T,
                          NbDistribution = 1000,SignifThreshold = 0.1)
          if(is.na(res$`Main MR results`[2,"P-value"])){
            v = c(tr1,tr2,unlist(res$`Main MR results`[1,]))
          }
          else{
            v = c(tr1,tr2,unlist(res$`Main MR results`[2,]))
          }
          v["GlobalTestP"] = res$`MR-PRESSO results`$`Global Test`$Pvalue
          mrpresso_res = rbind(mrpresso_res,v)
        })
      }
    }
    if(!is.null(dim(mrpresso_res))){
      mrpresso_res = as.data.frame(mrpresso_res)
      for(j in 3:ncol(mrpresso_res)){
        mrpresso_res[[j]] = as.numeric(as.character(mrpresso_res[[j]]))
      }  
    }
  })
  
  # Add cGAUGE estimates - 
  # we treat U as unobserved, so there is no trait skeleton involved
  cor.test(df$X,df$Y)$p.value
  phenos = c("X","Y")
  trait_pair_pvals = list()
  for(pheno1 in phenos){
    trait_pair_pvals[[pheno1]] = list()
    for(pheno2 in phenos){
      if(pheno1==pheno2){next}
      gwas_res = t(sapply(ivs,run_lm,x=pheno1,z=pheno2,df = df))
      gwas_res = gwas_res[,4:1]
      trait_pair_pvals[[pheno1]][[pheno2]] = gwas_res
    }
  }
  G_vt = extract_skeleton_G_VT(GWAS_Ps,trait_pair_pvals,P1=p1,P2=p2,test_columns = 1)[[1]]
  G_t = matrix(T,2,2,dimnames = list(phenos,phenos))
  iv_sets = list()
  uniquely_mapped_ivs = rownames(G_vt)[rowSums(G_vt)==1]
  for(tr1 in phenos){
    iv_sets[[tr1]] = list()
    for(tr2 in phenos){
      iv_sets[[tr1]][[tr2]] = rownames(GWAS_Ps)[GWAS_Ps[,tr1]<p1]
      iv_sets[[tr1]][[tr2]] = intersect(iv_sets[[tr1]][[tr2]],
                                        uniquely_mapped_ivs)
    }
  }
  cgauge_mr_anal_res = list(
    "Egger" = run_pairwise_mr_analyses_with_iv_sets(GWAS_effects,GWAS_ses,iv_sets,
                                                    func=mr_egger,robust=T),
    "IVW" = run_pairwise_mr_analyses_with_iv_sets(GWAS_effects,GWAS_ses,iv_sets,
                                                  func=mr_ivw,robust=T)
  )
  
  try({
    # Add MRPRESSO
    for(tr1 in phenos){
      for(tr2 in phenos){
        if(tr1==tr2){next}
        currivs = iv_sets[[tr1]][[tr2]]
        if(length(currivs)<5){next}
        X = data.frame(E1b=GWAS_effects[currivs,tr1],O1b=GWAS_effects[currivs,tr2],
                       E1sd=GWAS_ses[currivs,tr1],O1sd=GWAS_ses[currivs,tr2])
        try({
          res = mr_presso(BetaOutcome = "O1b", BetaExposure = "E1b", 
                          SdOutcome = "O1sd", SdExposure = "E1sd",data=X,
                          OUTLIERtest=T,
                          DISTORTIONtest = T,
                          NbDistribution = 100,SignifThreshold = 0.1)
          if(is.na(res$`Main MR results`[2,"P-value"])){
            cgauge_mrpresso_res = rbind(cgauge_mrpresso_res,
                                        c(tr1,tr2,unlist(res$`Main MR results`[1,])))
          }
          else{
            cgauge_mrpresso_res = rbind(cgauge_mrpresso_res,
                                        c(tr1,tr2,unlist(res$`Main MR results`[2,])))
          }
        })
      }
    }
    # if(!is.null(dim(cgauge_mrpresso_res))){
    #   cgauge_mrpresso_res = as.data.frame(cgauge_mrpresso_res)
    #   for(j in 3:ncol(cgauge_mrpresso_res)){
    #     cgauge_mrpresso_res[[j]] = as.numeric(as.character(cgauge_mrpresso_res[[j]]))
    #   }
    # }
  })
  
  # Add the current results to our data matrix
  mr_anal_res$Egger$rep = rep
  mr_anal_res$IVW$rep = rep
  mrpresso_res$rep = rep
  cgauge_mr_anal_res$Egger$rep = rep
  cgauge_mr_anal_res$IVW$rep = rep
  cgauge_mrpresso_res$rep = rep
  all_mr_res[["egger"]] = rbind(all_mr_res[["egger"]],mr_anal_res$Egger)
  all_mr_res[["ivw"]] = rbind(all_mr_res[["ivw"]],mr_anal_res$IVW)
  all_mr_res[["mrpresso"]] = rbind(all_mr_res[["mrpresso"]],mrpresso_res)
  try({all_mr_res[["c_egger"]] = rbind(all_mr_res[["c_egger"]],cgauge_mr_anal_res$Egger)})
  try({all_mr_res[["c_ivw"]] = rbind(all_mr_res[["c_ivw"]],cgauge_mr_anal_res$IVW)})
  try({all_mr_res[["c_mrpresso"]] = rbind(all_mr_res[["c_mrpresso"]],cgauge_mrpresso_res)})
}

# Put the -log10 p-value results in a single data frame
eggerXY = all_mr_res[["egger"]][all_mr_res[["egger"]][,1]=="X" &
                                  all_mr_res[["egger"]][,2] == "Y",]
eggerXY$name = "Egger X->Y"
eggerYX = all_mr_res[["egger"]][all_mr_res[["egger"]][,1]=="Y" &
                                  all_mr_res[["egger"]][,2] == "X",]
eggerYX$name = "Egger Y->X"

ivwXY = all_mr_res[["ivw"]][all_mr_res[["ivw"]][,1]=="X" &
                                  all_mr_res[["ivw"]][,2] == "Y",]
ivwXY$name = "IVW X->Y"
ivwYX = all_mr_res[["ivw"]][all_mr_res[["ivw"]][,1]=="Y" &
                                  all_mr_res[["ivw"]][,2] == "X",]
ivwYX$name = "IVW Y->X"

pressoXY = all_mr_res[["mrpresso"]][all_mr_res[["mrpresso"]][,1]=="X" &
                              all_mr_res[["mrpresso"]][,2] == "Y",]
pressoXY$name = "MRPRESSO X->Y"
pressoYX = all_mr_res[["mrpresso"]][all_mr_res[["mrpresso"]][,1]=="Y" &
                              all_mr_res[["mrpresso"]][,2] == "X",]
pressoYX$name = "MRPRESSO Y->X"
names(pressoYX)[grepl("P-value",names(pressoYX))] = "p"
names(pressoXY)[grepl("P-value",names(pressoXY))] = "p"

mr_res_p_df = rbind(
  eggerXY[,c("name","p")],eggerYX[,c("name","p")],
  ivwXY[,c("name","p")],ivwYX[,c("name","p")],
  pressoXY[,c("name","p")],pressoYX[,c("name","p")]
)
mr_res_p_df$log10p = -log10(mr_res_p_df$p)
mr_res_p_df$log10p = pmin(mr_res_p_df$log10p,10,na.rm=T)

library(ggplot2)
ggplot(mr_res_p_df, aes(x=name, y=log10p)) + 
  geom_boxplot(notch=TRUE)+ coord_flip() + 
  ylab("-log10 P-value") + xlab("") + 
  geom_hline(yintercept = 2, linetype="dotted",color = "red", size=1.5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size=14)
  )

save.image(file = "~/Desktop/causal_inference_projects/ms3/fig2a_session.R")

