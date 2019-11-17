library("igraph")
library("MRPRESSO")

N = 2000
minBeta = 0.1
maxBeta = 0.9
minMAF = 0.05
maxMAF = 0.4
p=3
p1 = 1e-4

# Case 1: a mixture of real instruments and confounder instruments
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
plot(Bg)
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
plot(Bg)

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
  print(pheno)
  gwas_res = sapply(ivs,run_lm,x=pheno,z=NULL,df = df)
  GWAS_Ps[,pheno] = gwas_res[4,]
  GWAS_effects[,pheno] = gwas_res[1,]
  GWAS_ses[,pheno] = gwas_res[2,]
}

G_it = GWAS_Ps < p1
print(colSums(G_it))
# Run MR
# pleio size is set to 100 - no filtering of variants (a standard MR analysis)
mr_anal_res = list(
  "Egger" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
              pleio_size=100,pruned_lists=NULL,func=mr_egger,robust=T),
  "IVW" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
              pleio_size=100,pruned_lists=NULL,func=mr_ivw,robust=T)
)
# Add MRPRESSO
print("Done with Egger and IVW, running MRPRESSO")
mrpresso_res = c()
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
                        NbDistribution = 100,SignifThreshold = 0.1)
        if(is.na(res$`Main MR results`[2,"P-value"])){
          mrpresso_res = rbind(mrpresso_res,
                               c(tr1,tr2,unlist(res$`Main MR results`[1,])))
        }
        else{
          mrpresso_res = rbind(mrpresso_res,
                               c(tr1,tr2,unlist(res$`Main MR results`[2,])))
        }
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




