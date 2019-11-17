# Functions for simulating data and obtaining the cGAUGR/MR results 
# Get the effect size, se, and p-value for x~y|z 
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
# get distances (ancestry graph)
igraph_directed_distances<-function(Bg){
  distances = igraph::distances(Bg) # undirected distances
  for(i in 1:nrow(distances)){
    for(j in 1:ncol(distances)){
      distances[i,j] = length(igraph::shortest_paths(Bg,from=j,to=i)$vpath[[1]])-1
    }
  }
  return(distances)
}
# a function to add the known distances to the results
add_distances<-function(m,dists,newcolname = "KnownDistance"){
  m[,1] = as.character(m[,1])
  m[,2] = as.character(m[,2])
  v = c()
  for(i in 1:nrow(m)){
    v = c(v,dists[m[i,2],m[i,1]]) # add dist from exposure to outcome
  }
  m = cbind(m,v)
  colnames(m)[ncol(m)] = newcolname
  return(m)
}

