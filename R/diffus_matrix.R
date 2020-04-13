diffus_matrix=function(s0,adjacency,alpha=0.75,iter=10,difference=1e-6){

  require(Rfast)

  #prepocessing
  gene=intersect(rownames(s0),rownames(adjacency))

  s0=s0[rownames(s0) %in% gene,]
  s0=s0[,colnames(s0) %in% gene]
  s0=s0[,order(colnames(s0))]
  s0=s0[order(rownames(s0)),]

  adjacency=adjacency[rownames(adjacency) %in% gene,]
  adjacency=adjacency[,colnames(adjacency) %in% gene]
  adjacency=adjacency[,order(colnames(adjacency))]
  adjacency=adjacency[order(rownames(adjacency)),]

  diag(adjacency)=0
  adjacency=t(t(adjacency)/colsums(adjacency))

  #initialize
  snet_1=s0
  snet=snet_1

  #diffusion on adjacency matrix
  for(kk in 1:iter){

    snet_1<-alpha*mat.mult(adjacency,snet)+(1-alpha)*(s0)
    diff=max(abs(snet_1-snet))
    print(c("iteration:",kk,"difference:",diff))
    if(diff<difference){return(snet_1)}
    snet=snet_1

  }

  return(snet_1)
}
