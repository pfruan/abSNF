diffus_vec=function(signals,snet,type,beta=0.75,iter=10,difference=1e-6,top=100){

  require(Rfast)

  #pre-precoss association signals
  if (type=="pvalue"){
    signals[,2]=-qnorm(signals[,2])
    colnames(signals)=c("gene","score")
  }else{
    colnames(signals)=c("gene","score")
  }

  #column nomorlize the network
  snet=t(t(snet)/(colsums(snet)+.Machine$double.eps))

  #signals: gene  score
  #intialize
  p=as.matrix(signals$score)
  p1=p
  j=1

  #diffusion on network
  repeat{

    p=p1
    p1<-beta*mat.mult(snet,p)+(1-beta)*(signals$score)
    p_diff=sum(abs((p1-p)))
    j=j+1

    if (j>iter){break}
    if (p_diff<1e-6) {break}

  }

  #prioritize genes by their enhanced score
  res=data.frame(gene=signals$gene,score=p1)
  res=res[order(res$score,decreasing = T),]
  res=res[1:top,]
  return(res)
}
