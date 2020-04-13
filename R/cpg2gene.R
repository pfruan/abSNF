cpg2gene =  function(cpg, method='smallest'){

  #cpg: gene, cpg, pvalue

  if(method=='smallest')
  {
    tapply(cpg$pvalue, cpg$gene, min) -> gene2weight
    data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
    gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
  }

  if(method=='fisher')
  {
    tapply(cpg$pvalue, cpg$gene, function(x){sum(-2*log(x))->chi2; p=1-pchisq(chi2, 2*length(x)); return(p);}) -> gene2weight
    data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
    gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
  }

  if(method=='bon')
  {
    tapply(cpg$pvalue, cpg$gene, function(x){tmp2 = p.adjust(x, method='bonferroni'); return(min(tmp2));}) -> gene2weight
    data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
    gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
  }

  return(gene2weight)
}

