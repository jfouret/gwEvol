LRT_tests=function(GENE){
  tests=c('b_free vs M0','bsA vs M1','bsA vs bsA1')
  TEST=as.data.frame(tests)
  TEST$tests=factor(TEST$tests,levels = tests)
  TEST$ddl=c(1,2,1)
  TEST$Alt=c(GENE$lnL[GENE$model=='b_free'],GENE$lnL[GENE$model=='bsA'],GENE$lnL[GENE$model=='bsA'])
  TEST$Null=c(GENE$lnL[GENE$model=='M0'],GENE$lnL[GENE$model=='M1'],GENE$lnL[GENE$model=='bsA1'])
  TEST$chi2=2*(TEST$Alt-TEST$Null)
  TEST$pval=1-pchisq(TEST$chi2,TEST$ddl)
  TEST$gene_name=GENE[1,"gene_name"]
  TEST$target=GENE[1,"target"]
  #TEST$gene_name=gene
  #TEST$target=specie
  return(TEST)
}

evol_tests=function(data,padj_meth="fdr"){
  genes=unique(data$gene_name)
  targets=unique(data$target)
  targets=targets[!is.na(targets)]
  init=0
  require(dplyr)
  allTEST=by(data,data[,c('gene_name','target')],LRT_tests) %>% bind_rows()
  tests=c('b_free vs M0','bsA vs M1','bsA vs bsA1')
  for (test_name in tests){
      allTEST[allTEST$tests==test_name,"padj"]=p.adjust(allTEST[allTEST$tests==test_name,"pval"],method = padj_meth)
  }
  return(allTEST)
}
