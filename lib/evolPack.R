LRT_tests=function(gene,specie,data,tests){
  GENE=subset(data,(gene_name==gene) & ( (target==specie)|(is.na(target)) )  )
  TEST=as.data.frame(tests)
  TEST$tests=factor(TEST$tests,levels = tests)
  TEST$ddl=c(1,2,1)
  TEST$Alt=c(GENE$lnL[GENE$model=='b_free'],GENE$lnL[GENE$model=='bsA'],GENE$lnL[GENE$model=='bsA'])
  TEST$Null=c(GENE$lnL[GENE$model=='M0'],GENE$lnL[GENE$model=='M1'],GENE$lnL[GENE$model=='bsA1'])
  TEST$chi2=2*(TEST$Alt-TEST$Null)
  TEST$pval=1-pchisq(TEST$chi2,TEST$ddl)
  TEST$gene_name=gene
  TEST$target=specie
  return(TEST)
}

evol_tests=function(data,alpha=0.05,padj_meth="fdr"){#TODO ajdFilter
  tests=c('b_free vs M0','bsA vs M1','bsA vs bsA1')
  genes=unique(data$gene_name)
  targets=unique(data$target)
  targets=targets[!is.na(targets)]
  init=0
  for (gene in genes){
    for (target in targets){
      if (init==0){
        allTEST=LRT_tests(gene,target,data,tests)
        init=1
      }else{
	#LRT_tests(gene,target,data)
        allTEST=rbind(allTEST,LRT_tests(gene,target,data,tests))
      }
    }
  }
  allTEST=allTEST[,c('tests','gene_name','target','pval')]
  #allTEST$padj=p.adjust(allTEST$pval,method = padj_meth)

  ###TODO ==> ONLY ONE PADJ OR MANY ???
  allTEST$padj=p.adjust(allTEST$pval,method = padj_meth)
  pvalTEST=spread(allTEST,tests,pval)
  pvalTEST$padj=NULL
  pvalTEST=getMethodsResults(pvalTEST,tests,alpha)
  padjTEST=spread(allTEST,tests,padj)
  padjTEST$pval=NULL
  padjTEST=getMethodsResults(padjTEST,tests,alpha)
  #init=0
  #for (test in tests){
    #New_name=paste(test,'adj',sep=' ')
    #if (init==0){
      #tests_adj=c(New_name)
      #init=1
    #}else{
      #tests_adj=c(tests_adj,New_name)
    #}

    #if (adj==T){
      #allTEST[New_name]=p.adjust(allTEST[test][,1],method = padj_meth)
    #}else{
      #allTEST[New_name]=allTEST[test][,1]
    #}
  #}
  return(list(allTEST,pvalTEST,padjTEST))
#TODO filter and use padj ! ! !
}


getMethodsResults=function(tests_df,tests,alpha){
  bfree_M0=subset(tests_df,!is.na(tests_df[,'b_free vs M0']))
  bsA_M1=subset(tests_df,!is.na(tests_df[,'bsA vs M1']))
  bsA_bsA1=subset(tests_df,!is.na(tests_df[,'bsA vs bsA1']))
  new_df=bfree_M0[,c('gene_name','target')]

  bfree_M0$pM0=FALSE
  #pM0 T si w1(b_free)>w0(M0)
  for (i in (1:dim(bfree_M0)[1])){
    sdata=subset(data,(gene_name==bfree_M0$gene_name[i])&((target==bfree_M0$target[i])|(is.na(target)))&((model=='M0')|(model=='b_free')))
    w0=subset(sdata,model=='M0')$w0
    w1=subset(sdata,model=='b_free')$w1
    if (w1>w0){
      bfree_M0$pM0[i]=T
    }
  }
#tests=c('b_free vs M0','bsA vs M1','bsA vs bsA1')
  new_df$Branch=F
  new_df$BranchSites=F


  for (i in 1:dim(new_df)[1]) {
    if (    ( bfree_M0[i,'b_free vs M0']<alpha )  &  ( bfree_M0$pM0[i]==T )  ){
      new_df$Branch[i]=T
    }
    if (   (bsA_M1[i,'bsA vs M1']<alpha)    & (bsA_bsA1[i,'bsA vs bsA1']<alpha) ) {
      new_df$BranchSites[i]=T
    }
  }
  return(new_df)
}

makePhyloGraph=function(data,gene,specie){
  require(ggthemes)
  require(ggplot2)
  require(gridExtra)
  require(cowplot)
  require(gplots)
  require(tidyr)
  data=unique(data)
  data$fmodel=factor(data$model,levels = c('b_free','M0','bsA','bsA1','M1'))
  GENE=subset(data,(gene_name==gene) & ( (target==specie)|(is.na(target)) )  )
  
  GENE$p2=1-GENE$p0-GENE$p1
  
  wGENE=gather(GENE,wclass,w,w0:w2)
  wGENE$wclass=factor(wGENE$wclass,levels=c('w0','w1','w2'))
  wGENE=subset(wGENE,!is.na(w))
  wGENE$wlim=wGENE$w
  wGENE$wlim[wGENE$w>2]=2
  
  wGENE$wsci=format(wGENE$w,scientific = T,digits = 2)
  wGENE$wsci[wGENE$wsci=="1.0e+00"]="1"
  
  omega=ggplot(wGENE,aes(x=wclass,y=wlim,fill=wlim)) + 
    geom_bar(stat='identity') +
    geom_text(stat='identity',aes(label=wsci),vjust=-0.4,fontface='bold')+
    theme_hc() + # a theme based on Highcharts JS.
    scale_fill_gradient(low = "grey60", high = "black") +
    guides(fill=FALSE) +
    scale_x_discrete(labels=c('w0' = expression(omega[0]),
                              'w1' = expression(omega[1]),
                              'w2' = expression(omega[2])))+
    labs(x='Categories',y='dN/dS ratio')+
    facet_grid(.~fmodel,space = "free",scales = "free")+
    theme(strip.text.y = element_text(size=14, angle=0))
  
  
  
  pGENE=gather(GENE,pclass,p,c(p0,p1,p2))
  pGENE$pclass=factor(pGENE$pclass,levels=c('p0','p1','p2'))
  pGENE=subset(pGENE,!is.na(p))
  
  proportion=ggplot(pGENE,aes(x='',y=100*p,fill=pclass))+
    geom_bar(width = 1,stat='identity')+
    coord_polar(theta = "y", start=0)+
    #geom_text(stat='identity',aes(label=p),vjust=-0.4,fontface='bold')+
    theme_hc() +
    scale_fill_grey(labels=c('p0' = expression("p"[0]),
                             'p1' = expression("p"[1]),
                             'p2' = expression("p"[2])),
                    name="Proportions:") +
    #guides(fill=FALSE) +
    scale_x_discrete()+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          strip.text.y = element_text(size=16, angle=0))+
    facet_grid(fmodel~.)
  
  tests=c('b_free vs M0','bsA vs M1','bsA vs bsA1')
  
  
  TEST=as.data.frame(tests)
  TEST$tests=factor(TEST$tests,levels = tests)
  TEST$ddl=c(1,2,1)
  TEST$Alt=c(GENE$lnL[GENE$model=='b_free'],GENE$lnL[GENE$model=='bsA'],GENE$lnL[GENE$model=='bsA'])
  TEST$Null=c(GENE$lnL[GENE$model=='M0'],GENE$lnL[GENE$model=='M1'],GENE$lnL[GENE$model=='bsA1'])
  
  
  TEST$chi2=2*(TEST$Alt-TEST$Null)
  
  TEST$pval=1-pchisq(TEST$chi2,TEST$ddl)
  lTEST=gather(TEST,model,lnL,c(Null,Alt))
  lTEST$model=factor(lTEST$model,levels = c('Null','Alt'))
  max=max(lTEST$lnL)
  min=min(lTEST$lnL)
  
  diff=(max-min)/0.9
  min_plot=min
  max_plot=min+diff
  
  lnL=ggplot(lTEST,aes(x=model,y=lnL,fill=lnL))+
    coord_flip(ylim=c(min_plot,max_plot))+
    geom_bar(stat='identity')+
    scale_fill_gradient(low = "black", high = "grey60") +
    theme_hc() +
    guides(fill=FALSE) +
    labs(y="logLikelihood",x="Models")+
    facet_grid(tests~.)+
    theme(strip.text.y = element_text(size=12, angle=320))
  
  TEST$pvalsci=format(TEST$pval,scientific = T,digits = 2)
  
  max_pval=max(c(TEST$pval,0.05))
  
  pval=ggplot(TEST,aes(x=factor(tests,levels=c('bsA vs bsA1','bsA vs M1','b_free vs M0')),y=pval,fill=pval))+
    coord_flip()+
    geom_bar(stat='identity')+
    geom_hline(yintercept=0.05,colour='red',size=1.1)+
    geom_label(stat='identity',aes(label=pvalsci,y=max_pval*(1/10),fontface='bold'),fill="white")+
    scale_fill_gradient(low = "black", high = "grey60") +
    guides(fill=FALSE) +
    labs(y="p-value",x="")+
    theme_hc() +
    theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  
  
  bottom=plot_grid(lnL,pval,nrow=1,rel_widths = c(0.55,0.45))
  left=plot_grid(omega,bottom,ncol=1,rel_heights = c(0.6,0.4))
  #labels=c('Estimation of dN/dS ratio for different evolutionary models','Likelihood ratio tests')
  all=plot_grid(left,proportion,nrow=1,rel_widths = c(0.7,0.3))
  title <- ggdraw() + draw_label(paste("Evolutionary analysis of ",gene," with ",specie," as forground branch",sep=''), fontface='bold')
  plot_grid(title, all, ncol=1, rel_heights=c(0.03, 1))
  
  return(plot_grid(title, all, ncol=1, rel_heights=c(0.03, 1))) # rel_heights values control title margins
}
