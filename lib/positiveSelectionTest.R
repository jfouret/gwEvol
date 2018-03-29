#!/usr/bin/Rscript
author='Julien FOURET'
contact='julien.fouret@fouret.me'
year='2016'
version='SEDMATCHGITVERSION'
gitRepository='SEDMATCHGITREPO'
startMessage=paste('Perform Log-likelihood ratio tests for positive selection analysis '," ; Version : ",version," ; ",year," ; Author : ",author," for more informations or enquiries please contact ",contact," (WARNING: pval and FC column names must begin respectively with 'p.val' and 'FC')",sep='')
suppressMessages(library("argparse"))
parser <- ArgumentParser(description=startMessage)
parser$add_argument("-outDir",
                    type="character",
                    metavar="/path",
                    help="output directory for the whole analysis")
parser$add_argument("-alpha",
                    type="character",
                    default='0.05',
                    help="1st type error cut-off for positively selected genes")
args <- parser$parse_args()

# LOADING LIBRARIES
suppressMessages(library(ggthemes))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(gplots))
suppressMessages(library(tidyr))
suppressMessages(library(R.utils))
suppressMessages(library(dplyr))
suppressMessages(library(VennDiagram))

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

#DEFINE PARAMETERS
outDir=getAbsolutePath(args$outDir)
resDir=paste(outDir,'results',sep='/')
statDir=resDir
setwd(statDir)
parametersFileName=paste(resDir,'parameters.tab',sep='/')
testResultsFileName=paste(statDir,'testsResults_',sep='/')

#IMPORT DATAFRAME
param=read.table(file = parametersFileName,header = T,sep = "\t",na.strings = "NA",dec = "." )
print("Test computation")
#put in parameters adjustments #TODO
data=evol_tests(param)
print("Test results writing")
write.table(data,file=paste(testResultsFileName,'raw','.tab',sep=''),quote=F,row.names=F,sep="\t")

alpha=as.numeric(args$alpha)

print("Summary graphs")

for (chosen_stat in c("pval","padj")){
  sdata=data[,c("tests","gene_name","target",chosen_stat)]
  names(sdata)[4]='stat'
  sdata=spread(sdata,key=tests,value = stat)
  sdata$branch_ratio=mutate(subset(param,model=="b_free"),b_free_ratio=w1/w0)$b_free_ratio
  branch=subset(sdata,(`b_free vs M0`<alpha)&(`branch_ratio`>1))$gene_name
  branch_site=subset(sdata,(`bsA vs bsA1`<alpha)&(`bsA vs M1`<alpha))$gene_name
  vpval=venn.diagram(list("Branch"=branch,"Branch-site"=branch_site),height = 8, width = 8,fill = c("red", "green"),
               alpha = c(0.5, 0.5), cat.pos=c(0,45),cex = 2.5,cat.cex=2,cat.fontface = 4,lty =2, fontfamily =3, 
               filename = paste(chosen_stat,"_venn.svg",sep=''),imagetype='svg',main=paste("<",args$alpha," ",chosen_stat,sep=''),main.cex=3);
  write(paste(branch,collapse = "\n"),paste("branch_",chosen_stat,".txt",sep=''))
  write(paste(branch_site,collapse = "\n"),paste("branch_site_",chosen_stat,".txt",sep=''))

  # Unite

  bbsTests=subset(data,gene_name %in% intersect(branch,branch_site))
  bonlyTests=subset(data,gene_name %in% setdiff(branch,branch_site))
  bsonlyTests=subset(data,gene_name %in% setdiff(branch_site,branch))

  bbsTests$type="Both"
  bonlyTests$type="Branch only"
  bsonlyTests$type="Branch-site only"
  all=unique(rbind(bbsTests,bonlyTests,bsonlyTests)[,c("type","gene_name")])
  all$type=factor(all$type,levels =c("Branch only","Both","Branch-site only"),ordered = T)

  rownames(all)=all$gene_name
  rownames(param)=paste(param$gene_name,param$model)

  all[,"bsa_w2"]=param[paste(all$gene_name,"bsA"),"w2"]
  all[,"bfree_w0"]=param[paste(all$gene_name,"b_free"),"w0"]
  all[,"bfree_w1"]=param[paste(all$gene_name,"b_free"),"w1"]
  all[,"bsa_p1"]=param[paste(all$gene_name,"bsA"),"p1"]
  all[,"bsa_p2"]=param[paste(all$gene_name,"bsA"),"p2a"]+param[paste(all$gene_name,"bsA"),"p2b"]
  all[,"bsa_p1"]=param[paste(all$gene_name,"bsA"),"p1"]

  A=length(setdiff(branch_site,branch))
  B=length(setdiff(branch,branch_site))

  AB=length(intersect(branch,branch_site))

  base=seq(1,9)
  mbreaks=c(0.001*base,0.01*base,0.1*base,base,10*base,100*base,1000*base,10000*base,100000*base)

  venn=draw.pairwise.venn(A+AB,
                          B+AB,
                          AB, 
                          category = c("Branch-site", "Branch"), lty = rep("blank", 2),
                          fill = c("red", "blue"), 
                          alpha = rep(0.5, 2),
                          cat.pos = c(45,-45),
                          cat.dist = rep(0.1, 2))

  w2_w1=ggplot(all,aes(x=bfree_w1,y=bsa_w2,color=type))+ 
    #stat_density_2d(aes(fill = ..level..), geom = "polygon")+
    geom_point(size=0.8)+
    scale_y_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #geom_errorbar()+
    scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_fill_gradientn(colours =rev(RColorBrewer::brewer.pal(7,"RdYlBu")))+
    scale_color_manual(values =rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    theme_bw() + # a theme based on Highcharts JS.
    #labs(x='bfree_w1',y='w2(bsA)/w1(b_free)')+
    #facet_grid(.~type)+
    guides(color=F)+
    coord_cartesian(ylim=c(1,1000),xlim=c(0.01,10))+
    theme(strip.text.y = element_text(size=14, angle=0),
          panel.background = element_rect(fill='gray80'))

  w0_w1=ggplot(all,aes(x=bfree_w1,y=bfree_w0,color=type))+ 
    #stat_density_2d(aes(fill = ..level..), geom = "polygon")+
    geom_point(size=0.8)+
    guides(color=F)+
    scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #geom_errorbar()+
    scale_color_manual(values =rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    scale_fill_gradientn(colours =rev(RColorBrewer::brewer.pal(7,"RdYlBu")))+
    theme_bw() + # a theme based on Highcharts JS.
    #labs(x='bfree_w0',y='bfree_w0')+
    coord_cartesian(ylim=c(0,max(all$bfree_w0)),xlim=c(0.01,10))+
    #facet_grid(.~type)+
    theme(strip.text.y = element_text(size=14, angle=0),
          panel.background = element_rect(fill='gray80'))
  w2w1_p2=ggplot(all,aes(x=bsa_p2,y=bsa_w2/bfree_w1,color=type))+
    #stat_density_2d(aes(fill = ..level..), geom = "polygon")+
    geom_point(size=0.8)+
    guides(color=F)+
    scale_y_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_color_manual(values =rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    #geom_errorbar()+
    #scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_fill_gradientn(colours =rev(RColorBrewer::brewer.pal(7,"RdYlBu")))+
    theme_bw() + # a theme based on Highcharts JS.
    #labs(x='p2',y='bfree_w0')+
    coord_cartesian(ylim=c(1,25000),xlim=c(0,0.4))+
    #facet_grid(.~type)+
    theme(strip.text.y = element_text(size=14, angle=0),
          panel.background = element_rect(fill='gray80'))

  w0=ggplot(all,aes(x=bfree_w0,fill=type))+
    stat_density(position=position_dodge(0),alpha=0.8,color="black")+
    #scale_y_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #geom_errorbar()+
    #scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    theme_void()+ 
    guides(fill=F)+
    #labs(x='p2',y='bfree_w0')+
    #facet_grid(.~type)+
    coord_flip(xlim=c(0,max(all$bfree_w0)))+
    scale_y_reverse()+
    theme(strip.text.y = element_text(size=14, angle=0),
          plot.margin = margin(15, 2, 30, 10))

  w2=ggplot(all,aes(x=bsa_w2,fill=type))+
    stat_density(position=position_dodge(0),alpha=0.8,color="black")+
    #scale_y_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #geom_errorbar()+
    scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    theme_void() + 
    guides(fill=F)+
    #labs(x='p2',y='bfree_w0')+
    #facet_grid(.~type)+
    coord_flip(xlim=c(1,1000))+
    scale_y_reverse()+
    theme(strip.text.y = element_text(size=14, angle=0),
          plot.margin = margin(5, 2, 30, 10)) 

  bsa_p2=ggplot(all,aes(x=bsa_p2,fill=type))+
    stat_density(position=position_dodge(0),alpha=0.8,color="black")+
    #scale_y_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #geom_errorbar()+
    #scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    theme_classic() + 
    guides(fill=F,x=F,y=F)+
    #labs(x='p2',y='bfree_w0')+
    #facet_grid(.~type)+
    scale_y_continuous(expand = c(0, 0))+
    coord_cartesian(xlim=c(0,0.4))+
    theme(strip.text.y = element_text(size=14, angle=0),
          plot.margin = margin(10, 5, 0, 45),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.line = element_blank(),
          axis.ticks.y=element_blank()) 


  bfree_w1=ggplot(all,aes(x=bfree_w1,fill=type))+
    geom_density(position=position_dodge(0),
                 alpha=0.8,color="black")+
    #(position=position_dodge(0),alpha=0.8,color="black")+
    #scale_y_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    #geom_errorbar()+
    scale_x_log10(minor_breaks=mbreaks,breaks=c(0.001,0.01,0.1,1,10,100,1000,10000,100000,1000000))+
    scale_fill_manual(values = rev(RColorBrewer::brewer.pal(3,"PRGn")))+
    theme_classic() + 
    guides(fill=F,x=F,y=F)+
    #labs(x='p2',y='bfree_w0')+
    #facet_grid(.~type)+
    scale_y_continuous(expand = c(0, 0))+
    coord_cartesian(xlim=c(0.01,10))+
    theme(strip.text.y = element_text(size=14, angle=0),
          plot.margin = margin(10,5, 0, 45),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.line = element_blank(),
          axis.ticks.y=element_blank()) 

  lay <- rbind(c(10,9,9,9,9,9,9,8,8,8,8,8,8),
               c(2,4,4,4,4,4,4,6,6,6,6,6,6),
               c(2,4,4,4,4,4,4,6,6,6,6,6,6),
               c(2,4,4,4,4,4,4,6,6,6,6,6,6),
               c(2,4,4,4,4,4,4,6,6,6,6,6,6),
               c(1,3,3,3,3,3,3,11,7,7,7,7,7),
               c(1,3,3,3,3,3,3,5,5,5,5,5,5),
               c(1,3,3,3,3,3,3,5,5,5,5,5,5),
               c(1,3,3,3,3,3,3,5,5,5,5,5,5)
  )
  venn=
    draw.pairwise.venn(A,
                       B,
                       AB, 
                       lwd=0.1,
                       category = c("Branch-site", "Branch"),
                       fill = c("#af8dc3","#7fbf7b"),
                       col=rep("black",2),
                       alpha = rep(0.8, 2),
                       cat.pos = c(180,-180),
                       cat.dist = rep(0.2, 2))

  g=arrangeGrob(w0,w2,w0_w1,w2_w1,textGrob(""),w2w1_p2,grobTree(venn),bsa_p2,bfree_w1,textGrob(""),textGrob(""),layout_matrix = lay)
  ggsave(filename = paste("../reports/summary_",chosen_stat,".pdf"),g,width=7, height=5)
}