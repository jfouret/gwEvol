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

#put in parameters adjustments #TODO
data=evol_tests(param)

write.table(data,file=paste(testResultsFileName,'raw','.tab',sep=''),quote=F,row.names=F,sep="\t")



### Plot results with pval
chosen_stat="pval"
sdata=data[,c("tests","gene_name","target",chosen_stat)]
names(sdata)[4]='stat'
sdata=spread(sdata,key=tests,value = stat)
sdata$branch_ratio=mutate(subset(param,model=="b_free"),b_free_ratio=w1/w0)$b_free_ratio
branch=subset(sdata,(`b_free vs M0`<0.05)&(`branch_ratio`>1))$gene_name
branch_site=subset(sdata,(`bsA vs bsA1`<0.05)&(`bsA vs M1`<0.05))$gene_name
vpval=venn.diagram(list("Branch"=branch,"Branch-site"=branch_site),height = 8, width = 8,fill = c("red", "green"),
             alpha = c(0.5, 0.5), cat.pos=c(0,45),cex = 2.5,cat.cex=2,cat.fontface = 4,lty =2, fontfamily =3, 
             filename = "pval_venn.svg",imagetype='svg',main="<0.05 pval",main.cex=3);
write(paste(branch,collapse = "\n"),"branch_pval.txt")
write(paste(branch_site,collapse = "\n"),"branch_site_pval.txt")

### Plot results with padj
chosen_stat="padj"
sdata=data[,c("tests","gene_name","target",chosen_stat)]
names(sdata)[4]='stat'
sdata=spread(sdata,key=tests,value = stat)
sdata$branch_ratio=mutate(subset(param,model=="b_free"),b_free_ratio=w1/w0)$b_free_ratio
branch=subset(sdata,(`b_free vs M0`<0.05)&(`branch_ratio`>1))$gene_name
branch_site=subset(sdata,(`bsA vs bsA1`<0.05)&(`bsA vs M1`<0.05))$gene_name
vpadj=venn.diagram(list("Branch"=branch,"Branch-site"=branch_site), height = 8, width = 8,fill = c("red", "green"),
                   alpha = c(0.5, 0.5), cat.pos=c(0,45),cex = 2.5,cat.cex=2,cat.fontface = 4,lty =2, fontfamily =3, 
                   filename = "padj_venn.svg",imagetype='svg',main="<0.05 padj (Benjamini-Hochberg) ",main.cex=3);
write(paste(branch,collapse = "\n"),"branch_padj.txt")
write(paste(branch_site,collapse = "\n"),"branch_site_padj.txt")
