#!/usr/bin/Rscript
author='Julien FOURET'
contact='julien.fouret@fouret.me'
year='2016'
version='SEDMATCHGITVERSION'

startMessage=paste('Perform Log-likelihood ratio tests for positive selection analysis '," ; Version : ",version," ; ",year," ; Author : ",author," for more informations or enquiries please contact ",contact," (WARNING: pval and FC column names must begin respectively with 'p.val' and 'FC')",sep='')
suppressMessages(library("argparse"))
parser <- ArgumentParser(description=startMessage)
parser$add_argument("-outDir",
                    type="character",
                    metavar="/path",
                    help="output directory for the whole analysis")
parser$add_argument("-evolPack",
                    type="character",
		    default='/export/scripts/phylogeny/evolPack.R',
                    metavar="/path",
                    help="Path of the evolPack package")
parser$add_argument("-vennPack",
                    type="character",
                    default='/export/scripts/utils/venn.R',
  		    metavar="/path",
                    help="Path of the vennPack package")
parser$add_argument("-pval",
                    type="character",
                    default='0.05',
                    help="P-value cut-off for positively selected genes")
parser$add_argument("-graphSign",
                    action="store_true",
                    default='FALSE',
                    help="graph parameters distribution for all significant results")
args <- parser$parse_args()

# LOADING LIBRARIES
suppressMessages(library(ggthemes))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(gplots))
suppressMessages(library(tidyr))
suppressMessages(source(args$evolPack))
suppressMessages(source(args$vennPack))
suppressMessages(library(R.utils))

#DEFINE PARAMETERS
outDir=getAbsolutePath(args$outDir)
resDir=paste(outDir,'results',sep='/')
statDir=resDir
setwd(statDir)
parametersFileName=paste(resDir,'parameters.tab',sep='/')
testResultsFileName=paste(statDir,'testsResults_',sep='/')

#IMPORT DATAFRAME
data=read.table(file = parametersFileName,header = T,sep = "\t",na.strings = "NA",dec = "." )
data=unique(data)
data$fmodel=factor(data$model,levels = c('b_free','M0','bsA','bsA1','M1'))

#put in parameters adjustments #TODO
listdf=evol_tests(data,alpha=as.numeric(args$pval))
#### listdf 1:tests pval and padj   2: branch and branch-site results for pval     3: branch and branch-site results for padj
i=1
namefile=c('pval','padj')
for (allTEST in listdf[2:3]){
  bANDbs=subset(allTEST,(Branch==T)&(BranchSites==T))
  b=subset(allTEST,(Branch==T))
  bs=subset(allTEST,(BranchSites==T))
  bsONLY=subset(allTEST,(Branch==F)&(BranchSites==T))
  bONLY=subset(allTEST,(Branch==T)&(BranchSites==F))
  venn=venn2(names = c('Branch','Branch-sites'),col=c('blue',"orangered"),weigths = c(dim(bONLY)[1],dim(bANDbs)[1],dim(bsONLY)[1]))
  ggsave(paste('venn_',namefile[i],'.pdf',sep=''),venn,width = 8.27, height = 5.83, units = "in")
  write.table(allTEST,file=paste(testResultsFileName,namefile[i],'.tab',sep=''),quote=F,row.names=F,sep="\t")
  i=i+1
}
write.table(listdf[[1]],file=paste(testResultsFileName,'raw','.tab',sep=''),quote=F,row.names=F,sep="\t")

# graph all parameters if needed
if (args$graphSign==TRUE){
	dir.create('branch_specific',recursive=T)
	dir.create('branch-site_specific',recursive=T)
	dir.create('aspecific',recursive=T)
	gene_target=unique(cbind(bONLY$gene_name,bONLY$target))
	for (i in 1:dim(gene_target)[1]){
    		file=paste('branch_specific/',gene_target[i,1],'_',gene_target[i,2],'.pdf',sep='')
    		graph=makePhyloGraph(data,gene_target[i,1],gene_target[i,2])
    		ggsave(file,graph, width = 11.69, height = 8.27, units = "in")
	}

	gene_target=unique(cbind(bsONLY$gene_name,bsONLY$target))
	for (i in 1:dim(gene_target)[1]){
  		file=paste('branch-site_specific/',gene_target[i,1],'_',gene_target[i,2],'.pdf',sep='')
  		graph=makePhyloGraph(data,gene_target[i,1],gene_target[i,2])
  		ggsave(file,graph, width = 11.69, height = 8.27, units = "in")
	}

	gene_target=unique(cbind(bANDbs$gene_name,bANDbs$target))
	for (i in 1:dim(gene_target)[1]){
  		file=paste('aspecific/',gene_target[i,1],'_',gene_target[i,2],'.pdf',sep='')
  		graph=makePhyloGraph(data,gene_target[i,1],gene_target[i,2])
  		ggsave(file,graph, width = 11.69, height = 8.27, units = "in")
	}
}
