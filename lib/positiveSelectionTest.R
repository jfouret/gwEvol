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
parser$add_argument("-pval",
                    type="character",
                    default='0.05',
                    help="P-value cut-off for positively selected genes")
parser$add_argument("-graphSign",
                    action="store_true",
                    default='FALSE',
                    help="graph parameters distribution for all significant results")
parser$add_argument("-evolPack",
                    type="character",
                    help="evolution package")
parser$add_argument("-venn",
                    type="character",
                    help="venn package")
args <- parser$parse_args()

# LOADING LIBRARIES
suppressMessages(library(ggthemes))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(gplots))
suppressMessages(library(tidyr))
suppressMessages(source(args$evolPack))
suppressMessages(source(args$venn))
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
data_evol=evol_tests(data)

write.table(data_evol,file=paste(testResultsFileName,'raw','.tab',sep=''),quote=F,row.names=F,sep="\t")
