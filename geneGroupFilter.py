#!/usr/bin/python
import argparse
version=1.0
year=2016
author='Julien Fouret'
contact='julien@fouret.me'

#parse argument
parser = argparse.ArgumentParser(description='filter kegg.tab and go.tab from allInOne table to keep only genes with a good alignment',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="path of the output directory from positive selection analysis")
args=parser.parse_args()

#Import libraries
import sys
import os
import re
sys.path.append('/export/home/jfouret/lib/')
from myfunctions import *

#Load rootDir
rootedDir=loadRoot(args.outDir)
os.chdir(rootedDir.reports)

#Variable definition
keggIn='kegg.tab'
keggOut='kegg_filter.tab'
goIn='go.tab'
goOut='go_filter.tab'
allInOneFileName='allInOne.tab'

## step 1 - get the list of genes

geneList=list()
with open(allInOneFileName,'r') as allInOneFile:
	allInOneFile.readline()
	for line in allInOneFile.readlines():
		line=line.rstrip()
		lineList=line.split("\t")
		geneName=lineList[1]
		if lineList[3]=='no':
			dup=''
		else:
			dup='dup'+lineList[3]+'_'
		prefix=lineList[2]
		gene=dup+prefix+geneName
		if lineList[8]=='YES':
			geneList.append(gene)

# step 2.0 - Function to filter a commat separated list
def filterGenes(genesString,refList=geneList):
	oldList=genesString.split(',')
	newList=list()
	for gene in oldList:
		if gene in refList:
			newList.append(gene)
	if len(newList)==0:
		return '__NO__'
	else:
		return ','.join(newList)

# step 2.1 - filter kegg
keggOutFile=open(keggOut,'w')
with open(keggIn,'r') as keggInFile:
	header=keggInFile.readline()
	keggOutFile.write(header)
	for line in keggInFile.readlines():
		line=line.rstrip()
		lineList=line.split("\t")
		newGenes=filterGenes(lineList[2])
		if newGenes!='__NO__':
			keggOutFile.write(lineList[0]+"\t"+lineList[1]+"\t"+newGenes+"\n")
keggOutFile.close()

# step 2.2 - filter go
goOutFile=open(goOut,'w')
with open(goIn,'r') as goInFile:
        header=goInFile.readline()
        goOutFile.write(header)
        for line in goInFile.readlines():
                line=line.rstrip()
                lineList=line.split("\t")
                newGenes=filterGenes(lineList[3])
                if newGenes!='__NO__':
                        goOutFile.write(lineList[0]+"\t"+lineList[1]+"\t"+lineList[2]+"\t"+newGenes+"\n")
goOutFile.close()

# save and quit
saveRoot(rootedDir)
sys.exit()