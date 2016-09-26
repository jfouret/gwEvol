#!/usr/bin/python
import argparse
version=1.0
year=2016
author='Julien Fouret'
contact='julien@fouret.me'
##parse argument
parser = argparse.ArgumentParser(description='put all results in allInOne tableparser',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="path of the output directory from positive selection analysis")
args=parser.parse_args()

import sys
import os
import re
sys.path.append('/export/home/jfouret/lib/')
from myfunctions import *


rootedDir=loadRoot(args.outDir)

os.chdir(rootedDir.reports)

reGene=re.compile('^(kg_|sp_|rs_)*(.*)$')
geneFileName='geneMethods.tab'
goSuffix='/go.bed'
keggSuffix='/kegg.bed'
alnSuffix='/algn_quality.tab'
nameSuffix='/name.bed'
outFileName='allInOne.tab'

outFile=open(outFileName,'w')
headList=['kgID','geneName','prefix','kegg','go','branch','branchsite','aln']
outFile.write("\t".join(headList)+"\n")
geneFile=open(geneFileName,'r')
skipHead=True
for line in geneFile.readlines():
	if skipHead:
		skipHead=False	
	else:
		line=line.rstrip()
		lineList=line.split("\t")
		mGene=reGene.match(lineList[0])
		geneName=mGene.group(2)
		if mGene.group(1)==None:
			prefix=''
		else:
			prefix=mGene.group(1)
		branch=lineList[2]
		branchSite=lineList[3]
		if branch=='FALSE' and  branchSite=='FALSE':
			continue
		with open(lineList[0]+alnSuffix,'r') as alnFile:
			alnFile.readline()
			alnLine=alnFile.readline()
			alnLine=alnLine.rstrip()
			alnLineList=alnLine.split("\t")
			aln=alnLineList[9]
		with open(lineList[0]+nameSuffix,'r') as nameFile:
			nameFile.readline()
			nameLine=nameFile.readline()
			nameLine=nameLine.rstrip()
			nameLineList=nameLine.split("\t")
			kgID=nameLineList[1]
		with open(lineList[0]+keggSuffix,'r') as keggFile:
			keggSkipHead=True
			keggIDList=list()
			for keggLine in keggFile.readlines():
				if keggSkipHead:
					keggSkipHead=False
				else:
					keggLine=keggLine.rstrip()
					keggLineList=keggLine.split("\t")
					keggIDList.append(keggLineList[1])
		keggIDs=';'.join(keggIDList)
		with open(lineList[0]+goSuffix,'r') as goFile:
			goSkipHead=True
			goIDList=list()
			for goLine in goFile.readlines():
				if goSkipHead:
					goSkipHead=False
				else:
					goLine=goLine.rstrip()
					goLineList=goLine.split("\t")
					goIDList.append(goLineList[1])
		goIDs=';'.join(goIDList)
		writeList=[kgID,geneName,prefix,keggIDs,goIDs,branch,branchSite,aln]
		outFile.write("\t".join(writeList)+"\n")
outFile.close()
geneFile.close()


