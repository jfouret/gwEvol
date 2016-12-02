#!/usr/bin/python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien@fouret.me'
##parse argument
parser = argparse.ArgumentParser(description='put all results in allInOne tableparser',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="path of the output directory from positive selection analysis")
parser.add_argument('-aln_dir', metavar='/path', required=True, help="Aln repository in with annotation")
args=parser.parse_args()

import sys
import os
import re
from jupype import *

rootedDir=loadRoot(args.outDir)
os.chdir(rootedDir.reports)
reGene=re.compile('^(dup)*([0-9]*)(_)*(kg_|sp_|rs_)*(.*)$')
geneFileName='geneMethods.tab'

goPath='go.tab'
keggPath='kegg.tab'
namePath='name.tab'

outFileName='allInOne.tab'
outFile=open(outFileName,'w')
headList=['kgID','geneName','prefix','duplicate','kegg','go','branch','branchsite','aln']
outFile.write("\t".join(headList)+"\n")

geneFile=open(geneFileName,'r')
skipHead=True

#TODO DEPOSITORY
# + ln kegg et go tab
alnRepo=os.path.abspath(args.aln_dir)
submitOneShell('ln -s '+alnRepo+'/../reports/kegg.tab .')
submitOneShell('ln -s '+alnRepo+'/../reports/go.tab .')

for line in geneFile.readlines():
	if skipHead:
		skipHead=False	
	else:
		line=line.rstrip()
		lineList=line.split("\t")
		mGene=reGene.match(lineList[0])
		geneName=mGene.group(5)
		if mGene.group(4)==None:
			prefix=''
		else:
			prefix=mGene.group(4)
		if mGene.group(1)==None:
			duplicate='no'
		else:
			duplicate=mGene.group(2)
		branch=lineList[2]
		branchSite=lineList[3]
		print(lineList[0])
		if branch=='FALSE' and  branchSite=='FALSE':
			sourceAlnNameFile=submitOneShell('ls '+alnRepo+'/'+lineList[0]+'-uc*.fa')['out'].rstrip()
			alnFileName=alnSuffix
			cmdList=['cd '+rootedDir.reports]
			cmdList.append('ln -sf '+sourceAlnNameFile+' ./algn.fa')
			checkOpt={
				'-repDir':rootedDir.reports,
				'-fore':'25',
				'-back':'50',
				'-targets':rootedDir.logs.read('scripts')['positiveSelection.pymark'].value,
			}
			cmdList.append(checkAlign.create(rootedDir.reports,checkOpt))
			submitShell(cmdList)
			nameFileName=submitOneShell('ls '+alnRepo+'/'+lineList[0]+'-uc*'+nameSuffix)['out'].rstrip()
			keggFileName=submitOneShell('ls '+alnRepo+'/'+lineList[0]+'-uc*'+keggSuffix)['out'].rstrip()
			goFileName=submitOneShell('ls '+alnRepo+'/'+lineList[0]+'-uc*'+goSuffix)['out'].rstrip()
		else:
			alnFileName=lineList[0]+'/'+alnSuffix
			nameFileName=lineList[0]+'/'+nameSuffix
			keggFileName=lineList[0]+'/'+keggSuffix
			goFileName=lineList[0]+'/'+goSuffix

		with open(alnFileName,'r') as alnFile:
			alnFile.readline()
			alnLine=alnFile.readline()
			alnLine=alnLine.rstrip()
			alnLineList=alnLine.split("\t")
			aln=alnLineList[9]
		with open(nameFileName,'r') as nameFile:
			nameFile.readline()
			nameLine=nameFile.readline()
			nameLine=nameLine.rstrip()
			nameLineList=nameLine.split("\t")
			kgID=nameLineList[1]
		with open(keggFileName,'r') as keggFile:
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
		with open(goFileName,'r') as goFile:
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
		writeList=[kgID,geneName,prefix,duplicate,keggIDs,goIDs,branch,branchSite,aln]
		outFile.write("\t".join(writeList)+"\n")
outFile.close()
geneFile.close()

# exit
saveRoot(rootedDir)
sys.exit()
