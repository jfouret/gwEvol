#!/usr/bin/python

import argparse

version=1.0
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'

parser = argparse.ArgumentParser(description='perform positive selection analysis via qsub submitting',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#parser.add_argument('-', metavar='', required=True, help="")
#parser.add_argument('-', metavar='',default='', required=False, help="")

parser.add_argument('-alnRepo', metavar='/path', required=True, help="alignement repositories (a folder with all alignments for each gene name named 'GeneName-ID.fa' or '.fna' or '.fasta')")
parser.add_argument('-mark', metavar='spec1,...,specN', required=True, help="branch(es) to mark as target for positive selection analysis (using names in the tree)")
parser.add_argument('-outDir', metavar='/path', required=True, help="Output directory")
parser.add_argument('-pipelineRoot', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git', required=False, help="Root folder of the pipeline for positive selection analysis")
parser.add_argument('-target', metavar='target',default='automatic', required=False, help="Names of target group")
parser.add_argument('-queue', metavar='name',default='long7', required=False, help="Queue for qsub")
parser.add_argument('-subset', metavar='/path',default='None', required=False, help="file with a list of gene name to consider (only)")
parser.add_argument('-tree', metavar='tree.nh', required=True, help="Tree file with newick format")

args=parser.parse_args()

import shutil
import re
import os
import subprocess

def mkdirp(path):
	if not os.path.exists(path):
		os.makedirs(path)

outDir=os.path.abspath(args.outDir)
treeFile=os.path.abspath(args.tree)
treeMeta=list('No metadatas')
if os.path.exists(args.tree+'.metadata'):
	with open (args.tree+'.metadata', "r") as metadataFile:
		treeMeta=metadataFile.readlines()
alnRepo=os.path.abspath(args.alnRepo)
alnRepoMeta=list('No metadatas')
if os.path.exists(args.alnRepo+'.metadata'):
	with open (args.alnRepo+'.metadata', "r") as metadataFile:
		alnRepoMeta=metadataFile.readlines()
mkdirp(outDir+'/logs')
errorFile=open(outDir+'/logs/error','w')
markList=args.mark.split(',')
mark='bbb'.join(markList)
#bbb used on mathods script for as separator
if args.target=='automatic':
	prefix='_'.join(markList)

#list of gene for subsetting
subset=False
if args.subset!='None':
	subset=True
	subsetList=list()
	subsetFile=open(args.subset,'r')
	for line in subsetFile.readlines():
		geneToSubset=line.rstrip("\n")
		subsetList.append(geneToSubset)
	subsetFile.close()

#Creer un dico pour les alnFile
rName='/([a-zA-Z0-9_.@ :]*)-[^/]*.$'
alnFileDict=dict()
for file in os.listdir(alnRepo):
	if file.endswith(".fa") or file.endswith(".fasta") or file.endswith(".fna"):
        	fileAbs=os.path.abspath(file)
		m=re.search(rName,fileAbs)
		if m:
			key=m.group(1)
			if subset:
				if key in subsetList:
					alnFileDict[key]=fileAbs
			else:
				alnFileDict[key]=fileAbs
		else:
			errorFile.write("WARNING: "+fileAbs+" do not allow retrieval of gene Name '"+rName+"'\n")


#Subset in alnFileDict based on keys
#

# definir les processeurs dans la methodes
def runAnalysis(alnFile,method):
	toSubmit='qsub '+args.pipelineRoot+'/methods/'+method+'.pbs -v prot="'+alnFile+'",mark="'+mark+'",tree="'+treeFile+'",prefix="'+prefix+'" -d $PWD -q '+args.queue+' -N '+alnFile+'_'+method
	child=subprocess.Popen(toSubmit,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	output,error=child.communicate()
	if error!='':
		errorFile.write('WARNING: when submitting '+alnFile+ ' with '+method+' method:'+"\n"+error)
		return(error)
	else:
		return(output)

#Se mettre dans le dossier output
mkdirp(outDir)
os.chdir(outDir)

#Launch jobs
bsJobs=dict()
bJobs=dict()
for keys in alnFileDict:
	mkdirp(keys)
	os.chdir(keys)
	bsJobs[keys]=runAnalysis(alnFileDict[keys],'branch')
	bJobs[keys]=runAnalysis(alnFileDict[keys],'branch_site')
	os.chdir('..')

#print report & metadatas
os.chdir('logs')
with open ('tree.metadata', "w") as logFile:
	logFile.write("\n".join(treeMeta))
shutil.copy(treeFile,'tree')
with open ('alignments.metadata', "w") as logFile:
	logFile.write("\n".join(alnRepoMeta))	
with open ('branchJobs.txt', "w") as logFile:
	for keys in bJobs:
		logFile.write(keys+": "+bJobs[keys]+"\n")
with open ('branchSiteJobs.txt', "w") as logFile:
        for keys in bsJobs:
                logFile.write(keys+": "+bsJobs[keys]+"\n")
with open ('target.txt','w') as logFile:
	logFile.write('Target species: '+args.mark)



