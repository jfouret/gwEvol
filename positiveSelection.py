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
parser.add_argument('-only_branch',default=False, action='store_true', help="perform only branch model (if you have already done others)")
parser.add_argument('-pipelineRoot', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git', required=False, help="Root folder of the pipeline for positive selection analysis")
parser.add_argument('-target', metavar='target',default='automatic', required=False, help="Names of target group")
parser.add_argument('-batch', metavar='N',default='20', required=False, help="number of genes per batch")
parser.add_argument('-queue', metavar='name',default='long7', required=False, help="Queue for qsub")
parser.add_argument('-subset', metavar='/path',default='None', required=False, help="file with a list of gene name to consider (only)")
parser.add_argument('-tree', metavar='tree.nh', required=True, help="Tree file with newick format")

args=parser.parse_args()

import re
import os
import subprocess
import sys
sys.path.append('/export/home/jfouret/lib/')
from myfunctions import *

rootedDir=RootDir(args.outDir,True,{'paml':str})
rootedDir.logs.writeArgs(args)
batch=int(args.batch)
treeFile=os.path.abspath(args.tree)
alnRepo=os.path.abspath(args.alnRepo)
rootedDir.logs.addMeta('Tree',treeFile.rstrip('/')+'.metadata')
rootedDir.logs.addMeta('Alignments',alnRepo.rstrip('/')+'.metadata')
markList=args.mark.split(',')
mark=',,'.join(markList)
if args.target=='automatic':
	prefix='_'.join(markList)
else:
	prefix=args.target
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
os.chdir(alnRepo)
rName='/([a-zA-Z0-9_.@ :\-]*)-(uc[^/]*)\.fa$'
alnFileDict=dict()
for file in os.listdir(alnRepo):
        fileAbs=os.path.abspath(file)
	m=re.search(rName,fileAbs)
	if m:
		key=m.group(1)
		if subset:
			if key in subsetList:
				alnFileDict[key]=fileAbs
		else:
			alnFileDict[key]=fileAbs

ete3=Command('ete3 evol','ete3 version')

def submitAnalysis(batchDict,batchName):
	global rootedDir
	pbsFilePath=rootedDir.pbs+'/'+batchName+'.pbs'
	pbsErr=rootedDir.pbsLogs+'/'+batchName+'.e'
	pbsOut=rootedDir.pbsLogs+'/'+batchName+'.o'
	end='echo "Job finished on $HOSTNAME at time : $(date)"'
	cmd=list()
	for keys in batchDict:
		mkdirp(rootedDir.paml+'/'+keys)
		opt={
			'--noimg':'',
			'-t':treeFile,
			'--alg':batchDict[keys],
			'-o':prefix,
			'--mark':mark,
			'--cpu':'1'
		}
		pos=['1> '+prefix+'.out','2> '+prefix+'.err']
		for model in modelList:
			opt['--models']=model
			cmd.append(ete3.create(rootedDir.paml+'/'+keys,opt,pos).replace('@','\@'))
	#submitQsubWithPBS(createPBS(pbsFilePath,cmd,batchName,pbsErr,pbsOut,queue=args.queue,workdir=rootedDir.paml))
	createPBS(pbsFilePath,cmd,batchName,pbsErr,pbsOut,queue=args.queue,workdir=rootedDir.paml)
os.chdir(rootedDir.paml)

if args.only_branch:
	modelList=list(['b_free','bsA1','bsA'])
else:
	modelList=list(['M0','b_free','M1','bsA1','bsA'])

batchJobs=dict()
count=0
batchNumber=1
batchDict=dict()
lastKey=alnFileDict.keys()[-1]
for keys in alnFileDict:
	count+=1
	batchDict[keys]=alnFileDict[keys]
	if count==batch or keys==lastKey:
		batchName='Batch_'+str(batchNumber)
		count=0
		batchNumber+=1
		batchJobs[batchName]=submitAnalysis(batchDict,batchName)
		batchDict=dict()
saveRoot(rootedDir)

