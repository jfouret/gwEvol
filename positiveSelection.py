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
parser.add_argument('-pbsServer', metavar='domain',default='illumina-00-priv', required=False, help="domain name or IP of the pbs server")
parser.add_argument('-pipelineRoot', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git', required=False, help="Root folder of the pipeline for positive selection analysis")
parser.add_argument('-target', metavar='target',default='automatic', required=False, help="Names of target group")
parser.add_argument('-batch', metavar='N',default='20', required=False, help="number of genes per batch")
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
batch=int(args.batch)
treeFile=os.path.abspath(args.tree)
treeMeta=list()
treeMeta.append('No metadatas')
if os.path.exists(args.tree+'.metadata'):
	with open (args.tree+'.metadata', "r") as metadataFile:
		treeMeta=metadataFile.readlines()
alnRepo=os.path.abspath(args.alnRepo)
alnRepoMeta=list()
alnRepoMeta.append('No metadatas')
if os.path.exists(args.alnRepo+'.metadata'):
	with open (args.alnRepo+'.metadata', "r") as metadataFile:
		alnRepoMeta=metadataFile.readlines()

pbsDir=outDir+'/pbs'
logDir=outDir+'/logs'
pbsLogDir=outDir+'/pbsLogs'
pamlDir=outDir+'/paml'

mkdirp(outDir)
mkdirp(pamlDir)
mkdirp(pbsDir)
mkdirp(logDir)
mkdirp(pbsLogDir)

errorFile=open(logDir+'/error','w')
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

#subset
#def runAnalysis(keys,method,byModel=True):
#	alnFile=alnFileDict[keys]
#	if byModel:
#		toSubmit='qsub '+args.pipelineRoot+'/methods/model.pbs -v prot="'+alnFile+'",mark="'+mark+'",tree="'+treeFile+'",prefix="'+prefix+'",model="'+method+'" -d "$PWD" -q '+args.queue+' -N '+keys+'_'+method
#	else:
#		toSubmit='qsub '+args.pipelineRoot+'/methods/'+method+'.pbs -v prot="'+alnFile+'",mark="'+mark+'",tree="'+treeFile+'",prefix="'+prefix+'" -d $PWD -q '+args.queue+' -N '+keys+'_'+method
#
#	child=subprocess.Popen(toSubmit,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#	output,error=child.communicate()
#	if error!='':
#		errorFile.write('WARNING: when submitting '+alnFile+ ' with '+method+' method:'+"\n"+error)
#		return(error)
#	else:
#		return(output)

def submitAnalysis(batchDict,batchName):
	pbsFilePath=pbsDir+'/'+batchName+'.pbs'
	pbsErr=args.pbsServer+':'+pbsLogDir+'/'+batchName+'.pbsErr'
	pbsOut=args.pbsServer+':'+pbsLogDir+'/'+batchName+'.pbsOut'
	qsub='qsub '+pbsFilePath+' -d "$PWD" -q '+args.queue+' -N '+batchName+' -l nodes=1:ppn=1 -o '+pbsOut+' -e '+pbsErr
	start='echo "Job started on $HOSTNAME at time : $(date)"'
	end='echo "Job finished on $HOSTNAME at time : $(date)"'
	cmd=list()
	for keys in batchDict:
		mkdirp(keys)
                cmd.append('cd '+keys)
		for model in modelList:
			cmd.append('ete3 evol --noimg -t '+treeFile+' --alg '+batchDict[keys]+' -o '+prefix+' --mark '+mark+' --models '+model+' --cpu 1 1> '+prefix+'.out 2> '+prefix+'.err'+"\n\n")
		cmd.append('cd ..')
	with open (pbsFilePath, "w") as pbsFile:
		pbsFile.write(start+"\n"+"\n".join(cmd)+"\n"+end+"\n")
	child=subprocess.Popen(qsub,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        output,error=child.communicate()
	if error!='':
                errorFile.write('WARNING: when submitting '+batchName+':'+"\n"+error)
                return(error)
        else:
                return(output)

#se mettre dans le dossier output
os.chdir(pamlDir)

if args.only_branch:
	modelList=list(['b_free','bsA1','bsA'])
else:
	modelList=list(['M0','b_free','M1','bsA1','bsA'])

#Launch jobs
#bsJobs=dict()
#bJobs=dict()
#modelJobs=dict()
#for model in modelList:
	#for keys in alnFileDict:
		#mkdirp(keys)
		#os.chdir(keys)
		#newKeys=keys+"\t"+model
		#bsJobs[keys]=runAnalysis(alnFileDict[keys],'branch')
		#bJobs[keys]=runAnalysis(alnFileDict[keys],'branch_site')
		#modelJobs[newKeys]=runAnalysis(keys,model)
		#os.chdir('..')

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

#print report & metadatas
os.chdir(logDir)
with open ('tree.metadata', "w") as logFile:
	logFile.write("path of the tree: "+treeFile+"\n")
	logFile.write("".join(treeMeta))
shutil.copy(treeFile,'tree')
with open ('alignments.metadata', "w") as logFile:
	logFile.write("path of the repository for alignments: "+alnRepo+"\n")
	logFile.write("".join(alnRepoMeta))
#with open ('branchJobs.txt', "w") as logFile:
#	for keys in bJobs:
#		logFile.write(keys+": "+bJobs[keys])
#with open ('branchSiteJobs.txt', "w") as logFile:
#        for keys in bsJobs:
#                logFile.write(keys+": "+bsJobs[keys])
#with open ('modelJobs.txt','w') as logFile:
#	logFile.write("#GeneName\tModel\tJobId")
#	for keys in modelJobs:
#		logFile.write(keys+"\t"+modelJobs[keys])

with open ('batchJobs.txt','w') as logFile:
	for keys in batchJobs:
		logFile.write(keys+"\t"+batchJobs[keys])
with open ('target.txt','w') as logFile:
	logFile.write("Target species: "+args.mark)

