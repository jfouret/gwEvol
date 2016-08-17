#!/usr/bin/python

import argparse

version=1.0
year=2016
author='Julien Fouret'
contact='julien@fouret.me'

##parse argument
parser = argparse.ArgumentParser(description='Make reports for phylo',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="path of the output directory from positive selection analysis")
parser.add_argument('-general_rmd', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git/report.rmd', required=False, help="path of general rmd file")
parser.add_argument('-gene_rmd', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git/summary.rmd', required=False, help="path of gene rmd file")
parser.add_argument('-repDir', metavar='name', default='report',required=False, help="report dir")
parser.add_argument('-cpu', metavar='N', default='16',required=False, help="number of threads to use")
parser.add_argument('-adj', metavar='/path', default='yes',required=False, help="yes if adjusment, anything else for no adjustment")
parser.add_argument('-genePythia', metavar='/path', default='/export/scripts/biblio/genePythia.py',required=False, help="script for bibliography")
args=parser.parse_args()


### import
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import os
import sys
import subprocess
import multiprocessing as mp
import copy
import errno

#variable in args
mainDir=os.path.abspath(args.outDir)

###Put in a separate function
def getOpt(fileName,module,option):
	head = None
	with open(fileName,'r') as file:
		for line in file.readlines():
			if head is None:
				head=dict()
				index=0
				for var in line.rstrip().split(';'):
					head[var]=index
					index+=1
				continue
			lineList=line.rstrip().split(';')
			if lineList[head['module']]==module and lineList[head['option']]==option:
				return(lineList[head['value']])
				break
###

def viewAlgn(alnFileName,prefix,tree=mainDir+'/logs/tree',type='compactseq'):
	global mainDir
	#'+alnDir+'/'+geneName+'*.fa'
	cmd='ete3 view -t '+tree+' --alg '+alnFileName+' --alg_format fasta --alg_type '+type+' -i '+prefix+'.pdf'
	try:
		child=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		output,error=child.communicate()
        #print(error)
	except:
		cmd="echo 'error drawing the alignments' > t0."+prefix+".pdf"
		child=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
		output2,error2=child.communicate()
		

def Report(rmdFile):
	cmdList=list()
	cmdList.append('cp '+rmdFile+' .')
	cmdList.append('rmd=$(ls *rmd);Rscript -e "rmarkdown::render(\'${rmd}\')"')
	cmdList.append('rm *.rmd')
	child=subprocess.Popen(';'.join(cmdList),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	output,error=child.communicate()
	#print(error)

def mkdirp(path):
	try:
		os.makedirs(path)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise

def readBEB(geneName):
	global mainDir
	child=subprocess.Popen('ls '+mainDir+'/paml/'+geneName+'/*/bsA.*/out',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	output,error=child.communicate()
	fileName=output.rstrip()
	file=open(fileName,'r')
	BEB=""
	level=0
	for line in file.readlines():
		line=line.rstrip()
		if level==3:
			break
		if level==2:
			if line=='':
				level=3
			BEB=BEB+"\t".join(line.lstrip().split(' '))+"\n"
		if level==1:
			level=2
		if line=='Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)':
			level=1
	file.close()
	with open('BEB.tab','w') as file:
		file.write(BEB)

def byGeneStudy(geneName):
	global args
	global mainDir
	os.chdir(mainDir+'/'+args.repDir)
	mkdirp(geneName)
	#print(geneName)
	os.chdir(geneName)
	cmdList=list()
	cmdList.append('cp '+getOpt(mainDir+'/logs/metadata','positiveSelection.py','alnRepo')+'/'+geneName+'-*uc*uniprot.bed uniprot.bed')
	cmdList.append('cp '+getOpt(mainDir+'/logs/metadata','positiveSelection.py','alnRepo')+'/'+geneName+'-*uc*name.bed name.bed')
	cmdList.append('cp '+getOpt(mainDir+'/logs/metadata','positiveSelection.py','alnRepo')+'/'+geneName+'-*uc*fa algn.fa')
	cmdList.append('cp '+getOpt(mainDir+'/logs/metadata','positiveSelection.py','alnRepo')+'/'+geneName+'-*uc*kegg.bed kegg.bed')
        cmdList.append('cp '+getOpt(mainDir+'/logs/metadata','positiveSelection.py','alnRepo')+'/'+geneName+'-*uc*go.bed go.bed')
        cmdList.append('cp '+getOpt(mainDir+'/logs/metadata','positiveSelection.py','alnRepo')+'/'+geneName+'-*uc*alias.bed alias.bed')
	cmdList.append(args.genePythia+' -gene '+geneName+' -alias alias.bed -prefix immunity -term immune,immunity,interferon,immuno,immunlogy'
	cmdList.append(args.genePythia+' -gene '+geneName+' -alias alias.bed -prefix virology -term nipah,virus,paramyxovirus,viral,virology'
	cmdList.append(args.genePythia+' -gene '+geneName+' -alias alias.bed -prefix metabolism -term ROS,metabolism'
	child=subprocess.Popen(';'.join(cmdList),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
	outpur,error=child.communicate()
	viewAlgn('algn.fa','algn')
	if geneDict[geneName]=='bs':
		mkdirp('bs')
		os.chdir('bs')
		readBEB(geneName)
		file=open('BEB.tab','r')
		algn=dict()
		fasta=open('../algn.fa','r')
		for line2 in fasta.readlines():
			line2=line2.rstrip()
			if line2[:1]=='>':
				key=line2.lstrip('>')
				algn[key]=''
			else:
				algn[key]+=line2
		for line in file.readlines():
			pos=line.rstrip().split("\t")[0]
			print(geneName+':'+pos)
			if pos=='':
				break
			posFile=open(pos+'.uniprot.bed','w')
			with open('../uniprot.bed','r') as uniprot:
				head=None
				for feature in uniprot.readlines():
					featureList=feature.rstrip().split(';')
					if head==None:
						posFile.write(feature)
						head=dict()
						index=0
						for var in featureList:
							head[var]=index
							index+=1
					elif eval(pos)>eval(featureList[head['start']]) and eval(pos)<=eval(featureList[head['end']]):
						posFile.write(feature)
			posFile.close()
			nuclFile=open(pos+'.nucl.fa','w')
			aaFile=open(pos+'.aa.fa','w')
			for key in algn:
				codon=algn[key][eval(pos)*3-3:eval(pos)*3]
				try: aa=str(Seq(codon,generic_dna).translate())
				except: aa='N'
				nuclFile.write('>'+key+"\n"+codon+"\n")
				aaFile.write('>'+key+"\n"+aa+"\n")
			nuclFile.close()
			fasta.close()
			aaFile.close()
			viewAlgn(pos+'.nucl.fa',pos+'.nucl',type='fullseq')
			viewAlgn(pos+'.aa.fa',pos+'.aa',type='fullseq')
		file.close()
		os.chdir('..')
	Report(args.gene_rmd)	
	os.chdir('..')

os.chdir(mainDir)
mkdirp(args.repDir)
os.chdir(args.repDir)
cmdGeneral=list()
cmdGeneral.append('cp '+mainDir+'/results/parameters.tab .')
if args.adj=='yes':
	cmdGeneral.append('cp '+mainDir+'/stats/testsResults_padj.tab geneMethods.tab')
else:
	cmdGeneral.append('cp '+mainDir+'/stats/testsResults_pval.tab geneMethods.tab')
cmdGeneral.append('cp '+mainDir+'/stats/testsResults_raw.tab tests.tab')

if args.adj=='yes':
	cmdGeneral.append('cp '+mainDir+'/stats/venn_padj.pdf venn.pdf')
else:
	cmdGeneral.append('cp '+mainDir+'/stats/venn_pval.pdf venn.pdf')

child=subprocess.Popen(';'.join(cmdGeneral),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
output,error=child.communicate()

geneDict=dict() # pass, bs, no
head=None
with open('geneMethods.tab','r') as statFile:
	for line in statFile.readlines():
		lineList=line.rstrip().split("\t")
		if head==None:
			head=dict()
			index=0
			for var in lineList:
				head[var]=index
				index+=1
		else:
			if lineList[head['BranchSites']]=='TRUE':
				geneDict[lineList[head['gene_name']]]='bs'
			elif lineList[head['Branch']]=='TRUE':
				geneDict[lineList[head['gene_name']]]='pass'
			else:
				geneDict[lineList[head['gene_name']]]='no'

Report(args.general_rmd)

geneDictCopy=copy.deepcopy(geneDict)

for key in geneDictCopy:
	if geneDict[key]=='no':
		del geneDict[key]



pool=mp.Pool(processes=eval(args.cpu))
res=pool.map_async(byGeneStudy,geneDict.keys())
res.wait()
pool.close()
pool.join()




