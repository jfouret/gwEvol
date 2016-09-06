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
parser.add_argument('-general_biblio', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git/positiveSelection.bibtex', required=False, help="path of general bibliography file")
parser.add_argument('-gene_rmd', metavar='/path',default='/export/work/batnipah/phylogeny/selection/git/summary.rmd', required=False, help="path of gene rmd file")
parser.add_argument('-adj', metavar='/path', default='yes',required=False, help="yes if adjusment, anything else for no adjustment")
parser.add_argument('-genePythia', metavar='/path', default='/export/scripts/biblio/',required=False, help="script for bibliography")
args=parser.parse_args()
### import
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import os
import sys
import subprocess
import copy
import re
sys.path.append('/export/home/jfouret/lib/')
from myfunctions import *
#variable in args
rootedDir=loadRoot(args.outDir)
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
ete3=Command('ete3 view','ete3 version')
ete3.log()
def viewAlgn(alnFileName,prefix,gene,tree=rootedDir.logs.read('scripts')['positiveSelection.pytree'].value,type='compactseq'):
	global rootedDir
	global ete3
	workdir=rootedDir.reports+'/'+gene
	options={
		'-t':tree,
		'--alg':alnFileName,
		'--alg_format':'fasta',
		'--alg_type':type,
		'-i':prefix+'.pdf'
	}
	cmd=ete3.create(workdir,options)
	submitOneShell(cmd)
def Report(rmdFile,gene,biblio='None',jobID=None):
	global rootedDir
	workdir=rootedDir.reports+'/'+gene
	cmdList=list()
	cmdList.append('cp '+rmdFile+' .')
	if biblio!='None':
		cmdList.append('cp '+biblio+' .')
	cmdList.append('rmd=$(ls *rmd);Rscript -e "rmarkdown::render(\'${rmd}\')"')
	cmdList.append('rm *.rmd')
	submitQsubWithPBS(createPBS(rootedDir.pbs+'/'+gene+'_rep.pbs',cmdList,gene+'_rep',rootedDir.pbsLogs+'/'+gene+'_rep.e',rootedDir.pbsLogs+'/'+gene+'_rep.o',workdir=workdir,waitList=jobID))
def readBEB(geneName):
	global rootedDir
	child=subprocess.Popen('ls '+rootedDir.paml+'/'+geneName+'/*/bsA.*/out',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
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
genepythiaGit=Git(args.genePythia)
genepythia=gitCommand(genepythiaGit,'genePythia.py')
genepythia.log()
reGeneName=re.compile('^(dup[0-9]+_)*([kg_]*[rs_]*[sp_]*)(.*)$')
def byGeneStudy(geneName):
	global args
	global rootedDir
	global genepythia
	global reGeneName
	alnRepo=rootedDir.logs.read('scripts')['positiveSelection.pyalnRepo'].value
	os.chdir(rootedDir.reports)
	workdir=rootedDir.reports+'/'+geneName
	mkdirp(workdir)
	cmdList=['cd '+workdir]
	cmdList.append('cp '+alnRepo+'/'+geneName+'-uc*uniprot.bed uniprot.bed')
	cmdList.append('cp '+alnRepo+'/'+geneName+'-uc*name.bed name.bed')
	cmdList.append('cp '+alnRepo+'/'+geneName+'-uc*fa algn.fa')
	cmdList.append('cp '+alnRepo+'/'+geneName+'-uc*kegg.bed kegg.bed')
        cmdList.append('cp '+alnRepo+'/'+geneName+'-uc*go.bed go.bed')
        cmdList.append('cp '+alnRepo+'/'+geneName+'-uc*alias.bed alias.bed')
	prefixTerm={
		'nipah':'\'nipah,paramyxovirus,paramyxoviridae,hendra,hennipavirus\'',
		'immunoviro':'\'innate immunity,innate immune,virus,viral,immune response\'',
		'metabo':'\'ROS,reactive oxygen species,oxydative stress,free radicals,reactive nitrogen species\'',
		'neuro':'\'brain damage,encephalitis,neuropathology,neurological disorders,neurological deseases\''
	}
	cmdListPBS=['cd '+workdir]
	for prefix,term in prefixTerm.iteritems():
		options={
			'-gene':reGeneName.match(geneName).group(3),
			'-alias':'alias.bed',
			'-prefix':prefix,
			'-term':term
		}
		cmdListPBS.append(genepythia.create(workdir,options))
	submitShell(cmdList)
	jobID=submitQsubWithPBS(createPBS(rootedDir.pbs+'/'+geneName+'_rep.biblio',cmdListPBS,geneName+'_biblio',rootedDir.pbsLogs+'/'+geneName+'_biblio.e',rootedDir.pbsLogs+'/'+geneName+'_biblio.o',workdir=workdir))
	viewAlgn('algn.fa','algn',geneName)

	if geneDict[geneName]=='bs':
		mkdirp(workdir+'/bs')
		os.chdir(workdir+'/bs')
		readBEB(geneName)
		bebfile=open('BEB.tab','r')
		algn=dict()
		fasta=open('../algn.fa','r')
		for line2 in fasta.readlines():
			line2=line2.rstrip()
			if line2[:1]=='>':
				key=line2.lstrip('>')
				algn[key]=''
			else:
				algn[key]+=line2
		fasta.close()
		for line in bebfile.readlines():
			pos=line.rstrip().split("\t")[0]
			print(geneName+':'+pos)
			if pos=='':
				break
			posFile=open(pos+'.uniprot.bed','w')
			with open('../uniprot.bed','r') as uniprot:
				head=None
				for feature in uniprot.readlines():
					featureList=feature.rstrip().split("\t")
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
			alnFile=open(pos+'.aln.tab','w')
			alnFile.write("species\tcodons\tamino acids\n")
			for key in algn:
				codon=algn[key][eval(pos)*3-3:eval(pos)*3]
				try: aa=str(Seq(codon,generic_dna).translate())
				except: aa='N'
				alnFile.write(key+"\t"+codon+"\t"+aa+"\n")
			alnFile.close()
		bebfile.close()
		os.chdir('..')
	Report(args.gene_rmd,geneName,jobID=[jobID])
	os.chdir('..')
cmdGeneral=['cd '+rootedDir.reports]
cmdGeneral.append('cp '+rootedDir.results+'/parameters.tab .')
if args.adj=='yes':
	cmdGeneral.append('cp '+rootedDir.results+'/testsResults_padj.tab geneMethods.tab')
else:
	cmdGeneral.append('cp '+rootedDir.results+'/testsResults_pval.tab geneMethods.tab')
cmdGeneral.append('cp '+rootedDir.results+'/testsResults_raw.tab tests.tab')

if args.adj=='yes':
	cmdGeneral.append('cp '+rootedDir.results+'/venn_padj.pdf venn.pdf')
else:
	cmdGeneral.append('cp '+rootedDir.results+'/venn_pval.pdf venn.pdf')


submitShell(cmdGeneral)

geneDict=dict() # pass, bs, no
head=None
with open(rootedDir.reports+'/geneMethods.tab','r') as statFile:
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
				print(lineList[head['gene_name']])
			elif lineList[head['Branch']]=='TRUE':
				geneDict[lineList[head['gene_name']]]='pass'
			else:
				geneDict[lineList[head['gene_name']]]='no'
Report(args.general_rmd,'',args.general_biblio)
for key in geneDict:
	if geneDict[key]!='no':
		byGeneStudy(key)

