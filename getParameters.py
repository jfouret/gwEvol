#!/usr/bin/python

import argparse

version=1.1
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'

##parse argument
parser = argparse.ArgumentParser(description='allow',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="path of the output directory from positive selection analysis")

args=parser.parse_args()

import re
import sys
import glob
import os

#function
def mkdirp(path):
        if not os.path.exists(path):
                os.makedirs(path)

def read_model(path,model,target,name):
	file=open(path,'r')
	[w0,w1,w2,p0,p1,p2a,p2b,lnL]=['err','err','err','err','err','err','err','err']
	for line in file:

		newline=line.rstrip('\n')
		rlnL = re.compile('^lnL\(.*\):[\s]+([^\s]*)')
		rkappa = re.compile('^kappa.*=[\s]+([^\s]*)')

		if rlnL.match(newline):
			lnL=rlnL.match(newline).group(1)
		
		if rkappa.match(newline):
			kappa=rkappa.match(newline).group(1)
		
		if model == 'b_free':
			rw = re.compile('^w\s*\(dN.*for branches:\s*([^\s]*)\s([^\s]*)')
			if rw.match(newline):
				w0=rw.match(newline).group(1)
				w1=rw.match(newline).group(2)
		
		elif model == 'M1':
			rw = re.compile('^w:\s+([^\s]+)\s+([^\s]+)')
			if rw.match(newline):
				w0=rw.match(newline).group(1)
				w1=rw.match(newline).group(2)

			rp = re.compile('^p:\s+([^\s]+)\s+([^\s]+)')
			if rp.match(newline):
				p0=rp.match(newline).group(1)
				p1=rp.match(newline).group(2)
			
		elif model == 'bsA':
			rw = re.compile('^foreground\s+w\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)')
			if rw.match(newline):
				w0=rw.match(newline).group(1)
				w1=rw.match(newline).group(2)
				w2=rw.match(newline).group(3)
				
			rp = re.compile('^proportion\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)')
			if rp.match(newline):
				p0=rp.match(newline).group(1)
				p1=rp.match(newline).group(2)
				p2a=rp.match(newline).group(3)
				p2b=rp.match(newline).group(4)
			
		elif model == 'bsA1':
			rw = re.compile('^foreground\s+w\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)')
			if rw.match(newline):
				w0=rw.match(newline).group(1)
				w1=rw.match(newline).group(2)
				w2=rw.match(newline).group(3)
				
			rp = re.compile('^proportion\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)')
			if rp.match(newline):
				p0=rp.match(newline).group(1)
				p1=rp.match(newline).group(2)
				p2a=rp.match(newline).group(3)
				p2b=rp.match(newline).group(4)
			
		else:
			romega = re.compile('^omega.*=[\s]*([^\s]*)')
			if romega.match(newline):
				omega=romega.match(newline).group(1)
	file.close
	if model == 'b_free':
		toReturn=[name,model,target,w0,w1,'NA','NA','NA','NA','NA',lnL]
	elif model == 'M1':
		toReturn=[name,model,'NA',w0,w1,'NA',p0,p1,'NA','NA',lnL]
	elif model == 'bsA':
		toReturn=[name,model,target,w0,w1,w2,p0,p1,p2a,p2b,lnL]
	elif model == 'bsA1':	
		toReturn=[name,model,target,w0,w1,w2,p0,p1,p2a,p2b,lnL]
	else:
		toReturn=[name,model,'NA',omega,'NA','NA','NA','NA','NA','NA',lnL]
	if 'err' in toReturn:
		errorFile.write('WARNING parsing results of gene '+name+' in target '+target+' with model '+model+"\n")
	return(toReturn)

#TODO add metadatas

outDir=os.path.abspath(args.outDir)

pamlDir=outDir+'/paml/'
resDir=outDir+'/results/'

mkdirp(resDir)

resultFile=resDir+'parameters.tab'
regExOutPath=re.compile('([^\/]+)\/([^\/]+)\/([^\.]+).*\/out')
resultFile=open(resultFile,'w')
resultFile.write('\t'.join(['gene_name','model','target','w0','w1','w2','p0','p1','p2a','p2b','lnL'])+"\n")


#error
errorFile=open(outDir+'/logs/error','a')

os.chdir(pamlDir)# in paml dir for glob ! ! !
for fileName in glob.glob('*/*/*/out'):
	m=regExOutPath.match(fileName)
	if m:
		name=m.group(1)
		target=m.group(2)
		model=m.group(3)
		if model == 'M0' or model == 'M1' or model == 'b_free' or model == 'bsA' or model == 'bsA1':
			line=read_model(fileName,model,target,name)
			resultFile.write('\t'.join(line)+"\n")
resultFile.close()
errorFile.close()


