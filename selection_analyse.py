#!/usr/bin/python

import argparse
import re 
import sys

version=1.0
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'

##parse argument
parser = argparse.ArgumentParser(description='allow',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-file', metavar='list of files', required=True, help="path of the 'out' files from ete3")
parser.add_argument('-results', metavar='results file', required=True, help="file to write all the results")

args=parser.parse_args()

#function
def read_model(path,model,target,name):

	file=open(path,'r')

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
		return [name,model,target,w0,w1,'NA','NA','NA','NA','NA',lnL]
	elif model == 'M1':
		return [name,model,'NA',w0,w1,'NA',p0,p1,'NA','NA',lnL]
	elif model == 'bsA':
		return [name,model,target,w0,w1,w2,p0,p1,p2a,p2b,lnL]
	elif model == 'bsA1':	
		return [name,model,target,w0,w1,w2,p0,p1,p2a,p2b,lnL]
	else:
		return [name,model,'NA',omega,'NA','NA','NA','NA','NA','NA',lnL]

rout=re.compile('([^\/]+)\/([^\/]+)_branch[^\/]*\/([^\.]+).*\/out')

results=open(args.results,'w')
list_out=open(args.file,'r')
results.write('\t'.join(['gene_name','model','target','w0','w1','w2','p0','p1','p2a','p2b','lnL'])+"\n")
for line in list_out.readlines():
	newline = line.rstrip('\n')
	if rout.match(newline):
		name=rout.match(newline).group(1)
		target=rout.match(newline).group(2)
		model=rout.match(newline).group(3)
		if model == 'M0' or model == 'M1' or model == 'b_free' or model == 'bsA' or model == 'bsA1':
			li=read_model(newline,model,target,name)
			results.write('\t'.join(li)+"\n")
			
list_out.close()
results.close()
