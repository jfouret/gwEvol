#!/usr/bin/python
import argparse
version=1.0
year=2016
author='Julien Fouret'
contact='julien@fouret.me'
##parse argument
parser = argparse.ArgumentParser(description='Perform statistics on alignments',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-repDir', metavar='/path', required=True, help="path of the report directory")
parser.add_argument('-fore', metavar='N', required=True, help="Maximum percentage of gap in the foreground ")
parser.add_argument('-back', metavar='N', required=True, help="Maximum percentage of gap in the background ")
parser.add_argument('-targets', metavar='spec1,spec2,...', required=True, help="Names of species part of the foreground")
args=parser.parse_args()

foreGap=0
foreN=0
foreSpec=0
backGap=0
backN=0
backSpec=0
length=None
foreList=args.targets.split(',')

alnDict=dict()

alnFileName=args.repDir+'/algn.fa'
resultFileName=args.repDir+'/algn_quality.tab'

alnFile=open(alnFileName,'r')
for line in alnFile.readlines():
	line=line.rstrip()
	if line[:1]=='>':
		key=line.lstrip('>')
		alnDict[key]=''
	else:
		alnDict[key]+=line
alnFile.close()
for spec,seq in alnDict.iteritems():
	if length==None:
		length=len(seq)
	if spec in foreList:
		foreN+=seq.count('N')
		foreGap+=seq.count('-')
		foreSpec+=1
	else:
		backN+=seq.count('N')
		backGap+=seq.count('-')
		backSpec+=1
validity='YES'
if ((backN+backGap)>(int(args.back)*length*backSpec/100)) or ((foreN+foreGap)>(int(args.fore)*length*foreSpec/100)):
	validity='NO'
with open(resultFileName,'w') as resultFile:
	header=['fore_gap','fore_N','fore_spec','fore_max','back_gap','back_N','back_spec','back_max','length','validity']
	line=[str(foreGap),str(foreN),str(foreSpec),args.fore,str(backGap),str(backN),str(backSpec),args.back,str(length),validity]
	resultFile.write("\t".join(header)+"\n"+"\t".join(line)+"\n")






