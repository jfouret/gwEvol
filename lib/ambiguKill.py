#!/usr/bin/env python
import argparse
#add gitRepository
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien.fouret@fouret.me'
parser = argparse.ArgumentParser(description="""

#-----------------------|ambiguKill|------------------------#
      random(reproducible) removing of ambiguous bases
      with the asumption that main source of polymorphism at
      DNA level is consequence of the genetic code redundancy 
      and not amino-acid polomorphism... to clarify
-------------------------------------------------------------

""",epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-input', metavar='/path', required=True, help="""
input file""")
parser.add_argument('-output', metavar='/path', required=True, help="""
output file""")
parser.add_argument('-seed', metavar='N', required=False,default="1991", help="""
""")

parser.add_argument('-v',action='store_true', help="verbose")
args=parser.parse_args()

import random
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO

iupacDict={
	"K":["G","T"],
	"M":["A","C"],
	"S":["C","G"],
	"W":["A","T"],
	"R":["A","G"],
	"Y":["T","C"],
	"A":["A"],
	"T":["T"],
	"C":["C"],
	"G":["G"],
	"N":["N"],
	"-":["-"]
}

from Bio.Data.CodonTable import TranslationError

def kept_codon(codon):
	global iupacDict
	list_codon=[""]
	for letter in codon:
		newListCodon=list()
		for possible_letter in iupacDict[letter]:
			for item in list_codon:
				newListCodon.append(item+possible_letter)
		list_codon=newListCodon
	list_aa=list()
	for cod in list_codon:
		try:
			list_aa.append(str(Seq(cod, generic_dna).translate()))
		except TranslationError:
			list_aa.append('?')

	aa_occur=dict()
	max_occur=0
	for aa in set(list_aa):
		aa_occur[aa]=list_aa.count(aa)
		if aa_occur[aa]>max_occur:
			max_occur=aa_occur[aa]

	keep_index=[]

	for index in range(0,len(list_codon)):
		if aa_occur[list_aa[index]]==max_occur:
			keep_index.append(index)

	return([list_codon[x] for x in keep_index])

random.seed(int(args.seed))
records=list( SeqIO.parse(args.input, "fasta"))
for index_record in range(0,len(records)):
	seq=str(records[index_record].seq)
	records[index_record].seq=Seq(''.join([random.sample(kept_codon(seq[3*x:3*x+3].upper()),1)[0] for x in range(0,(len(seq)/3))]),generic_dna)

SeqIO.write(records, args.output, "fasta")
