#!/usr/bin/env python
import argparse
#add gitRepository
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien.fouret@fouret.me'
parser = argparse.ArgumentParser(description="""

#-----------------------|gwEvol-paml|-----------------------#
     Genome-wide screening of coding sequences with PAML
-------------------------------------------------------------
This software serializes genome-wide analyses of molecular ev-
olution on HPC infrastructure. All coding sequences are expec-
ted to be produced using gwAlign-Unify. 

The ete3 software is needed to drive PAML analysis. 

All models available in ete3 can be runned. However default
arguments will only run models necessary for Branch and 
Branch-site positive selection analyses.

""",epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-aln_dir', metavar='/path', required=True, help="""
MSAs repository produced with gwAlign-Unify""")
parser.add_argument('-mark', metavar='ete3_mark', required=True, help="""
