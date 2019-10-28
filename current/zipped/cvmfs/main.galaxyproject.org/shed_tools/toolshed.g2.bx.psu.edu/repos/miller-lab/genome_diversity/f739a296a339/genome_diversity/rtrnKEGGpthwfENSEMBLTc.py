#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       calclenchange.py
#       
#       Copyright 2011 Oscar Bedoya-Reina <oscar@niska.bx.psu.edu>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

import argparse,os,sys


def main():
	parser = argparse.ArgumentParser(description='Adds the fields KEGG gene codes and KEGG pathways to an input table of ENSEMBL transcript codes.')
	parser.add_argument('--loc_file',metavar='correlational database',type=str,help='correlational database')
	parser.add_argument('--species',metavar='species name',type=str,help='the species of interest in loc_file')
	parser.add_argument('--output',metavar='output TXT file',type=str,help='the output file with the table in txt format. The output will have two more fields: KEGG gene codes and KEGG pathways of each ENSEMBL code' )
	parser.add_argument('--posENSEMBLclmn',metavar='column number',type=int,help='the column with the ENSEMBLE transcript code')
	parser.add_argument('--input',metavar='input TXT file',type=str,help='the input file with the table in txt format')
	#~ 
	#~Open arguments 
	class C(object):
		pass
	fulargs=C()
	parser.parse_args(sys.argv[1:],namespace=fulargs)
	#test input vars
	inputf,loc_file,species,output,posENSEMBLclmn=fulargs.input,fulargs.loc_file,fulargs.species,fulargs.output,fulargs.posENSEMBLclmn
	posENSEMBLclmn-=1#correct pos
	#~ Get the extra variables
	crDB=[x.split() for x in open(loc_file).read().splitlines() if x.split()[0]==species][0]
	sppPrefx,dinput=crDB[0],crDB[1]#X should be replaced by the position in which the Conversion Dictionary File (CDF) is placed
	#make a dictionary of the input CDF
	dKEGGcPthws=dict([(x.split('\t')[0],'\t'.join(x.split('\t')[1:])) for x in open(dinput).read().splitlines() if x.strip()])
	#~ add the two new columns
	sall=[]
	#lENSEMBLTc=[x.split('\t') for x in open(inputf).read().splitlines() if x.strip()]
	lENSEMBLTc = []
	with open(inputf) as fh:
	    for line in fh:
	        if line.startswith('#'):
	            continue
	        lENSEMBLTc.append(line.rstrip('\r\n').split('\t'))
	nLines=len(lENSEMBLTc)
	cLines=0
	sall=[]#the output list for with additional fields
	#~ 
	while cLines<nLines:
		cLines+=1
		lENSEMBLTcKEGGgKEGGpth=lENSEMBLTc.pop(0)
		ENSEMBLTc=lENSEMBLTcKEGGgKEGGpth[posENSEMBLclmn]
		try:
			KEGGgKEGGpth=dKEGGcPthws[ENSEMBLTc]
		except:
			KEGGgKEGGpth='\t'.join(['U','N'])
		sall.append('\t'.join(['\t'.join(lENSEMBLTcKEGGgKEGGpth),KEGGgKEGGpth]))
	#~ 
	salef=open(output,'w')
	salef.write('\n'.join(sall))
	salef.close()
	return 0
	

if __name__ == '__main__':
	main()

