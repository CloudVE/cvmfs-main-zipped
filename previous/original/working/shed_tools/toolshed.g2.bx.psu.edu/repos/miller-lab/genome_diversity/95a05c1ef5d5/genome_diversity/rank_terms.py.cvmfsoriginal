#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       GOFisher.py
#       
#       Copyright 2013 Oscar Reina <oscar@niska.bx.psu.edu>
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

import argparse
import os
import sys
from fisher import pvalue
from decimal import Decimal,getcontext

def rtrnGOcENSEMBLc(inExtnddfile,columnENSEMBLTExtndd,columnGOExtndd):
	"""
	"""
	dGOTENSEMBLT={}
	for eachl in open(inExtnddfile,'r'):
		if eachl.strip():
			ENSEMBLT=eachl.splitlines()[0].split('\t')[columnENSEMBLTExtndd]
			GOTs=set(eachl.splitlines()[0].split('\t')[columnGOExtndd].split('.'))
			GOTs=GOTs.difference(set(['','U','N']))
			for GOT in GOTs:
				try:
					dGOTENSEMBLT[GOT].add(ENSEMBLT)
				except:
					dGOTENSEMBLT[GOT]=set([ENSEMBLT])
	#~ 
	##dGOTENSEMBLT.pop('')
	ENSEMBLTGinGO=set.union(*dGOTENSEMBLT.values())
	return dGOTENSEMBLT,ENSEMBLTGinGO

def rtrnENSEMBLcSAPs(inSAPsfile,columnENSEMBLT,ENSEMBLTGinGO):
	"""
	returns a set of the ENSEMBLT codes present in the input list and
	in the GO file
	"""
	sENSEMBLTSAPsinGO=set()
	for eachl in open(inSAPsfile,'r'):
		ENSEMBLT=eachl.splitlines()[0].split('\t')[columnENSEMBLT]
		if ENSEMBLT in ENSEMBLTGinGO:
			sENSEMBLTSAPsinGO.add(ENSEMBLT)
	return sENSEMBLTSAPsinGO

def rtrnCounts(dGOTENSEMBLT,ENSEMBLTGinGO,sENSEMBLTSAPsinGO):
	"""
	returns a list of the ENSEMBLT codes present in the input list and
	in the GO file. The terms in this list are: 'Go Term','# Genes in
	the GO Term','# Genes in the list and in the GO Term','Enrichement
	of the GO Term for genes in the input list','Genes in the input list
	present in the GO term'
	"""
	getcontext().prec=2#set 2 decimal places
	SAPs_all=len(sENSEMBLTSAPsinGO)
	NoSAPs_all=len(ENSEMBLTGinGO)-SAPs_all
	#~ 
	lp=len(dGOTENSEMBLT)
	cnt=0
	#~ 
	ltfreqs=[]
	for echGOT in dGOTENSEMBLT:
		cnt+=1
		##print 'Running "%s", %s out of %s'%(echGOT,cnt,lp)
		CntGO_All=len(dGOTENSEMBLT[echGOT])
		SAPs_GO=len(dGOTENSEMBLT[echGOT].intersection(sENSEMBLTSAPsinGO))
		NoSAPs_GO=CntGO_All-SAPs_GO
		pval=pvalue(SAPs_GO,NoSAPs_GO,SAPs_all,NoSAPs_all)
		#~ outl.append('\t'.join([str(x) for x in(str(pval.two_tail),CntGO_All,SAPs_GO,'.'.join(sorted(dGOTENSEMBLT[echGOT].intersection(sENSEMBLTSAPsinGO))),echGOT)]))
		ltfreqs.append([(SAPs_GO/Decimal(CntGO_All)),SAPs_GO,Decimal(str(pval.two_tail))*1,echGOT])
	#~ 
	ltfreqs.sort()
	ltfreqs.reverse()
	outl=[]
	cper,crank=Decimal('2'),0
	#~ 
	for perc,cnt_go,pval,goTerm in ltfreqs:
		if perc<cper:
			crank+=1
			cper=perc
		outl.append('\t'.join([str(cnt_go),str(perc),str(crank),str(pval),goTerm]))
	#~ 
	return outl
	

def main():
	#~ 
	parser = argparse.ArgumentParser(description='Returns the count of genes in GO categories and their statistical overrrepresentation, from a list of genes and an extended file (i.e. plane text with ENSEMBLT and GO terms).')
	parser.add_argument('--input',metavar='input TXT file',type=str,help='the input file with the table in txt format.',required=True)
	parser.add_argument('--inExtnddfile',metavar='input TXT file',type=str,help='the input file with the extended table in txt format.',required=True)
	parser.add_argument('--output',metavar='output TXT file',type=str,help='the output file with the table in txt format.',required=True)
	parser.add_argument('--columnENSEMBLT',metavar='column number',type=int,help='column with the ENSEMBL transcript code in the input file.',required=True)
	parser.add_argument('--columnENSEMBLTExtndd',metavar='column number',type=int,help='column with the ENSEMBL transcript code in the extended file.',required=True)
	parser.add_argument('--columnGOExtndd',metavar='column number',type=int,help='column with the GO terms in the extended file.',required=True)

	args = parser.parse_args()

	inSAPsfile = args.input
	inExtnddfile = args.inExtnddfile
	saleGOPCount = args.output
	columnENSEMBLT = args.columnENSEMBLT
	columnENSEMBLTExtndd = args.columnENSEMBLTExtndd
	columnGOExtndd = args.columnGOExtndd

	#~ 
	dGOTENSEMBLT,ENSEMBLTGinGO=rtrnGOcENSEMBLc(inExtnddfile,columnENSEMBLTExtndd,columnGOExtndd)
	sENSEMBLTSAPsinGO=rtrnENSEMBLcSAPs(inSAPsfile,columnENSEMBLT,ENSEMBLTGinGO)
	outl=rtrnCounts(dGOTENSEMBLT,ENSEMBLTGinGO,sENSEMBLTSAPsinGO)
	#~ 
	saleGOPCount=open(saleGOPCount,'w')
	saleGOPCount.write('\n'.join(outl))
	saleGOPCount.close()
	#~ 
	return 0

if __name__ == '__main__':
	main()

