#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       calcfreq.py
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
from decimal import Decimal,getcontext
from LocationFile import LocationFile
from fisher import pvalue

#method to rank the the pthways by mut. freq.
def rankd(ltfreqs):
	ordvals=sorted(ltfreqs)#sort and reverse freqs.
	#~ 
	outrnk=[]
	tmpFreq0,tmpCount,pval,tmpPthw=ordvals.pop()#the highest possible value
	crank=1
	outrnk.append('\t'.join([str(tmpCount),str(tmpFreq0),str(crank),str(pval),tmpPthw]))
	totalnvals=len(ordvals)
	cnt=0
	while totalnvals>cnt:
		cnt+=1
		tmpFreq,tmpCount,pval,tmpPthw=ordvals.pop()
		if tmpFreq!=tmpFreq0:
			crank=len(outrnk)+1
			tmpFreq0=tmpFreq
		outrnk.append('\t'.join([str(tmpCount),str(tmpFreq),str(crank),str(pval),tmpPthw]))
	return outrnk
		

def main():
	parser = argparse.ArgumentParser(description='Obtain KEGG images from a list of genes.')
	parser.add_argument('--input',metavar='input TXT file',type=str,help='the input file with the table in txt format')
	parser.add_argument('--output',metavar='output TXT file',type=str,help='the output file with the table in txt format. Column 1 is the count of genes in the list, Column 2 is the percentage of the pathway genes present on the list. Column 3 is the rank based on column 2')
	parser.add_argument('--posKEGGclmn',metavar='column number',type=int,help='the column with the KEGG pathway code/name')
	parser.add_argument('--KEGGgeneposcolmn',metavar='column number',type=int,help='column with the KEGG gene code')
	parser.add_argument('--loc_file',metavar='location file',type=str,help='location file')
	parser.add_argument('--species',metavar='species',type=str,help='species')
	#~Open arguments 
	class C(object):
		pass
	fulargs=C()
	parser.parse_args(sys.argv[1:],namespace=fulargs)
	#test input vars
	inputf,outputf,posKEGGclmn,Kgeneposcolmn=fulargs.input,fulargs.output,fulargs.posKEGGclmn,fulargs.KEGGgeneposcolmn
	locf,species=fulargs.loc_file,fulargs.species
	#make a dictionary of valid genes
	posKEGGclmn-=1
	Kgeneposcolmn-=1
	dKEGGcPthws=dict([(x.split('\t')[Kgeneposcolmn],set(x.split('\t')[posKEGGclmn].split('.'))) for x in open(inputf).read().splitlines()[1:] if x.split('\t')[posKEGGclmn] not in set(['U','N'])])
	for u in ['U','N']:
		try:
			a=dKEGGcPthws.pop(u)
		except:
			pass
	getcontext().prec=2#set 2 decimal places
	sdGenes=set([x for x in dKEGGcPthws.keys() if x.find('.')>-1])
	while True:#to correct names with more than one gene
		try:
			mgenes=sdGenes.pop()
			pthwsAssotd=dKEGGcPthws.pop(mgenes)
			mgenes=mgenes.split('.')
			for eachg in mgenes:
				dKEGGcPthws[eachg]=pthwsAssotd
		except:
			break
	#~ Count genes

	location_file = LocationFile(locf)
	prefix, kxml_dir_path, dict_file = location_file.get_values(species)
	dPthContsTotls = {}
	try:
	    with open(dict_file) as fh:
	        for line in fh:
	            line = line.rstrip('\r\n')
	            value, key = line.split('\t')
	            dPthContsTotls[key] = int(value)
	except IOError, err:
	    print >> sys.stderr, 'Error opening dict file {0}: {1}'.format(dict_file, err.strerror)
	    sys.exit(1)
	
	dPthContsTmp=dict([(x,0) for x in dPthContsTotls.keys()])#create a list of genes
	sdGenes=set(dKEGGcPthws.keys())#list of all genes
	cntGens=0
	ltGens=len(sdGenes)
	while cntGens<ltGens:
		cGen=sdGenes.pop()
		sKEGGcPthws=dKEGGcPthws.pop(cGen)
		for eachP in sKEGGcPthws:
			if eachP!='N':
				if eachP in dPthContsTmp:
					dPthContsTmp[eachP]+=1
				else:
					print >> sys.stderr, "Error: pathway not found in database: '{0}'".format(eachP)
					sys.exit(1)
		cntGens+=1
	#~ Calculate Freqs.
	ltfreqs=[]
	cntAllKEGGinGnm=sum(dPthContsTotls.values())
	cntListKEGGinGnm=sum(dPthContsTmp.values())
	cntAllKEGGNOTinGnm=cntAllKEGGinGnm-cntListKEGGinGnm
	for pthw in dPthContsTotls:
		cntAllKEGGinpthw=Decimal(dPthContsTotls[pthw])
		try:
			cntListKEGGinpthw=Decimal(dPthContsTmp[pthw])
		except:
			cntListKEGGinpthw=Decimal('0')
		cntAllKEGGNOTinpthw=cntAllKEGGinpthw-cntListKEGGinpthw
		pval=pvalue(cntListKEGGinpthw,cntAllKEGGNOTinpthw,cntListKEGGinGnm,cntAllKEGGNOTinGnm)

		ltfreqs.append([(cntListKEGGinpthw/cntAllKEGGinpthw),cntListKEGGinpthw,Decimal(str(pval.two_tail))*1,pthw])
	#~ ltfreqs=[((Decimal(dPthContsTmp[x])/Decimal(dPthContsTotls[x])),Decimal(dPthContsTmp[x]),x) for x in dPthContsTotls]
	tabllfreqs='\n'.join(rankd(ltfreqs))
	salef=open(outputf,'w')
	salef.write(tabllfreqs)
	salef.close()
	return 0
	

if __name__ == '__main__':
	main()
