#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from read_file import *
from read_s_file import *
import random
import string

type_input = sys.argv[1]
seq_file = sys.argv[2]
ref_file = sys.argv[3]
shift = sys.argv[4]
length = sys.argv[5]
t_end = sys.argv[6]
map_type = sys.argv[7]
output_file = sys.argv[8]


if map_type!="default":
    s = ""
    sm = ""
    s = s+"-v "+sys.argv[9]
    sm = sm+"-v "+sys.argv[9]
    sm = sm+" -5 "+sys.argv[10]
    sm = sm+" -3 "+sys.argv[11]
    s = s+" -k "+sys.argv[12]
    sm = sm+" -k "+sys.argv[12]
    if sys.argv[13]:
        s = s+" -a"
        sm = sm+" -a"
    if int(sys.argv[14])>=1:
        s = s+" -m "+sys.argv[14]
        sm = sm+" -m "+sys.argv[14]
    if sys.argv[15]:
        s = s+" --best --strata "
        sm = sm+" --best --strata "
    
else:
    s = "-v 3 -a --best --strata "
    sm = "-v 3 -a --best --strata "

ospath = os.path.realpath(sys.argv[0])
ost = ospath.split('/')
syspath = ""
for i in range(len(ost)-1):
    syspath = syspath+ost[i].strip()
    syspath = syspath+'/'

syspathrs = os.getcwd()
syspathrs = syspathrs+'/'

os.system("bowtie-build -f "+ref_file+" "+syspathrs+"ref > "+syspathrs+"log.txt")

os.system("cp "+seq_file+" "+syspathrs+"seq0.fa")

if type_input == "fasta":
    tp = 'fasta'
if type_input == "fastq":
    tp = 'fastq'

k = 0

if type_input == "fasta":
    os.system("bowtie "+sm+"-f "+syspathrs+"ref"+" "+syspathrs+"seq"+str(k)+".fa --quiet -S > "+syspathrs+"map"+str(k)+".sam")
if type_input == "fastq":
    os.system("bowtie "+sm+"-q "+syspathrs+"ref"+" "+syspathrs+"seq"+str(k)+".fa --quiet -S > "+syspathrs+"map"+str(k)+".sam")

while(True):
    os.system("samtools view -Sb -F 0xfff "+syspathrs+"map"+str(k)+".sam > "+syspathrs+"mapped"+str(k)+".bam 2>"+syspathrs+"log.txt") #get mapped reads
    os.system("samtools view -Sb -f 0x4 "+syspathrs+"map"+str(k)+".sam > "+syspathrs+"umapped"+str(k)+".bam 2>"+syspathrs+"log.txt") #get unmapped reads
    os.system("samtools view -Sb -f 0x10 "+syspathrs+"map"+str(k)+".sam > "+syspathrs+"rmapped"+str(k)+".bam 2>"+syspathrs+"log.txt") #get reversed mapped reads
    os.system("samtools merge -f "+syspathrs+"unmapped"+str(k)+".bam "+syspathrs+"umapped"+str(k)+".bam "+syspathrs+"rmapped"+str(k)+".bam") #get reversed mapped reads
    os.system("samtools view -h -o "+syspathrs+"unmapped"+str(k)+".sam "+syspathrs+"unmapped"+str(k)+".bam") #get reversed mapped reads
    if k>0:
        os.system("samtools view -h -o "+syspathrs+"mapped"+str(k)+".sam "+syspathrs+"mapped"+str(k)+".bam") #get reversed mapped reads
        os.system("cut -f 1 "+syspathrs+"unmapped"+str(k)+".sam > "+syspathrs+"unmapped"+str(k)+".txt")
        os.system("cut -f 1 "+syspathrs+"mapped"+str(k)+".sam > "+syspathrs+"mapped"+str(k)+".txt")
        os.system("python "+syspath+"remove_map.py "+syspathrs+"unmapped"+str(k)+".txt "+syspathrs+"mapped"+str(k)+".txt "+syspathrs+"runmapped"+str(k)+".txt")
        os.system("rm "+syspathrs+"mapped"+str(k)+".sam")
        os.system("rm "+syspathrs+"mapped"+str(k)+".txt")
        os.system("rm "+syspathrs+"unmapped"+str(k)+".txt")
    else:
        os.system("cut -f 1 "+syspathrs+"unmapped"+str(k)+".sam > "+syspathrs+"runmapped"+str(k)+".txt")
    
    os.system("rm "+syspathrs+"unmapped"+str(k)+".bam")
    os.system("rm "+syspathrs+"umapped"+str(k)+".bam")
    os.system("rm "+syspathrs+"rmapped"+str(k)+".bam")
    os.system("python "+syspath+"seq_track.py "+syspathrs+"runmapped"+str(k)+".txt "+syspathrs+"seq"+str(k)+".fa "+syspathrs+"unmap_seq"+str(k)+".fa "+tp) #get unmapped sequence
    os.system("python "+syspath+"truncate.py "+syspathrs+"unmap_seq"+str(k)+".fa "+shift+" "+syspathrs+"seq"+str(k+1)+".fa "+length+" "+t_end) #truncate unmapped sequence
    os.system("rm "+syspathrs+"seq"+str(k)+".fa") #Remove sequences being mapped
    os.system("rm "+syspathrs+"map"+str(k)+".sam") #Remove mapping file
    os.system("rm "+syspathrs+"unmap_seq"+str(k)+".fa") #Remove unmapped sequnce
    os.system("rm "+syspathrs+"runmapped"+str(k)+".txt")
    os.system("rm "+syspathrs+"unmapped"+str(k)+".sam")
    
    os.system("wc -l "+syspathrs+"seq"+str(k+1)+".fa > "+syspathrs+"count"+str(k+1)+".txt")
    c = read_sp_file(syspathrs+"count"+str(k+1)+".txt")
    if c[0][0] == '0': #If no reads is in the sequence file, stop
        os.system("rm "+syspathrs+"count"+str(k+1)+".txt")
        os.system("rm "+syspathrs+"seq"+str(k+1)+".fa")
        break
    os.system("rm "+syspathrs+"count"+str(k+1)+".txt")
    k = k+1
    os.system("bowtie "+s+"-f "+syspathrs+"ref"+" "+syspathrs+"seq"+str(k)+".fa --quiet -S > "+syspathrs+"map"+str(k)+".sam")


ss = ""
for i in range(0,k+1):
    ss = ss+" "+syspathrs+"mapped"+str(i)+".bam"


os.system("samtools merge -f "+syspathrs+"combine.bam"+" "+ss)
os.system("samtools sort "+syspathrs+"combine.bam sorted")
os.system("samtools view -b -h sorted.bam > " + output_file)
#print("samtools merge mapped_all.bam"+ss)
os.system("rm "+syspathrs+"mapped*.bam")
os.system("rm "+syspathrs+"combine.bam")
os.system("rm "+syspathrs+"sorted.bam")
os.system("rm "+syspathrs+"ref*")
#os.system("rm -r "+syspathrs)


