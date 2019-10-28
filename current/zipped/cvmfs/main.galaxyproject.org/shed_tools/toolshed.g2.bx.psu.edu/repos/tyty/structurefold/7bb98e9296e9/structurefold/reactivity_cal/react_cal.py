#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO
import math
from parse_dis_react import *
from react_norm_function import *
import os
import random
import string


dist_file1 = sys.argv[1] #plus library
dist_file2 = sys.argv[2] #minus library
seq_file = sys.argv[3] #Reference library(genome/cDNA)
nt_spec = sys.argv[4] #only show reactivity for AC or ATCG
flag_in = sys.argv[5] # perform 2-8% normalization (1) or not (0)
threshold = sys.argv[6] #Threshold to cap the reactivities
output_file = sys.argv[7]


distri_p = parse_dist(dist_file1)
distri_m = parse_dist(dist_file2)
threshold = float(threshold)


syspathrs = os.getcwd()

h = file(syspathrs+"react.txt",'w')
flag_in = int(flag_in)

seqs = SeqIO.parse(open(seq_file),'fasta');
nt_s = set()
for i in range(len(nt_spec)):
    nt_s.add(nt_spec[i])

flag = 0
trans = []
distri_p = distri_p[1]
distri_m = distri_m[1]

#thres = int(threshold)


transcripts = {}
for seq in seqs:
    n = seq.id
    trans.append(n)
    transcripts[n] = seq.seq.tostring()
    

#print(distri_p)
        

for i in range(0, len(trans)):
    h.write(trans[i])
    h.write('\n')       
    for j in range(len(distri_p[trans[i]])):
        distri_p[trans[i]][j] = math.log((int(distri_p[trans[i]][j])+1),math.e)
    for j in range(len(distri_m[trans[i]])):
        distri_m[trans[i]][j] = math.log((int(distri_m[trans[i]][j])+1),math.e)       
    s_p = sum(distri_p[trans[i]])
    s_m = sum(distri_m[trans[i]])
    length = len(distri_p[trans[i]])
    if s_p!= 0 and s_m!= 0:
        r = []
        for j in range(0, len(distri_p[trans[i]])):
            f_p = (float(distri_p[trans[i]][j]))/float(s_p)*length
            f_m = (float(distri_m[trans[i]][j]))/float(s_m)*length
            raw_react = f_p-f_m
            r.append(max(0, raw_react))
                
    if s_p!= 0 and s_m!= 0:    
        for k in range(1,(len(r)-1)):
            if transcripts[trans[i]][k-1] in nt_s:
                h.write(str(float('%.3f'%r[k])))
                h.write('\t')
            else:
                h.write('NA')
                h.write('\t')
        k = k+1
        if transcripts[trans[i]][k-1] in nt_s:
            h.write(str(float('%.3f'%r[k])))
            h.write('\n')
        else:
            h.write('NA')
            h.write('\n')
            

h.close()

if flag_in:
    react_norm((syspathrs+"react.txt"),output_file, threshold)
else:
    h_o = file(output_file, 'w')
    f_i = open(syspathrs+"react.txt")
    for aline in f_i.readlines():
        h_o.write(aline.strip())
        h_o.write('\n')
os.system("rm -f "+syspathrs+"react.txt")

#os.system("rm -r "+syspathrs)
    
     
            
    
    
        





















        





