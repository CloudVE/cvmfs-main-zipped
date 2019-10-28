#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
import os
from read_file import *
import random
import string

fasta_file = sys.argv[1]
map_file = sys.argv[2]
result_file = sys.argv[3]

syspathrs = os.getcwd()

os.system("samtools view -F 0xfff "+map_file+"|cut -f 3,4 > "+syspathrs+"map_info.txt") 

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta');
length_seq = {};
for seq in fasta_sequences:
        nuc = seq.id;
        length_seq[nuc] = len(seq.seq.tostring());



mapping = {}
transcripts = []

f = open(syspathrs+"map_info.txt");
for aline in f.readlines():
    tline = aline.strip();
    tl = tline.split('\t');
    if tl[0].strip() not in transcripts:
        transcripts.append(tl[0].strip());
        mapping[tl[0].strip()] = [];

    mapping[tl[0].strip()].append(tl[1].strip());

distribution = {};
coverage = {};
for transcript in length_seq:
    distribution[transcript] = [];
    for i in range(0, length_seq[transcript]):
        distribution[transcript].append(0);
    sum_count = float(0);
    if transcript in mapping:
        for j in range(0, len(mapping[transcript])):
            index = mapping[transcript][j];
            #count = reads[mapping[transcript][j][0]];
            sum_count = sum_count + 1;
            distribution[transcript][int(index)-1] = distribution[transcript][int(index)-1] + 1;
            coverage[transcript] = float(sum_count)/float(length_seq[transcript]);
    else:
        coverage[transcript] = 0

        
        
    

h = file(result_file, 'w')
for transcript in length_seq:
    h.write(transcript);
    h.write('\n')
    for i in range(0, length_seq[transcript]):
        h.write(str(distribution[transcript][i]))
        h.write('\t')
    h.write('\n')
    h.write('\n')

#os.system("rm -r "+syspathrs)

    

f.close();
h.close()




