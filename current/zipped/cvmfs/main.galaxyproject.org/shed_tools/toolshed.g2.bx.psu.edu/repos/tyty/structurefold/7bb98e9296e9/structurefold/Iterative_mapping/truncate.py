#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]
shift_in = sys.argv[2]
result_file = sys.argv[3]
length = sys.argv[4]
t_end = sys.argv[5]

shift = int(shift_in)
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta');
h = file(result_file,'w')
for seq in fasta_sequences:
        nuc = seq.id;
        sequence = seq.seq.tostring();
        if (len(sequence)-shift)>=int(length):
                h.write('>'+nuc)
                h.write('\n')
                if t_end == 'three_end':
                        h.write(sequence[0:(len(sequence)-shift)])
                if t_end == 'five_end':
                        h.write(sequence[(shift):(len(sequence))])
                h.write('\n')




h.close()




