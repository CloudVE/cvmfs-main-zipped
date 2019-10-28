#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from read_file import *
from Bio import SeqIO

unmap_file = sys.argv[1]
reads_file = sys.argv[2]
result_file = sys.argv[3]
tp = sys.argv[4]


unmap = read_t_file(unmap_file);

h = file(result_file, 'w')

reads = SeqIO.parse(reads_file,tp)
um = set()
for i in range(0, len(unmap)):
    id_r = unmap[i][0]
    um.add(id_r)

for read in reads:
    if read.id in um:
        h.write('>')
        h.write(read.id)
        h.write('\n')
        h.write(read.seq.tostring())
        h.write('\n')
    


h.close()




