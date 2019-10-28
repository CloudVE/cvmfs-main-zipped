#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from read_file import *
from Bio import SeqIO

map_file = sys.argv[1]
result_file = sys.argv[2]


#reads = read_t_file(read_file);

f = open(map_file);
h = file(result_file, 'w')

for aline in f.readlines():
    tline = aline.strip();
    tl = tline.split('\t');
    if len(tl)>4:
        if int(tl[1].strip())== 0:
            h.write(tline)
            h.write('\n')


f.close();
h.close()




