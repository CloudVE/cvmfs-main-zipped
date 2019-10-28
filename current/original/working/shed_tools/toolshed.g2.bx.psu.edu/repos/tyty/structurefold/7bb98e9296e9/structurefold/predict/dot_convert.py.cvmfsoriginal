#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

dot_file = sys.argv[1]
result_file = sys.argv[2]

h = file(result_file, 'w')
f = open(dot_file)



for aline in f.readlines():
    line = aline.strip()
    if line.find('>')!=-1:
        id_line = line
        idt = id_line.split('>')
        ids = idt[1].strip()
    else:
        if line.find('(')!=-1:
            structure_line = line
            st = structure_line.split(' ')
            structure = st[0].strip()
            enert = st[1].strip()
            if len(enert)>1:
                enertt = enert.split('(')
                enertt = enertt[1].strip()
            else:
                enertt = st[2].strip()
            enerttt = enertt.split(')')
            ener = enerttt[0].strip()
            h.write('>ENERGY = '+ener+'  '+ids+'\n')
            h.write(seq+'\n')
            h.write(structure+'\n')
        else:
            seq = line


    


f.close()
h.close()

