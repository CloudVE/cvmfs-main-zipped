#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO
import math
from parse_dis_react import *

def cap(a,value):
    if a>=value:
        return value
    else:
        return a

def react_norm(react_file, result_file, capped_value):
    print("Normalizing.....")
    react1 = parse_dist(react_file)
    react = react1[1]
    h = file(result_file, 'w')

    capped = int(capped_value)

    all_react = []


    for t in react:
        if react[t]!='null':
            for i in range(len(react[t])):
                if react[t][i]!='NA':                   
                    all_react.append(float(react[t][i]))


    all_react.sort(reverse = True)


    eight = all_react[int(len(all_react)*0.02):int(len(all_react)*0.1)]
    meight = sum(eight)/len(eight)

    for t in react:
        h.write(t)
        h.write('\n')
        if react[t]!='null':
            for i in range((len(react[t])-1)):
                if react[t][i]!='NA':
                    h.write(str(float('%.3f'%cap((float(react[t][i])/meight),capped))))
                else:
                    h.write('NA')
                h.write('\t')
            if react[t][i+1]!='NA':
                h.write(str(float('%.3f'%cap((float(react[t][i+1])/meight),capped))))
            else:
                h.write('NA')
            h.write('\n')

    h.close()
        





















        





