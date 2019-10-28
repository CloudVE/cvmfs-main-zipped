#!/usr/bin/env python
#Make a plot of reactivity distribution

import sys
import os
import numpy as np
import matplotlib
from pylab import *
import math

#Convert the reactivities (Make NA to 0)
def convert_react(a):
    r = []
    for i in range(len(a)):
        if a[i]!='NA':
            r.append(float(a[i]))
        else:
            r.append(float(0))
    return r
        

#Make a plot of the distribution
def make_plot(ar,id_s,path):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 16}
    matplotlib.rc('font', **font)
    N = len(ar)
    a = convert_react(ar)
    w = 1
    ind = np.arange(N)

    fig = figure()
    fig, ax = subplots()
    ax.bar(ind+w, a, width = w, color = 'black',edgecolor = 'black')
    ax.set_ylabel('Final Structural Reactivity (FSR)')
    ax.set_xlabel('Nucleotide Number')

    
    mag = int(math.log(N,10))-1
    tail = 10**mag

    intervel = int(math.ceil(float(N)/tail/5))
    tl = []
    k = 0
    upmax = int(math.ceil(float(N)/intervel/tail)*intervel*tail)+1
    ax.set_xticks(np.arange(0,upmax,intervel*tail))
    ax.set_xticklabels(np.arange(0,upmax,intervel*tail))
    savefig(os.path.join(path, id_s+'.tif'))



    
    
    


