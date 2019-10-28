#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

def parse_dist(in_file):
    result = []
    distribution = {}
    name = []
    f = open(in_file)
    flag = 0
    for aline in f.readlines():
        line = aline.strip()
        dis = line.strip()
        dist = dis.split('\t')
        if len(dist) > 0:
            if len(dist) == 1:
                if dist[0].strip().find('coverage')==-1:
                    if flag == 0:
                        name.append(line)
                        flag = 1
                        t_name = line
                    else:
                        distribution[t_name] = 'null'
                        name.append(line)
                        flag = 1
                        t_name = line
            else:
                distri = []
                for i in range(0, len(dist)):
                    distri.append(dist[i].strip())
                distribution[t_name] = distri
                flag = 0
    result.append(name)
    result.append(distribution)
    f.close()
    return result
                
                







        





