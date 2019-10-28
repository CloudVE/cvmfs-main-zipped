#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys



def read_sp_file(in_file):
    f = open(in_file);
    result = [];
    for aline in f.readlines():
        temp = [];
        tline = aline.strip();
        tl = tline.split(' ');
        for i in range(0, len(tl)):
            if len(tl[i].strip())>0:
                temp.append(tl[i].strip());
        result.append(temp);
    f.close();
    return result;


