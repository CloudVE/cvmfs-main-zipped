#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys



def read_t_file(in_file):
    f = open(in_file);
    result = [];
    for aline in f.readlines():
        temp = [];
        tline = aline.strip();
        tl = tline.split('\t');
        for i in range(0, len(tl)):
            temp.append(tl[i].strip());
        result.append(temp);
    f.close();
    return result;


