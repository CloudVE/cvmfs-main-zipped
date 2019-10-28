#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from read_file import *


unmap_file = sys.argv[1]
map_file = sys.argv[2]
result_file = sys.argv[3]


unmap = read_t_file(unmap_file)
mapped = read_t_file(map_file)
h = file(result_file, 'w')

maps = set()
for i in range(len(mapped)):
    maps.add(mapped[i][0])


for i in range(len(unmap)):
    name = unmap[i][0]
    if name not in maps:
        h.write(name)
        h.write('\n')


h.close()
