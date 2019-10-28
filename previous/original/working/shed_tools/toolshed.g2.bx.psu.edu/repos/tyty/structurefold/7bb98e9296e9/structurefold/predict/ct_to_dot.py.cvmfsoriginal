#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import shlex
import os
import subprocess
from read_file import *

ct_file = sys.argv[1]
path = sys.argv[2]
id_s = sys.argv[3]
result_file = sys.argv[4]

h = file(result_file, 'w')
os.system('grep "'+id_s+'" '+ct_file+' |wc -l > '+path+'/count.txt')
count = read_t_file(path+'/count.txt')
comm = ''
for i in range(int(count[0][0])):
    command = shlex.split('ct2dot %s %s %s' % (ct_file, str(i+1), os.path.join(path, 'db_file_%s.dbnn' % str(i+1))))
    subprocess.call(command)
    comm = comm +' '+path+'/db_file_'+str(i+1)+'.dbnn' 



os.system('cat'+comm+' > '+result_file)
for i in range(int(count[0][0])):
    command = shlex.split('rm %s' % (os.path.join(path, 'db_file_%s.dbnn' % str(i+1))))
    subprocess.call(command)
command = shlex.split('rm %s' % (os.path.join(path, 'count.txt')))
subprocess.call(command)
    


h.close()

