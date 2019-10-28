#RNA structure prediction & Output and illustrate reactivities

import sys
import shlex
import subprocess
import tarfile
from parse_dis_pac import *
from read_file import *
from Bio import SeqIO
import os
from rtts_plot import *
import random
import string


id_file = sys.argv[1]
seq_file = sys.argv[2]
predict_type = sys.argv[3]
temperature = sys.argv[4]
predict_program = sys.argv[5]
output_html = sys.argv[6]
output_directory = sys.argv[7]



flag = False
if predict_type!='silico': #input reactivity file if provided
    if predict_program == 'rs':
        react_file = sys.argv[8]
        slope = sys.argv[9]
        intercept = sys.argv[10]
    else:
        react_file = sys.argv[8]
        thres_h = sys.argv[9]
        thres_h = float(thres_h)
        thres_l = sys.argv[10]
        thres_l = float(thres_l)
        gqs = sys.argv[11]
        gqs = int(gqs)
        
    react = parse_dist(react_file)
    react = react[1]
    flag = True
else:
    if predict_program!='rs':
        gqs = sys.argv[8]
        gqs = int(gqs)


ospath = os.path.realpath(sys.argv[0])
ost = ospath.split('/')
syspathpt = ""
for i in range(len(ost)-1):
    syspathpt = syspathpt+ost[i].strip()
    syspathpt = syspathpt+'/'


syspath = os.getcwd()

ids = read_t_file(id_file)
sequences = SeqIO.parse(seq_file, 'fasta')


seqs = {}
for seq in sequences:
    seqs[seq.id] = seq.seq.tostring()

if len(ids)>100: #setup a limit of the number of sequence to be predicted
    print("Number of sequences exceeds limitation!")
    sys.exit(0)
    

#predict RNA structures

os.mkdir(output_directory)
flag3 = 0

id_predicted = set()
for i in range(len(ids)):
    flag2 = 0
    id_s = ids[i][0]
    #print(id_s)
    #Put RNA sequence and reactivities into files
    if id_s in seqs:
        fh = file(os.path.join(syspath,"temp.txt"), 'w')        
        fh.write('>'+id_s)
        fh.write('\n')
        fh.write(seqs[id_s])
        fh.close()
        if not flag:
            if predict_program == 'rs':
                command = shlex.split('Fold %s -T %s %s' % (os.path.join(syspath, 'temp.txt'), temperature, os.path.join(output_directory, '%s.ct' % id_s)))
                subprocess.call(command)
                command = shlex.split('python %s %s %s %s %s' % (os.path.join(syspathpt, 'ct_to_dot.py'), os.path.join(output_directory, '%s.ct' % id_s), output_directory, id_s, os.path.join(output_directory, '%s.dbn' % id_s)))
                subprocess.call(command)               
            else:
                if gqs:
                    os.system('RNAfold < '+syspath+'/temp.txt -T '+str(float(temperature)-273.15)+' --noconv -g > '+output_directory+'/'+id_s+'.dbnb')
                    
                else:
                    os.system('RNAfold < '+syspath+'/temp.txt -T '+str(float(temperature)-273.15)+' --noconv --noPS > '+output_directory+'/'+id_s+'.dbnb')
                command = shlex.split('python %s %s %s' % (os.path.join(syspathpt, 'dot_convert.py'), os.path.join(output_directory, '%s.dbnb' % id_s), os.path.join(output_directory, '%s.dbn' % id_s)))
                subprocess.call(command)
                if not gqs:
                    command = shlex.split('dot2ct %s %s' % (os.path.join(output_directory, '%s.dbn' % id_s), os.path.join(output_directory, '%s.ct' % id_s)))
                else:
                    command = shlex.split('mv -f %s %s' % (os.path.join(syspath, '%s_ss.ps' % id_s), os.path.join(output_directory, '%s.ps' % id_s)))
                subprocess.call(command)
                command = shlex.split('rm %s' % (os.path.join(output_directory, '%s.dbnb' % id_s)))
                subprocess.call(command)
        else:
            if id_s in react:
                fh = file(os.path.join(syspath, "constraint.txt"), 'w')
                make_plot(react[id_s], id_s, output_directory) #make a plot of the distribution of the reactivites of the input RNA
                if predict_program == 'rs': 
                    for j in range(0, (len(react[id_s]))):
                        if react[id_s][j]!='NA':
                            fh.write(str(j+1))
                            fh.write('\t')
                            fh.write(str(react[id_s][j]))
                            fh.write('\n')
                    fh.close()
                    command = shlex.split("Fold %s -sh %s -si %s -sm %s -T %s %s" % (os.path.join(syspath, "temp.txt"), 
                                                                 os.path.join(syspath, "constraint.txt"), intercept, slope, temperature, 
                                                                 os.path.join(output_directory, "%s.ct" % id_s)))
                    subprocess.call(command)
                    command = shlex.split('python %s %s %s %s %s' % (os.path.join(syspathpt, 'ct_to_dot.py'), os.path.join(output_directory, '%s.ct' % id_s), output_directory, id_s, os.path.join(output_directory, '%s.dbn' % id_s)))
                    subprocess.call(command)
                else:
                    fh.write('>'+id_s)
                    fh.write('\n')
                    fh.write(seqs[id_s])
                    fh.write('\n')
                    for j in range(0, (len(react[id_s]))):
                        if react[id_s][j]!='NA':
                            re = float(react[id_s][j])
                            if re>thres_h:
                                fh.write('x')
                            else:
                                if re<thres_l:
                                    fh.write('|')
                                else:
                                    fh.write('.')
                        else:
                            fh.write('.')
                    fh.write('.')
                    fh.close()
                    if gqs:
                        os.system('RNAfold < '+syspath+'/constraint.txt -T '+str(float(temperature)-273.15)+' -C --noconv -g > '+output_directory+'/'+id_s+'.dbnb')
                        
                    else:
                        os.system('RNAfold < '+syspath+'/constraint.txt -T '+str(float(temperature)-273.15)+' -C --noconv --noPS > '+output_directory+'/'+id_s+'.dbnb')
                    command = shlex.split('python %s %s %s' % (os.path.join(syspathpt, 'dot_convert.py'), os.path.join(output_directory, '%s.dbnb' % id_s), os.path.join(output_directory, '%s.dbn' % id_s)))
                    subprocess.call(command)
                    if not gqs:
                        command = shlex.split('dot2ct %s %s' % (os.path.join(output_directory, '%s.dbn' % id_s), os.path.join(output_directory, '%s.ct' % id_s)))
                    else:
                        command = shlex.split('mv -f %s %s' % (os.path.join(syspath, '%s_ss.ps' % id_s), os.path.join(output_directory, '%s.ps' % id_s)))
                    subprocess.call(command)
                    command = shlex.split('rm %s' % (os.path.join(output_directory, '%s.dbnb' % id_s)))
                    subprocess.call(command)                  

            else:
                print(id_s+" not in the data of react!")
                flag2 = 1
        if flag2 == 0:
            if predict_program == 'rs':
                command = shlex.split('draw %s.ct %s.ps' % (os.path.join(output_directory, id_s), os.path.join(output_directory, id_s)))
                subprocess.call(command)
                command = shlex.split('rm %s' % (os.path.join(output_directory, '%s.ct' % id_s)))
                subprocess.call(command)
            else:
                if not gqs:
                    command = shlex.split('draw %s.ct %s.ps' % (os.path.join(output_directory, id_s), os.path.join(output_directory, id_s)))
                    subprocess.call(command)
                    command = shlex.split('rm %s' % (os.path.join(output_directory, '%s.ct' % id_s)))
                    subprocess.call(command)
            flag3 = 1
        id_predicted.add(id_s)
    else:
        print(id_s+" not in the data of sequences!")

#Remove the unnecessary files
if flag3 == 1:

    tarball = tarfile.open(os.path.join(output_directory,'prediction_results.tar'), 'w:')
    for filename in os.listdir(output_directory):
        filepath = os.path.join(output_directory, filename)
        print filepath
        tarball.add(filepath, arcname=filename)
    #print os.listdir(syspath)
    #print os.listdir(output_directory)
    # tarball.add('%s.tif' % os.path.join(syspath, id_s), arcname='%s.tif' % id_s)
    tarball.close()
 
    h = open(output_html, 'wb' )
    h.write('<html><head><title><h1>Results of RNA structure prediction</h1></title></head><body>\n')

    h.write('<p>\n')
    h.write('Click <a href="%s">here</a> to download the compressed file containing all prediction results.\n' % (('prediction_results.tar')))
    #h.write('<\p>\n')
    h.write('<hr>\n')

    
    for id_p in id_predicted:
        h.write('<h4>'+id_p+'</h4><p><ul>\n')
        h.write('<li><a href="%s">%s</a></li>\n' % (('%s.dbn' % id_p), ('%s.dbn' % id_p)))
        h.write('<li><a href="%s">%s</a></li>\n' % (('%s.ps' % id_p), ('%s.ps' % id_p)))
        if flag:
            h.write('<li><a href="%s">%s</a></li>\n' % (('%s.tif' % id_p), ('%s.tif' % id_p)))
        h.write( '</ul></p>\n' )
        h.write('<hr>\n')
    h.write( '</body></html>\n' )
    h.close()
                
                
        
    

