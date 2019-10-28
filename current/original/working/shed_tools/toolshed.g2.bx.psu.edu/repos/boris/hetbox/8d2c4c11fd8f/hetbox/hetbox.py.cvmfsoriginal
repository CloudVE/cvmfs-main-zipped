#!/usr/bin/env python

# Code by Boris Rebolledo-Jaramillo
# (boris-at-bx.psu.edu)
# Edited by Nick Stoler
# (nick-at-bx.psu.edu)
# New in this version:
# - Add in proper header line if not present


import os
import sys
import array
import numpy
from rpy2.robjects import Formula
from rpy2.robjects.packages import importr
from rpy2 import robjects

def fail(message):
  sys.stderr.write(message+'\n')
  sys.exit(1)

COLUMN_LABELS = ['SAMPLE', 'CHR',  'POS', 'A', 'C', 'G', 'T', 'CVRG',
                 'ALLELES', 'MAJOR', 'MINOR', 'MINOR.FREQ.PERC.'] #, 'STRAND.BIAS']

COLUMN_LABELS_STRANDED= ['SAMPLE', 'CHR',  'POS', '+A', '+C', '+G', '+T', '-A', '-C', '-G', '-T',
                         'CVRG','ALLELES', 'MAJOR', 'MINOR', 'MINOR.FREQ.PERC.']

args = sys.argv[1:]
if len(args) >= 1:
  infile = args[0]
else:
  fail('Error: No input filename provided (as argument 1).')
if len(args) >= 2:
  outfile = args[1]
else:
  fail('Error: No output filename provided (as argument 2).')
if len(args) >= 3:
  report = args[2]
else:
  report = ''

# Check input file
add_header = False
if not os.path.exists(infile):
  fail('Error: Input file '+infile+' could not be found.')
with open(infile, 'r') as lines:
  line = lines.readline()
  if not line:
    fail('Error: Input file seems to be empty')
  line = line.strip().lstrip('#') # rm whitespace, comment chars
  labels = line.split("\t")
  if 'SAMPLE' not in labels:
    sys.stderr.write("Error: Input file does not seem to have a proper header "
      +"line.\nAdding an artificial header..")
    add_header = True

r = robjects.r
base = importr('base')
utils = importr('utils')
stats = importr('stats')
rprint = robjects.globalenv.get("print")
graphics = importr('graphics')
grdevices = importr('grDevices')
grdevices.png(file=outfile, width=1024, height=768, type="cairo")

# Read file into a data frame
if add_header:
    # add header line manually if not present
    DATA = utils.read_delim(infile, header=False)
    labels = robjects.r.names(DATA)
    for i in range(len(labels)):
        try:
            labels[i] = COLUMN_LABELS[i]
        except IndexError, e:
            try:
                labels[i] = COLUMN_LABELS_EXTENDED[i]
            except:
                fail("Error in input file: Too many columns (does not match hardcoded "
                +"column labels).")
else:
  DATA = utils.read_delim(infile)
  # Remove comment from header, if present
  labels = robjects.r.names(DATA)
  if labels[0][0:2] == 'X.':
    labels[0] = labels[0][2:]

# Multiply minor allele frequencies by 100 to get percentage
#  .rx2() looks up a column by its label and returns it as a vector
#  .ro turns the returned object into one that can be operated on per-element
minor_freq = DATA.rx2('MINOR.FREQ.PERC.').ro * 100
samples    = DATA.rx2('SAMPLE')

# Formula() creates a Python object representing the R object returned by x ~ y
formula = Formula('minor_freq ~ samples')
# The "environment" in .getenvironment() is the entire R workspace in which the
# Formula object exists. The R workspace meaning all the defined variables.
# Here, the .getenvironment() method is being used to set some variables in the
# R workspace

formula.getenvironment()['minor_freq'] = minor_freq
formula.getenvironment()['samples']    = samples


r.par(oma=array.array('i', [0,0,0,0]))
r.par(mar=array.array('i', [10,4,4,2]))
ylimit = array.array('i',[-5,50])

# create boxplot - fill kwargs1 with the options for the boxplot function
kwargs1 = {'ylab':"Minor allele frequency (%)", 'col':"gray", 'xaxt':"n",
           'outpch':"*",'main':"Distribution of minor allele frequencies",
           'cex.lab':"1.5"}
p = graphics.boxplot(formula, axes=0,ylim=ylimit, lty=1,**kwargs1)

table = base.table(DATA.rx2('SAMPLE'))
graphics.text(0, -1, 'N:', font=2)
for i in range(1, base.length(table)[0]+1, 1):
    graphics.text(i, -1, table[i-1], font=2)

graphlabels = base.names(table)
kwargs3 = {'pos':"-2", 'las':"2", 'cex.axis':"1"}
graphics.axis(1, at=range(1, len(graphlabels)+1, 1),labels=graphlabels, **kwargs3)
graphics.axis(2,at=(range(0,60,10)),pos=0,font=2)
grdevices.dev_off()

if not report:
    sys.exit(0)


SAMPLES=[]
for i in range(len(table)):
    SAMPLES.append(base.names(table)[i])

def boxstats(data,sample):
    VALUES = [100*float(x.strip().split('\t')[11]) for x in list(open(data)) if x.strip().split('\t')[0]==sample]
    NoHET  = len(VALUES)
    MEDIAN = numpy.median(VALUES)
    MAD    = numpy.median([abs(i - MEDIAN) for i in VALUES]) # Median absolute distance (robust spread statistic)
    return [NoHET,MEDIAN, MAD]

boxreport = open(report, "w+")
boxreport.write("#sample\tNo.sites\tmedian.freq\tMAD.freq\n")

for sample in SAMPLES:
    ENTRY = [sample] + boxstats(infile,sample)
    boxreport.write ('%s\t%d\t%.1f\t%.1f\n' % tuple([ENTRY[i] for i in [0,1,2,3]]))
boxreport.close()




