#!/usr/bin/env python

import sys
import subprocess
from Population import Population

################################################################################

if len(sys.argv) < 9:
    print >> sys.stderr, "Usage"
    sys.exit(1)

#input, p1_input, output, lo, hi, lo_ind, lo_ind_qual = sys.argv[1:8]
#individual_metadata = sys.argv[8:]
input, p1_input, output,  = sys.argv[1:4]
individual_metadata = sys.argv[4:]

p_total = Population()
p_total.from_tag_list(individual_metadata)

p1 = Population()
p1.from_population_file(p1_input)

if not p_total.is_superset(p1):
    print >> sys.stderr, 'There is an individual in the population that is not in the SNP table'
    sys.exit(1)

################################################################################

prog = 'aggregate'

args = []
args.append(prog)
args.append(input)

columns = p1.column_list()

for column in sorted(columns):
    args.append(column)

fh = open(output, 'w')

#print "args:", ' '.join(args)
p = subprocess.Popen(args, bufsize=-1, stdin=None, stdout=fh, stderr=sys.stderr)
rc = p.wait()
fh.close()

sys.exit(0)

