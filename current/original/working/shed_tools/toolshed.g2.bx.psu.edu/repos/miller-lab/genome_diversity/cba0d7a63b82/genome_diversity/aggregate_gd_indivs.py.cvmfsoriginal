#!/usr/bin/env python

import sys
import subprocess
from Population import Population

################################################################################

if len(sys.argv) < 6:
    print >> sys.stderr, "Usage"
    sys.exit(1)

input, p1_input, output, input_type  = sys.argv[1:5]
individual_metadata = sys.argv[5:]

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

if input_type == 'gd_snp':
    args.append('1')
elif input_type == 'gd_genotype':
    args.append('0')
else:
    print >> sys.stderr, "unknown input type:", input_type
    sys.exit(1)

columns = p1.column_list()

for column in sorted(columns):
    if input_type == 'gd_genotype':
        column = str(int(column) - 2)
    args.append(column)

fh = open(output, 'w')

#print "args:", ' '.join(args)
p = subprocess.Popen(args, bufsize=-1, stdin=None, stdout=fh, stderr=sys.stderr)
rc = p.wait()
fh.close()

sys.exit(0)

