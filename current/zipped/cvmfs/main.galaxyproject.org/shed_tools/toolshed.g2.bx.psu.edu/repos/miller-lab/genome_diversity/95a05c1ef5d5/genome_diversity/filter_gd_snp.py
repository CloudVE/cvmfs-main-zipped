#!/usr/bin/env python

import sys
import subprocess
from Population import Population

################################################################################

def convert_non_negative_int(string_value):
    try:
        val = int(string_value)
    except:
        print >> sys.stderr, '"%s" is not an integer' % string_value
        sys.exit(1)

    if val < 0:
        print >> sys.stderr, '"%d" is negative' % val
        sys.exit(1)

    return val
    

def convert_percent(string_value):
    if string_value.endswith('%'):
        val = convert_non_negative_int(string_value[:-1])
        if val > 100:
            print >> sys.stderr, 'percentage: "%d" > 100' % val
            sys.exit(1)
        val = val * -1
    else:
        val = convert_non_negative_int(string_value)

    return str(val)

################################################################################

if len(sys.argv) < 9:
    print >> sys.stderr, "Usage"
    sys.exit(1)

input, p1_input, output, lo, hi, lo_ind, lo_ind_qual = sys.argv[1:8]
individual_metadata = sys.argv[8:]

p_total = Population()
p_total.from_tag_list(individual_metadata)

p1 = Population()
p1.from_population_file(p1_input)

if not p_total.is_superset(p1):
    print >> sys.stderr, 'There is an individual in the population that is not in the SNP table'
    sys.exit(1)

lo = convert_percent(lo)
hi = convert_percent(hi)

################################################################################

prog = 'filter_snps'

args = []
args.append(prog)
args.append(input)
args.append(lo)
args.append(hi)
args.append(lo_ind)
args.append(lo_ind_qual)

columns = p1.column_list()

for column in sorted(columns):
    args.append(column)

fh = open(output, 'w')

#print "args:", ' '.join(args)
p = subprocess.Popen(args, bufsize=-1, stdin=None, stdout=fh, stderr=sys.stderr)
rc = p.wait()
fh.close()

sys.exit(0)

