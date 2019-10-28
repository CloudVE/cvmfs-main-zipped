#!/usr/bin/env python

#  <command interpreter="python">
#    add_fst_column.py "$input" "$p1_input" "$p2_input" "$data_source.choice" "$data_source.min_value" "$retain" "$discard_fixed" "$biased" "$output"
#    #for $individual, $individual_col in zip($input.dataset.metadata.individual_names, $input.dataset.metadata.individual_columns)
#        #set $arg = '%s:%s' % ($individual_col, $individual)
#        "$arg"
#    #end for
#  </command>

import sys
import subprocess
from Population import Population

################################################################################

if len(sys.argv) < 12:
    print >> sys.stderr, "Usage"
    sys.exit(1)

input, p1_input, p2_input, genotypes, min_reads, min_qual, retain, discard_fixed, biased, output = sys.argv[1:11]
individual_metadata = sys.argv[11:]

p_total = Population()
p_total.from_tag_list(individual_metadata)

p1 = Population()
p1.from_population_file(p1_input)
if not p_total.is_superset(p1):
    print >> sys.stderr, 'There is an individual in population 1 that is not in the SNP table'
    sys.exit(1)

p2 = Population()
p2.from_population_file(p2_input)
if not p_total.is_superset(p2):
    print >> sys.stderr, 'There is an individual in population 2 that is not in the SNP table'
    sys.exit(1)

################################################################################

prog = 'Fst_column'

args = []
args.append(prog)
args.append(input)
args.append(genotypes)
args.append(min_reads)
args.append(min_qual)
args.append(retain)
args.append(discard_fixed)
args.append(biased)

columns = p1.column_list()
for column in columns:
    args.append('{0}:1'.format(column))

columns = p2.column_list()
for column in columns:
    args.append('{0}:2'.format(column))

fh = open(output, 'w')

#print "args:", ' '.join(args)
p = subprocess.Popen(args, bufsize=-1, stdin=None, stdout=fh, stderr=sys.stderr)
rc = p.wait()
fh.close()

sys.exit(0)

