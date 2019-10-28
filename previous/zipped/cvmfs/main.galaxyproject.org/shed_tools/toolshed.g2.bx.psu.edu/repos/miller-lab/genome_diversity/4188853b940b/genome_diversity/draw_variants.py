#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

if len(sys.argv) != 10:
    gd_util.die('Usage')

snp_input, indel_input, coverage_input, annotation_input, indiv_input, ref_name, min_coverage, output, ind_arg = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(indiv_input)
if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population individuals that is not in the SNP table')

################################################################################

prog = 'mk_Ji'

args = [ prog ]
args.append(snp_input)
args.append(indel_input)
args.append(coverage_input)
args.append(annotation_input)
args.append(min_coverage)
args.append(ref_name)

for tag in p1.tag_list():
    args.append(tag)

with open('mk_Ji.out', 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################

prog = 'varplot'

args = [ prog ]
args.append('-w')
args.append(3)
args.append('-s')
args.append(0.3)
args.append('-g')
args.append(0.2)
args.append('mk_Ji.out')

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)
