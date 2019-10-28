#!/usr/bin/env python

import gd_util
import sys
from Population import Population

################################################################################

if len(sys.argv) != 7:
    gd_util.die('Usage')

snp_input, coverage_input, indiv_input, min_coverage, output, ind_arg = sys.argv[1:]

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

p1 = Population()
p1.from_population_file(indiv_input)
if not p_total.is_superset(p1):
    gd_util.die('There is an individual in the population individuals that is not in the SNP table')

################################################################################

prog = 'mt_pi'

args = [ prog ]
args.append(snp_input)
args.append(coverage_input)
args.append(min_coverage)

for tag in p1.tag_list():
    args.append(tag)

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

sys.exit(0)

