#!/usr/bin/env python

import sys
import subprocess
from Population import Population

def run_program(prog, args, stdout_file=None, space_to_tab=False):
    #print "args: ", ' '.join(args)
    p = subprocess.Popen(args, bufsize=-1, executable=prog, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) = p.communicate()
    rc = p.returncode

    if stdout_file is not None:
        with open(stdout_file, 'w') as ofh:
            lines = stdoutdata.split('\n')
            for line in lines:
                line = line.strip()
                if line:
                    if space_to_tab:
                        line = line.replace(' ', '\t')
                    print >> ofh, line

    if rc != 0:
        print >> sys.stderr, "FAILED: rc={0}: {1}".format(rc, ' '.join(args))
        print >> sys.stderr, stderrdata
        sys.exit(1)

################################################################################

if len(sys.argv) < 10:
    print >> sys.stderr, "Usage"
    sys.exit(1)

snp_input, indel_input, coverage_input, annotation_input, indiv_input, ref_name, min_coverage, output = sys.argv[1:9]
individual_metadata = sys.argv[9:]

p_total = Population()
p_total.from_tag_list(individual_metadata)

p1 = Population()
p1.from_population_file(indiv_input)
if not p_total.is_superset(p1):
    print >> sys.stderr, 'There is an individual in the population individuals that is not in the SNP table'
    sys.exit(1)

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

run_program(prog, args, stdout_file='mk_Ji.out')

################################################################################

prog = 'varplot'

args = [ prog ]
args.append('-w')
args.append('3')
args.append('-s')
args.append('0.3')
args.append('-g')
args.append('0.2')
args.append('mk_Ji.out')

run_program(prog, args, stdout_file=output)
sys.exit(0)
