#!/usr/bin/env python

import sys
import subprocess
from Population import Population

################################################################################

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

if len(sys.argv) < 8:
    print >> sys.stderr, "Usage"
    sys.exit(1)

gd_saps_file, gd_snps_file, covered_intervals_file, gd_indivs_file, output_file = sys.argv[1:6]
individual_metadata = sys.argv[6:]

p_total = Population()
p_total.from_tag_list(individual_metadata)

p1 = Population()
p1.from_population_file(gd_indivs_file)
if not p_total.is_superset(p1):
    print >> sys.stderr, 'There is an individual in the population individuals that is not in the SNP table'
    sys.exit(1)

################################################################################

prog = 'get_pi'

args = [ prog ]
args.append(gd_saps_file)
args.append(gd_snps_file)
args.append(covered_intervals_file)

columns = p1.column_list()
for column in columns:
    args.append('{0}'.format(column))

run_program(None, args, stdout_file=output_file)

sys.exit(0)

