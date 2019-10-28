#!/usr/bin/env python

import errno
import sys
import os
import subprocess
from Population import Population
import gd_composite
from dpmix_plot import make_dpmix_plot
from LocationFile import LocationFile

################################################################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno <> errno.EEXIST:
            raise

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

if len(sys.argv) < 16:
    print "usage"
    sys.exit(1)

input, input_type, data_source, switch_penalty, ap1_input, ap2_input, p_input, output, output2, output2_dir, dbkey, ref_column, galaxy_data_index_dir, heterochromatin_loc_file = sys.argv[1:15]
individual_metadata = sys.argv[15:]

chrom = 'all'
add_logs = '0'

loc_path = os.path.join(galaxy_data_index_dir, heterochromatin_loc_file)
location_file = LocationFile(loc_path)
heterochrom_path = location_file.get_values_if_exists(dbkey)
if heterochrom_path is None:
    heterochrom_path = '/dev/null'

population_list = []

p_total = Population()
p_total.from_tag_list(individual_metadata)

ap1 = Population(name='Ancestral population 1')
ap1.from_population_file(ap1_input)
population_list.append(ap1)
if not p_total.is_superset(ap1):
    print >> sys.stderr, 'There is an individual in ancestral population 1 that is not in the SNP table'
    sys.exit(1)

ap2 = Population(name='Ancestral population 2')
ap2.from_population_file(ap2_input)
population_list.append(ap2)
if not p_total.is_superset(ap2):
    print >> sys.stderr, 'There is an individual in ancestral population 2 that is not in the SNP table'
    sys.exit(1)

p = Population(name='Potentially admixed')
p.from_population_file(p_input)
population_list.append(p)
if not p_total.is_superset(p):
    print >> sys.stderr, 'There is an individual in the population that is not in the SNP table'
    sys.exit(1)

mkdir_p(output2_dir)

################################################################################
# Create tabular file
################################################################################

misc_file = os.path.join(output2_dir, 'misc.txt')

prog = 'dpmix'
args = [ prog ]
args.append(input)
args.append(ref_column)
args.append(chrom)
args.append(data_source)
args.append(add_logs)
args.append(switch_penalty)
args.append(heterochrom_path)
args.append(misc_file)

columns = ap1.column_list()
for column in columns:
    if input_type == 'gd_genotype':
        args.append('{0}:1:{1}'.format(int(column) - 2, ap1.individual_with_column(column).name))
    else:
        args.append('{0}:1:{1}'.format(column, ap1.individual_with_column(column).name))

columns = ap2.column_list()
for column in columns:
    if input_type == 'gd_genotype':
        args.append('{0}:2:{1}'.format(int(column) - 2, ap2.individual_with_column(column).name))
    else:
        args.append('{0}:2:{1}'.format(column, ap2.individual_with_column(column).name))

columns = p.column_list()
for column in columns:
    if input_type == 'gd_genotype':
        args.append('{0}:0:{1}'.format(int(column) - 2, p.individual_with_column(column).name))
    else:
        args.append('{0}:0:{1}'.format(column, p.individual_with_column(column).name))

run_program(None, args, stdout_file=output, space_to_tab=True)

################################################################################
# Create pdf file
################################################################################

pdf_file = os.path.join(output2_dir, 'dpmix.pdf')
make_dpmix_plot(dbkey, output, pdf_file, galaxy_data_index_dir)

################################################################################
# Create html
################################################################################

info_page = gd_composite.InfoPage()
info_page.set_title('dpmix Galaxy Composite Dataset')

display_file = gd_composite.DisplayFile()
display_value = gd_composite.DisplayValue()

out_pdf = gd_composite.Parameter(name='dpmix.pdf', value='dpmix.pdf', display_type=display_file)
out_misc = gd_composite.Parameter(name='misc.txt', value='misc.txt', display_type=display_file)

info_page.add_output_parameter(out_pdf)
info_page.add_output_parameter(out_misc)

if data_source == '0':
    data_source_value = 'sequence coverage'
elif data_source == '1':
    data_source_value = 'estimated genotype'

in_data_source = gd_composite.Parameter(description='Data source', value=data_source_value, display_type=display_value)
in_switch_penalty = gd_composite.Parameter(description='Switch penalty', value=switch_penalty, display_type=display_value)

info_page.add_input_parameter(in_data_source)
info_page.add_input_parameter(in_switch_penalty)

misc_populations =  gd_composite.Parameter(name='Populations', value=population_list, display_type=gd_composite.DisplayPopulationList())

info_page.add_misc(misc_populations)

with open(output2, 'w') as ofh:
    print >> ofh, info_page.render()

sys.exit(0)


