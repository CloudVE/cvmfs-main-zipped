#!/usr/bin/env python

import gd_util
import sys
import os
from Population import Population
import gd_composite
from dpmix_plot import make_dpmix_plot
from LocationFile import LocationFile

################################################################################

if len(sys.argv) != 22:
    print "usage"
    sys.exit(1)

input, input_type, data_source, switch_penalty, ap1_input, ap1_name, ap2_input, ap2_name, ap3_input, ap3_name, p_input, output, output2, output2_dir, dbkey, ref_column, galaxy_data_index_dir, heterochromatin_loc_file, ind_arg, het_arg, add_logs = sys.argv[1:]

if ap3_input == '/dev/null':
    populations = 2
else:
    populations = 3

chrom = 'all'

if het_arg == 'use_installed':
    loc_path = os.path.join(galaxy_data_index_dir, heterochromatin_loc_file)
    location_file = LocationFile(loc_path)
    heterochrom_path = location_file.get_values_if_exists(dbkey)
    if heterochrom_path is None:
        heterochrom_path = '/dev/null'
elif het_arg == 'use_none':
    heterochrom_path = '/dev/null'
else:
    heterochrom_path = het_arg

population_list = []

p_total = Population()
p_total.from_wrapped_dict(ind_arg)

ap1 = Population(name='Ancestral population 1')
ap1.from_population_file(ap1_input)
population_list.append(ap1)
if not p_total.is_superset(ap1):
    gd_util.die('There is an individual in ancestral population 1 that is not in the SNP table')

ap2 = Population(name='Ancestral population 2')
ap2.from_population_file(ap2_input)
population_list.append(ap2)
if not p_total.is_superset(ap2):
    gd_util.die('There is an individual in ancestral population 2 that is not in the SNP table')

if populations == 3:
    ap3 = Population(name='Ancestral population 3')
    ap3.from_population_file(ap3_input)
    population_list.append(ap3)
    if not p_total.is_superset(ap3):
        gd_util.die('There is an individual in ancestral population 3 that is not in the SNP table')

p = Population(name='Potentially admixed')
p.from_population_file(p_input)
population_list.append(p)
if not p_total.is_superset(p):
    gd_util.die('There is an individual in the population that is not in the SNP table')

gd_util.mkdir_p(output2_dir)

################################################################################
# Create tabular file
################################################################################

misc_file = os.path.join(output2_dir, 'summary.txt')

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
    col = int(column)
    name = ap1.individual_with_column(column).name
    first_token = name.split()[0]
    if input_type == 'gd_genotype':
        col -= 2
    args.append('{0}:1:{1}'.format(col, first_token))

columns = ap2.column_list()
for column in columns:
    col = int(column)
    name = ap2.individual_with_column(column).name
    first_token = name.split()[0]
    if input_type == 'gd_genotype':
        col -= 2
    args.append('{0}:2:{1}'.format(col, first_token))

if populations == 3:
    columns = ap3.column_list()
    for column in columns:
        col = int(column)
        name = ap3.individual_with_column(column).name
        first_token = name.split()[0]
        if input_type == 'gd_genotype':
            col -= 2
        args.append('{0}:3:{1}'.format(col, first_token))

columns = p.column_list()
for column in columns:
    col = int(column)
    name = p.individual_with_column(column).name
    first_token = name.split()[0]
    if input_type == 'gd_genotype':
        col -= 2
    args.append('{0}:0:{1}'.format(col, first_token))

with open(output, 'w') as fh:
    gd_util.run_program(prog, args, stdout=fh)

################################################################################
# Create pdf file
################################################################################

if populations == 3:
    state2name = {
        0:'heterochromatin',
        1:ap1_name,
        2:ap2_name,
        3:ap3_name
    }
else:
    state2name = {
        0:'heterochromatin',
        1:ap1_name,
        2:ap2_name
    }

pdf_file = os.path.join(output2_dir, 'picture.pdf')
make_dpmix_plot(dbkey, output, pdf_file, galaxy_data_index_dir, state2name=state2name, populations=populations)

################################################################################
# Create html
################################################################################

info_page = gd_composite.InfoPage()
info_page.set_title('dpmix Galaxy Composite Dataset')

display_file = gd_composite.DisplayFile()
display_value = gd_composite.DisplayValue()

out_pdf = gd_composite.Parameter(name='picture.pdf', value='picture.pdf', display_type=display_file)
out_misc = gd_composite.Parameter(name='summary.txt', value='summary.txt', display_type=display_file)

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

