#!/usr/bin/env python

import errno
import os
import shutil
import subprocess
import sys
from Population import Population
import gd_composite

################################################################################

def do_import(filename, files_path, min_reads, min_qual, min_spacing, tags, using_info, population_list):
    info_page = gd_composite.InfoPage()
    info_page.set_title('Prepare to look for population structure Galaxy Composite Dataset')

    display_file = gd_composite.DisplayFile()
    display_value = gd_composite.DisplayValue()

    out_ped = gd_composite.Parameter(name='admix.ped', value='admix.ped', display_type=display_file)
    out_map = gd_composite.Parameter(name='admix.map', value='admix.map', display_type=display_file)
    out_use = gd_composite.Parameter(description=using_info, display_type=display_value)

    info_page.add_output_parameter(out_ped)
    info_page.add_output_parameter(out_map)
    info_page.add_output_parameter(out_use)

    in_min_reads = gd_composite.Parameter(description='Minimum reads covering a SNP, per individual', value=min_reads, display_type=display_value)
    in_min_qual = gd_composite.Parameter(description='Minimum quality value, per individual', value=min_qual, display_type=display_value)
    in_min_spacing = gd_composite.Parameter(description='Minimum spacing between SNPs on the same scaffold', value=min_spacing, display_type=display_value)

    info_page.add_input_parameter(in_min_reads)
    info_page.add_input_parameter(in_min_qual)
    info_page.add_input_parameter(in_min_spacing)

    misc_populations = gd_composite.Parameter(name='Populations', value=population_list, display_type=gd_composite.DisplayPopulationList())
    info_page.add_misc(misc_populations)

    with open(filename, 'w') as ofh:
        print >> ofh, info_page.render()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno <> errno.EEXIST:
            raise

def die(message, exit=True):
    print >> sys.stderr, message
    if exit:
        sys.exit(1)

################################################################################

if len(sys.argv) < 9:
    die("Usage")

# parse command line
input_snp_filename, min_reads, min_qual, min_spacing, output_filename, output_files_path = sys.argv[1:7]
args = sys.argv[7:]

individual_metadata = []
population_files = []
population_names = []
all_individuals = False

for arg in args:
    if arg == 'all_individuals':
        all_individuals = True
    elif len(arg) > 11:
        tag = arg[:11]
        value = arg[11:]
        if tag == 'individual:':
            individual_metadata.append(value)
        elif tag == 'population:':
            filename, name = value.split(':', 1)
            population_files.append(filename)
            population_names.append(name)

p_total = Population()
p_total.from_tag_list(individual_metadata)

individual_population = {}

population_list = []

if all_individuals:
    p1 = p_total
    p1.name = 'All Individuals'
    population_list.append(p1)
else:
    p1 = Population()
    for idx in range(len(population_files)):
        population_file = population_files[idx]
        population_name = population_names[idx]
        this_pop = Population(population_name)
        this_pop.from_population_file(population_file)
        population_list.append(this_pop)
        p1.from_population_file(population_file)
        tags = p1.tag_list()
        for tag in tags:
            if tag not in individual_population:
                individual_population[tag] = population_name

if not p_total.is_superset(p1):
    print >> sys.stderr, 'There is an individual in the population that is not in the SNP table'
    sys.exit(1)

# run tool
prog = 'admix_prep'

args = []
args.append(prog)
args.append(input_snp_filename)
args.append(min_reads)
args.append(min_qual)
args.append(min_spacing)

tags = p1.tag_list()
for tag in tags:
    args.append(tag)

#print "args:", ' '.join(args)
p = subprocess.Popen(args, bufsize=-1, stdin=None, stdout=subprocess.PIPE, stderr=sys.stderr)
(stdoutdata, stderrdata) = p.communicate()
rc = p.returncode

if rc != 0:
    die('admix_prep failed: rc={0}'.format(rc))

using_info = stdoutdata.rstrip('\r\n')
mkdir_p(output_files_path)
output_ped_filename = os.path.join(output_files_path, 'admix.ped')
output_map_filename = os.path.join(output_files_path, 'admix.map')
shutil.copy2('admix.ped', output_ped_filename)
shutil.copy2('admix.map', output_map_filename)
do_import(output_filename, output_files_path, min_reads, min_qual, min_spacing, tags, using_info, population_list)

os.unlink('admix.ped')
os.unlink('admix.map')

sys.exit(0)

