#!/usr/bin/python

import sys
import optparse
import csv
import re

def main():
	parser = optparse.OptionParser()
	parser.add_option('-s', '--snpeff_file', dest = 'snpeff_file', action = 'store', type = 'string', default = None, help = "Path to the snpEff file")
	parser.add_option('-c', '--candidate_list', dest = 'candidate_list', action = 'store', type = 'string', default = None, help = "Two column tabular list of candidate Gene ID, Type")
	parser.add_option('-o', '--output', dest = 'output', action = 'store', type = 'string', default = None, help = "Output file name")	
	(options, args) = parser.parse_args()
	
	snpeff_file = options.snpeff_file
	candidate_list = options.candidate_list

	candidates = parse_candidate_list(candidate_list = candidate_list)
	mark_snpeff_file(snpeff_file = snpeff_file, output = options.output, candidates = candidates)

def skip_and_write_headers(writer = None, reader = None, i_file = None):
	# count headers
	comment = 0
	while reader.next()[0].startswith('#'):
		comment = comment + 1
	
	# skip and write headers
	i_file.seek(0)
	for i in range(0, comment):
		row = reader.next()
		writer.writerow(row)

def parse_candidate_list(candidate_list = ""):
	i_file = open(candidate_list, 'rU')
	reader  = csv.reader(i_file, delimiter = '\t',)

	candidates = {}
	for row in reader:
		gene_id = row[0]
		gene_type = row[1]
		candidates[gene_id] = gene_type

	i_file.close()
	
	return candidates

def mark_snpeff_file(snpeff_file = "", output = "", candidates = None):
	i_file = open(snpeff_file, 'rU')
	reader = csv.reader(i_file, delimiter = '\t')

	o_file = open(output, 'wb')
	writer = csv.writer(o_file, delimiter = '\t')

	skip_and_write_headers(writer = writer, reader = reader, i_file = i_file)

	for row in reader:
		gene_id = row[9]
		if gene_id in candidates:
			gene_type = candidates[gene_id]
			row.append(gene_type)
		else:
			row.append('')

		writer.writerow(row)

	o_file.close()
	i_file.close()

if __name__ == "__main__":
    main()
