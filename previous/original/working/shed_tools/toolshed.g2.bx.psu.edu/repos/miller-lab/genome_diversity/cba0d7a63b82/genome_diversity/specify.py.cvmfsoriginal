#!/usr/bin/env python

import sys
import base64

def parse_args(args):
    if len(args) < 3:
        usage()

    input_file, output_file = args[1:3]

    individuals = []
    checkboxes = []
    strings = []

    for arg in args[3:]:
        if ':' in arg:
            arg_type, arg = arg.split(':', 1)
        else:
            print >> sys.stderr, "unknown argument:", arg
            usage()

        if arg_type == 'individual':
            individuals.append(arg)
        elif arg_type == 'checkbox':
            checkboxes.append(arg)
        elif arg_type == 'string':
            strings.append(arg)
        else:
            print >> sys.stderr, "unknown argument:", arg
            usage()

    return input_file, output_file, individuals, checkboxes, strings

def usage():
    print >> sys.stderr, "Usage: %s <input> <output> [<individual:col:name> ...] [<checkbox:col:name> ...] [<string:base64> ...]" % (sys.argv[0])
    sys.exit(1)

def parse_individuals(individuals):
    ind_col2name = {}
    ind_name2col = {}

    for individual in individuals:
        if ':' in individual:
            column, name = individual.split(':', 1)
        else:
            print >> sys.stderr, "invalid individual specification:", individual
            usage()

        try:
            column = int(column)
        except:
            print "individual column is not an integer:", individual
            usage()

        if column not in ind_col2name:
            ind_col2name[column] = name
        else:
            if ind_col2name[column] != name:
                print "duplicate individual column:", name, column, ind_col2name[column]
                usage()

        if name not in ind_name2col:
            ind_name2col[name] = [column]
        elif column not in ind_name2col[name]:
            ind_name2col[name].append(column)

    return ind_col2name, ind_name2col

def parse_checkboxes(checkboxes, ind_col2name):
    columns = []

    for checkbox in checkboxes:
        if ':' in checkbox:
            column, name = checkbox.split(':', 1)
        else:
            print >> sys.stderr, "invalid checkbox specification:", checkbox
            usage()

        try:
            column = int(column)
        except:
            print "checkbox column is not an integer:", checkbox
            usage()

        if column not in ind_col2name:
            print "individual not in SNP table:", name
            usage()

        if column not in columns:
            columns.append(column)

    return columns

def parse_strings(strings, ind_col2name, ind_name2col):
    columns = []

    for string in strings:
        try:
            decoded = base64.b64decode(string)
        except:
            print >> sys.stderr, "invalid base64 string:", string
            usage()

        names = find_names(decoded, ind_name2col.keys())
        for name in names:
            cols = ind_name2col[name]
            if len(cols) == 1:
                col = cols[0]
                if col not in columns:
                    columns.append(col)
            else:
                print >> sys.stderr, "name with multiple columns:", name
                usage()

    return columns

def find_names(string, names):
    rv = []
    for name in names:
        if name in string:
            if name not in rv:
                rv.append(name)
    return rv




input_file, output_file, individuals, checkboxes, strings = parse_args(sys.argv)
ind_col2name, ind_name2col = parse_individuals(individuals)
cb_cols = parse_checkboxes(checkboxes, ind_col2name)
str_cols = parse_strings(strings, ind_col2name, ind_name2col)

out_cols = cb_cols
for col in str_cols:
    if col not in out_cols:
        out_cols.append(col)

with open(output_file, 'w') as fh:
    for col in sorted(out_cols):
        print >> fh, '\t'.join([str(x) for x in [col, ind_col2name[col], '']])

sys.exit(0)




