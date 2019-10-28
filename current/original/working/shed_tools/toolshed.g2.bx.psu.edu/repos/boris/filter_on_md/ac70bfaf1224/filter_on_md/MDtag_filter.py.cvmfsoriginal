#!/usr/bin/env python
# Boris Rebolledo Jaramillo
# berebolledo@gmail.com
# v1.0.2

"""Filter mapped reads on MD tag string.

This tool reads the MD tag of mapped reads (see SAM format specification).
The user defines the 5' and 3' windows n and m (in bp), respectively.
The mapped read is discarded if it contains any number of mismatches within n
bases of the read 5' end and within m bases of the read 3' end.
Option: save discarded reads in an additional SAM file.

usage: %prog in.sam n m out.sam saveDiscarded?[yes/no] [discarded.sam]

"""

import sys,re

if len(sys.argv) < 6:
    sys.exit('ERROR! Usage: %prog <in.sam> <n> <m> <out.sam> '
            '<saveDiscarded?[yes/no]> [<discarded.sam>]')
if sys.argv[5] == 'yes' and len(sys.argv) < 7:
    sys.exit('ERROR! Usage: %prog <in.sam> <n> <m> <out.sam> '
            '<saveDiscarded?[yes/no]> [<discarded.sam>]')      
      
sam_file=list(open(sys.argv[1]))
sam_mdfiltered=open(sys.argv[4],'w+')
if sys.argv[5] == 'yes':
    sam_discarded=open(sys.argv[6],'w+')


# The MD tag is generated out of the alignment operations, regardless the strand
# the read was mapped to. It describes the alignment from the leftmost aligned
# base, which might be either the 5' or 3' end of the original read.
# Luckily, the read orientation in the alignment can be extracted from the
# alignment flag.
# The code used here to identify the mapped reads orientation was taken from
# the 'explain_sam_flags.py' script of the Picard-Tools.
# It defines the meaning of the sam flags in plain English, so they can be used
# to identify reads mapped to the reverse strand: flag = 16
# http://picard-tools.sourcearchive.com/documentation/1.25-1/

lstFlags = [
    ("read paired", 0x1),
    ("read mapped in proper pair", 0x2),
    ("read unmapped", 0x4),
    ("mate unmapped", 0x8),
    ("read reverse strand", 0x10),
    ("mate reverse strand", 0x20),
    ("first in pair", 0x40),
    ("second in pair", 0x80),
    ("not primary alignment", 0x100),
    ("read fails platform/vendor quality checks", 0x200),
    ("read is PCR or optical duplicate", 0x400)
    ]

for line in sam_file:
    if line.split("\t")[0].startswith('@'):
        if sys.argv[5] == 'yes':
            sam_discarded.write(line)
            sam_mdfiltered.write(line)
        else:
            sam_mdfiltered.write(line)
    else:
        if len(line.split("\t")) > 11:        # Any SAM optional fields? 
            for field in line.split("\t")[11:]:
                if re.match('MD:Z:[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*', field):        # is MD tag field present?
                    flag = int(line.split("\t")[1]) 
                    for strFlagName, iMask in lstFlags:
                        if flag & iMask:
                            if iMask == 16: 
                                strand = 'reverse'
                                break
                            else:
                                strand = 'forward'
                        else:        # flag = 0 is not part of the definitions,
                                     # but it corresponds to forward strand mapping
                            strand = 'forward'
                    # The position of the optional fields in the SAM format is variable. Finds the location of the MD tag:
                    mdtag_idx = [i for i, item in enumerate(line.split("\t")[11:]) if re.match('MD:Z:[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*', item)][0]
                    mdtag=line.split("\t")[mdtag_idx+11].strip().split(":")[2]
                    md_list = re.split(r'(\D)',mdtag)
                    if strand == 'forward':
                        if int(md_list[0]) >= int(sys.argv[2]) and int(md_list[-1]) >= int(sys.argv[3]):
                            sam_mdfiltered.write(line)
                        else:
                            if sys.argv[5] == 'yes':
                                sam_discarded.write(line)
                    elif strand == 'reverse':
                        if int(md_list[0]) >= int(sys.argv[3]) and int(md_list[-1]) >= int(sys.argv[2]):
                            sam_mdfiltered.write(line)
                        else:
                            if sys.argv[5] == 'yes':
                                sam_discarded.write(line)
                    break

sam_mdfiltered.close()

if sys.argv[5] == 'yes':
    sam_discarded.close()