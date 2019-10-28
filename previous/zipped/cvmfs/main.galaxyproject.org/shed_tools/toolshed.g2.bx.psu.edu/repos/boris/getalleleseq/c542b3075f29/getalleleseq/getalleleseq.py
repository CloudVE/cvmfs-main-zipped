#!/usr/bin/env python
# Boris Rebolledo-Jaramillo (boris-at-bx.psu.edu)
#
#usage: getalleleseq.py [-h] [-l INT] [-j FILE] [-d DIR] alleles
#
#Given a table with minor and major alleles per position, it generates the
#minor and major allele sequences in FASTA format
#
#positional arguments:
#  alleles               Table containing minor and major allele base per
#                        position. cols: [id, chr, pos, A, C, G, T, cvrg,
#                        plody, major, minor, freq_minor]
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -l INT, --seq-length INT
#                        Background sequence length. Bases in an artifical
#                        all-N-sequence of length INT will be replaced by
#                        either the major or minor allele base accordingly
#  -j FILE, --major-seq FILE
#                        File to write major allele sequences in FASTA multiple
#                        alignment format.
#  -d DIR, --minor-dir DIR
#                        Per sample minor allele sequences will be written to
#                        this directory
#
# The expected columns in the alleles table follow Nicholas Stoler's
# Variant Annotator tool format.  See Variant Annotator in Galaxy's tool shed
# http://testtoolshed.g2.bx.psu.edu/repos/nick/allele_counts_1 for more details
#
# Expected columns:
# 1.  sample_id
# 2.  chr
# 3.  position
# 4   counts for A's
# 5.  counts for C's
# 6.  counts for G's
# 7.  counts for T's
# (8.  counts for a's)
# (9.  counts for c's)
# (10. counts for g's)
# (11. counts for t's)
# 8.  (12.) Coverage
# 9.  (13.) Number of alleles passing a given criteria
# 10. (14.) Major allele
# 11. (15.) Minor allele
# 12. (16.) Minor allele frequency in position

import sys
import os
import argparse

def createseq(sample, allele, seq_size, table):
    """Generate major or minor allele sequence"""
    out_sequence = ['N' for i in range(seq_size)]
    sample_data  = [line for line in table if line[0] == sample]

    for entry in sample_data:
        position = int(entry[2])
        if len(entry)==12:
            number_of_alleles = int(entry[8])
            major_allele = entry[9].strip()
            minor_allele = entry[10].strip()
        else:
            number_of_alleles = int(entry[12])
            major_allele = entry[13].strip()
            minor_allele = entry[14].strip()

        if allele == 'major':
            out_sequence[position-1] = major_allele 
        elif allele == 'minor':
            if number_of_alleles >= 2: 
                out_sequence[position-1] = minor_allele 
            else:
                out_sequence[position-1] = major_allele 
    return out_sequence    

def printseq(sample,allele,seq,output):
    """Print out sequence"""
    #print >> output, '>{0}_{1}'.format(sample,allele)
    print >> output, '>{0}{1}'.format(sample,allele)
    for i in range(0,len(seq),70):
        print >> output, ''.join(seq[i:i+70])

def main():
    parser = argparse.ArgumentParser(description='Given a table with minor and major alleles per position, it generates the minor and major allele sequences in FASTA format', epilog='Boris Rebolledo-Jaramillo (boris-at-bx.psu.edu)')
    parser.add_argument('alleles', type=str, help='Table containing minor and major allele base per position. cols: [id, chr, pos, A, C, G, T, cvrg, plody, major, minor, freq_minor] ')
    parser.add_argument('-l','--seq-length', type=int, metavar='INT', help='Background sequence length. Bases in an artifical all-N-sequence of length INT will be replaced by either the major or minor allele base accordingly')
    parser.add_argument('-j','--major-seq', type=str, metavar='FILE', help='File to write major allele sequences in FASTA multiple alignment format.')
    parser.add_argument('-d', '--minor-dir', type=str, metavar='DIR', default='.', help="Per sample minor allele sequences will be written to this directory (Default: current directory)")
    parser.add_argument('-p', '--minor-prefix', type=str, metavar='STR', nargs='?', const='', default='', help=argparse.SUPPRESS) #Galaxy compatibility
    args = parser.parse_args()
    
    
    try:
        table = [line.strip().split('\t') for line in list(open(args.alleles)) if "#" not in line]
        samples = sorted(list(set([ line[0] for line in table ])))
    except:
        sys.exit('\nERROR: Could not open %s\n' % args.alleles)
    try:
        major_out = open(args.major_seq, 'w+')
    except:
        sys.exit('\nCould not create %s\n' % args.major_seq)

    # Single file for all major allele sequences in FASTA multiple alignment
    for sample in samples:
        sequence = createseq(sample,'major',args.seq_length,table)
        #printseq(sample,'major',sequence,major_out)
        printseq(sample,'',sequence,major_out)
    major_out.close()

    # Sample specific minor allele sequence in FASTA format
    try:
        os.makedirs(args.minor_dir)
    except:
        pass

    for sample in samples:
        if args.minor_prefix: # to fit Galaxy requirements
            name = sample.replace('_','')
            minor_name = "%s_%s_%s" % ('primary',args.minor_prefix,name+'-minor_visible_fasta')
        else: # for non-Galaxy 
            minor_name = sample+'-minor.fa'
        minor_out = open(os.path.join(args.minor_dir, minor_name), 'w+')
        sequence = createseq(sample,'minor',args.seq_length,table)
        #printseq(sample,'minor',sequence,minor_out)
        printseq(sample,'_minor',sequence,minor_out)
        minor_out.close()

if __name__ == "__main__": main()