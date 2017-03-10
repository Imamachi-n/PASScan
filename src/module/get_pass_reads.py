#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse

# Import third-party modules
import pysam

class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""

def main():
    parser = argparse.ArgumentParser(description='Extract a poly(A) site-supporting (PASS) reads from BAM file.')
    parser.add_argument('-i', '--ibam', dest='ibam', help='Input file in BAM file. [required]')
    parser.add_argument('-o', '--output-prefix', dest='output', help='Prefix of output file. [required]')
    parser.add_argument('-a', '--minimum-polya-length', dest='minimum_polya', type=int, help='Discard polyA site-supporting reads that are shorter than polyA LENGTH [INT]. Default: 6')
    args = parser.parse_args()

    # output
    output_file = open(args.output, 'w')
    minimum_polya_length = 6
    if args.minimum_polya:
        minimum_polya_length = args.minimum_polya

    # Input BAM file
    bam_file = pysam.AlignmentFile(args.ibam, 'rb')
    count = 0
    for line in bam_file:
        # Convert pysam.calignedsegment.AlignedSegment object -> string oject
        line = str(line).rstrip()
        print(line)
        data = line.split("\t")

        # Check polyA length
        polyA_length = int(data[0].split('_')[-1].split('.')[0].split('A')[0])
        print(polyA_length)
        # if polyA_length >= minimum_polya_length:
        count += 1
        if count >= 10:
            break

if __name__ == '__main__':
    main()
