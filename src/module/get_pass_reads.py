#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse

# Import third-party modules
# import pybedtools
import pysam
import re

class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""

def get_pass_read_run(ibam, output_file, minimum_polya, polya_direction):
    # Input BAM file
    bam_file = pysam.AlignmentFile(ibam, 'rb')
    count = 0
    for line in bam_file:
        chrom = bam_file.get_reference_name(line.reference_id)
        ed = int(line.reference_end)
        # Convert pysam.calignedsegment.AlignedSegment object -> string oject
        line = str(line).rstrip()
        data = line.split("\t")

        # Check polyA length
        polyA_length = int(data[0].split('_')[-1].split('.')[0].split('A')[0])
        if polyA_length < minimum_polya:
            continue

        st = int(data[3])
        strand = data[1]
        seq_length = int(data[8])
        polya_st = ''
        polya_ed = ''
        tail_st = ''
        tail_ed = ''
        if polya_direction == 'three':
            if strand == "0":    # '+'
                polya_st = ed
                polya_ed = ed + 1
                tail_st = polya_st
                tail_ed = polya_st + polyA_length
                strand = '+'
            elif strand == "16":    # '-'
                polya_st = st
                polya_ed = st + 1
                tail_st = polya_ed - polyA_length
                tail_ed = polya_ed
                strand = '-'
            else:
                continue
        elif polya_direction == 'five':
            if strand == "16":    # '-'
                polya_st = ed
                polya_ed = ed + 1
                tail_st = polya_st
                tail_ed = polya_st + polyA_length
                strand = '+'
            elif strand == "0":    # '+'
                polya_st = st
                polya_ed = st + 1
                tail_st = polya_ed - polyA_length
                tail_ed = polya_ed
                strand = '-'
            else:
                continue
        name = "{0}||{1}||{2}||{3}||{4}||{5}".format(chrom, polya_st, polya_ed, data[0], "0", strand)
        print(chrom, tail_st, tail_ed, name, "0", strand, sep="\t", end="\n", file=output_file)

        # Splicing: 18M320N22M => st+18+320+22 = ed_exon2, st+18+320 = st_exon2
        # D +, I - from reads
        # WARNING: Cannot Detect Soft-clipping!!
        # if re.match('.+D.+', data[5]):
        #     print(line)
        #     print(ed)
        #     count += 1
        #     if count >= 30:
        #         break
        # ciger_string_list = re.findall(r'(\d+\D+)', data[5])
        # if re.match('.+N.+', data[5]):
        #     ciger_string_list = re.findall(r'(\d+\D+)', data[5])
        #     for ciger in ciger_string_list:
        #         if re.match('\d+', ciger):
        #             continue
        #     print(line, end="\n", file=output_file)



if __name__ == '__main__':
    main()
