#!/usr/bin/env python

from __future__ import print_function
import sys
import string

# Import third-party modules
# import pybedtools
import pysam
import re

class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""

def internal_polyA_checker(seq, chrom, polya_st, polya_ed, name, strand):
    # Check internal polyA
    flg_internal_polya = 1
    output_dict = {}
    st_PASS_peak = ''
    ed_PASS_peak = ''
    polyA_count = 0
    for nucleotide in seq:
        if nucleotide == "A":
            polyA_count += 1
            continue
        else:
            if strand == '+':
                st_PASS_peak = int(polya_st) + polyA_count
                ed_PASS_peak = int(polya_ed) + polyA_count
            elif strand == '-':
                st_PASS_peak = int(polya_st) - polyA_count
                ed_PASS_peak = int(polya_ed) - polyA_count
            name = "{0}_NON-{1}A".format(name, str(polyA_count))

            # whether to pass >=2 nongenomic As
            remaining_polyA_length = len(seq[polyA_count:])
            if remaining_polyA_length >= 2:
                flg_internal_polya = 0
                break
            else:
                flg_internal_polya = 1
                break
    if flg_internal_polya == 1:
        output_dict["NG"] = 1
    elif flg_internal_polya == 0:
        output_dict["OK"] = [chrom, st_PASS_peak, ed_PASS_peak, name, strand, polyA_count]
    return output_dict

def get_pass_read_run(ibam, gfasta, output_file, log_file, internal_polyA_file, minimum_polya, polya_direction):
    # Input BAM file
    bam_file = pysam.AlignmentFile(ibam, 'rb')
    genome_fasta_file = pysam.FastaFile(gfasta)

    # Reverse sequence function
    rev_table = string.maketrans('ACGTacgt', 'TGCAtgca')
    def revcomp(seq, rev_table):
        return seq.translate(rev_table)[::-1]

    polyA_dict = {}
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

        # Determine internal polyA candidates
        if polya_direction == 'three':
            if strand == "0":    # '+'
                polya_st = ed
                polya_ed = ed + 1
                tail_st = polya_st
                tail_ed = polya_st + polyA_length
                strand = '+'
            elif strand == "16":    # '-'
                polya_st = st - 1
                polya_ed = st
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
                polya_st = st - 1
                polya_ed = st
                tail_st = polya_ed - polyA_length
                tail_ed = polya_ed
                strand = '-'
            else:
                continue
        name = data[0]

        # Get sequence
        seq = genome_fasta_file.fetch(chrom, tail_st, tail_ed)
        if strand == '-':
            seq = revcomp(seq, rev_table)

        # Check internal polyA
        result_dict = internal_polyA_checker(seq, chrom, polya_st, polya_ed, name, strand)
        st_PASS_peak = ''
        ed_PASS_peak = ''
        if "OK" in result_dict:
            info = result_dict["OK"]
            chrom = info[0]
            st_PASS_peak = info[1]
            ed_PASS_peak = info[2]
            name = info[3]
            polyA_count = info[5]
            print(chrom, st_PASS_peak, ed_PASS_peak, name, "0", strand, sep="\t", end="\n", file=output_file)

            # internal polyA count
            if not polyA_count in polyA_dict:
                polyA_dict[polyA_count] = 1
            else:
                polyA_dict[polyA_count] += 1
            count += 1    # PASS read count
        elif "NG" in result_dict:
            print(chrom, st, ed, name, seq, strand, sep="\t", end="\n", file=internal_polyA_file)

    for number in sorted(polyA_dict):
        polyA_count = polyA_dict[number]
        print(str(number) + "A", polyA_count, sep="\t", end="\n", file=log_file)

    print("PASS Reads: " + str(count), file=log_file)

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
