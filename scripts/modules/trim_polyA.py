#!/usr/bin/env python

from __future__ import print_function
import sys

class FormatError(Exception):
	"""
	Raised when an input file (FASTA or FASTQ) is malformatted.
	"""

def polyN_counter(seq, nucleotide_list, direction):
    '''
    seq: Each fastq sequence read
    nucleotide_list: Trimming candidate nucleotides such as [A], [A, C], [A, N].
    direction: Choose 'five' or 'three'.
               Determine the trimming direction.
    '''
    trimmed_seq = seq
    polyN_count = 0

    # Trimming from 3'end
    if direction == 'five':
        for nuc in seq:
            if nuc in nucleotide_list:
                # Counting pAs
                polyN_count += 1
            else:
                break
        # Get trimmed sequence
        trimmed_seq = seq[polyN_count:]

    # Trimming from 5'end
    elif direction == 'three':
        for nuc in seq[::-1]:
            if nuc in nucleotide_list:
                # Counting pAs
                polyN_count += 1
            else:
                break
        trim_position = len(trimmed_seq) - polyN_count
        # Get trimmed sequence
        trimmed_seq = seq[:trim_position]
    return [trimmed_seq, polyN_count]

def barcode_trimming(seq, barcode_length, barcode_direction):
    barcode_seq = ''
    if barcode_direction == "five":
        barcode_seq = seq[:barcode_length]
        seq = seq[barcode_length:]
    elif barcode_direction == "three":
        trim_position = len(seq) - barcode_length
        barcode_seq = seq[trim_position:]
        seq = seq[:trim_position]
    return [seq, barcode_seq]

def trim_polya_run(fastq_file, nucleotide_list, trimming_direction, minimum_length, fastq_output, log_output, barcode_seq_length, barcode_direction):
    '''
    fastq_file: FASTQ input file
    trimming_direction: Choose 'five' or 'three'.
                        Determine the trimming direction.
    barcode_seq_length: Barcode sequence length
    barcode_direction: Barcode sequence direction
    '''
    i = 3
    line_dict = {}
    polyA_length_dict = {}
    for (i, line) in enumerate(fastq_file):
        line = line.rstrip()

        # Read fastq 1st line
        if i % 4 == 0:
            if not line.startswith('@'):
                raise FormatError("Line {0} in FASTQ file is expected to start with '@', but found {1!r}".format(i+1, line[:10]))
            line_dict["name"] = line

        # Read fastq 2nd line
        elif i % 4 == 1:
            raw_sequence = line

            # Trimming barcode sequence
            barcode_seq = 'N'
            if barcode_seq_length != None and barcode_direction != None:
                tmp = barcode_trimming(raw_sequence, barcode_seq_length, barcode_direction)
                raw_sequence = tmp[0]
                barcode_seq = tmp[1]

            # Trimming polyA/T sequence
            tmp_info = polyN_counter(raw_sequence, nucleotide_list, trimming_direction)
            trimmed_seq = tmp_info[0]
            polyA_length = tmp_info[1]

            # Prepare random barcode & polyA trimmed reads
            line_dict["polyA_length"] = polyA_length
            line_dict["trimmed_seq"] = trimmed_seq
            line_dict["barcode_seq"] = barcode_seq

        # Read fastq 3rd line
        elif i % 4 == 2:
            if not line.startswith('+'):
                raise FormatError("Line {0} in FASTQ file is expected to start with '+', but found {1!r}".format(i+1, line[:10]))
            if line[1:] != line_dict["name"][1:]:
                raise FormatError("Line {0}: Sequence descriptions in the FASTQ file do not match ({1!r} != {2!r}).".format(i+1, line_dict["name"][1:], line[1:]))

        # Read fastq 4th line
        elif i % 4 == 3:
            # Discard trimmed reads that are shorter than LENGTH
            trimmed_seq = line_dict["trimmed_seq"]
            if len(trimmed_seq) < 18:
                continue

            # Prepare renamed name
            polyA_length = line_dict["polyA_length"]
            barcode_seq = line_dict["barcode_seq"]
            name_core = line_dict["name"].split('.')[0]
            name_other = line_dict["name"].split('.')[1]
            renamed_core = "{0}_{1}A_{2}.".format(name_core, polyA_length, barcode_seq)
            renamed_name = renamed_core + name_other

            # Option line
            option_line = "+{0}".format(renamed_name[1:])

            quality_score = line
            if barcode_seq_length != None and barcode_direction != None:
                tmp = barcode_trimming(quality_score, barcode_seq_length, barcode_direction)
                quality_score = tmp[0]
                tmp = barcode_trimming(quality_score, polyA_length, trimming_direction)
                quality_score = tmp[0]

            # Logging polyA length
            if not polyA_length in polyA_length_dict:
                polyA_length_dict[polyA_length] = 1
            else:
                polyA_length_dict[polyA_length] += 1

            # Trimmed FASTQ output
            print(renamed_name, end="\n", file=fastq_output)
            print(trimmed_seq, end="\n", file=fastq_output)
            print(option_line, end="\n", file=fastq_output)
            print(quality_score, end="\n", file=fastq_output)

    if i % 4 != 3:
        raise FormatError("FASTQ file ended prematurely")

    # Logging polyA length info
    for polyA_number in sorted(polyA_length_dict):
        polyA_num_count = polyA_length_dict[polyA_number]
        print(str(polyA_number), str(polyA_num_count), sep="\t", end="\n", file=log_output)

if __name__ == '__main__':
    main()
