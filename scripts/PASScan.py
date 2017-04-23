#!/home/akimitsu/miniconda2/bin/python
# -*- coding: utf-8 -*-
#
# PASScan: Determine polyA sites from 3'end sequence data.
# Copyright (c) 2017- Naoto Imamachi <naoto.imamachi@gmail.com>
#
# Citation:
#

from __future__ import print_function
import sys
import os
from datetime import datetime
import time

# Import third-party modules
import click
from modules.trim_polyA import trim_polya_run
from modules.get_pass_reads import get_pass_read_run

# Print out comment with now time
def now_time(comment):
    now = datetime.now()
    nowtime = "{0:%Y-%m-%d %H:%M:%S}".format(now)
    print ('[' + nowtime + ']' + ' ' + comment)

def get_absolute_path(filename):
    '''Get absolute path of a file'''
    return os.path.abspath(os.path.expanduser(filename))

def time_checker(end_time):
    end_h = int(end_time/3600)
    end_time -= 3600 * end_h
    if end_h < 10:
        end_h = "0" + str(end_h)
    else:
        end_h = str(end_h)

    end_m = int(end_time/60)
    end_time -= 60 * end_m
    if end_m < 10:
        end_m = "0" + str(end_m)
    else:
        end_m = str(end_m)

    end_s = int(end_time)
    if end_s < 10:
        end_s = "0" + str(end_s)
    else:
        end_s = str(end_s)
    return [end_h, end_m, end_s]

# Callback functions for value validation
def validate_ifastq(ctx, param, value):
    if not value:
        click.echo('Error: FASTQ file is not chosen.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    # Check isfile
    if not os.path.isfile(value):
        click.echo('Error: FASTQ file does not exist.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    return(value)

def validate_ibam(ctx, param, value):
    if not value:
        click.echo('Error: BAM file is not chosen.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    # Check isfile
    if not os.path.isfile(value):
        click.echo('Error: BAM file does not exist.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    return(value)

def validate_gfasta(ctx, param, value):
    if not value:
        click.echo('Error: Reference FASTA file is not chosen.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    # Check isfile
    if not os.path.isfile(value):
        click.echo('Error: Reference FASTA file does not exist.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    return(value)

def validate_output(ctx, param, value):
    if not value:
        click.echo('Error: Do not choose output PREFIX.\n')
        click.echo(ctx.get_help())
        ctx.exit()
    # Check isDirectory
    if not os.path.exists(os.path.dirname(get_absolute_path(value))):
        click.echo('WARNING: Directory {0} does not exist.'.format(value))
        click.echo('Creating the directory...\n')
        os.makedirs(os.path.dirname(get_absolute_path(value)))
    return(value)

def validate_log(ctx, param, value):
    if not value:
        pass
    # Check isDirectory
    elif not os.path.exists(os.path.dirname(get_absolute_path(value))):
        click.echo('WARNING: The directory to save log file does not exist.')
        click.echo('Creating the directory...\n')
        os.makedirs(os.path.dirname(get_absolute_path(value)))
    return(value)

def validate_trimmed_nucleotides(ctx, param, value):
    for x in value:
        if not x in ["A", "T", "G", "C", "N", ',']:
            click.echo('Error: Nucleotide_trimming_list "{0}" is invalid.\n'.format(value))
            click.echo(ctx.get_help())
            ctx.exit()
    value_list = list(set(value.rstrip(',\n').split(',')))
    return(value_list)

# Parameters
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def cmd():
    '''
    PASScan determine polyA sites from 3'end sequence data.

    Copyright (C) 2017- Naoto Imamachi <naoto.imamachi@gmail.com>
    '''
    pass

@cmd.command()
@click.option('-i', '--ifastq', 'ifastq', callback=validate_ifastq, help='Input file in FASTQ file. [required]')
@click.option('-o', '--output-prefix', 'output', callback=validate_output, help='Prefix of output file. [required]')
@click.option('-l', '--log-prefix', 'log', callback=validate_log, help='Prefix of log file. Default: None')
@click.option('-nl', '--nucleotide-trimming-list', 'nucleotide_trimming_list', callback=validate_trimmed_nucleotides, default='A,N', help='Comma-delimited trimmed nucleotide list such as "A,N" or "A,C,N". Default: A,N')
@click.option('-nd', '--nucleotide-trimming-direction', 'nucleotide_trimming_direction', type=click.Choice(['three', 'five']), default='three', help='Trimmed nucleotide direction. Default: three')
@click.option('-m', '--minimum-length', 'minimum_length', type=int, default=18, help='Discard trimmed reads that are shorter than LENGTH. Default: 18')
@click.option('-bl', '--barcode-seq-length', 'barcode_seq_length', type=int, default=None, help='Length of barcode sequence. Trim the Barcode sequence from reads. Default: None')
@click.option('-bd', '--barcode-direction', 'barcode_direction', type=click.Choice(['three', 'five']), default='five', help='Trimmed nucleotide direction. Default: five')
def trim_polya(ifastq, output, log, nucleotide_trimming_list, nucleotide_trimming_direction, minimum_length, barcode_seq_length, barcode_direction):
    start_time = time.time()
    now_time("Beginning PASScan run (v0.1.0)")
    print("-"*50)
    now_time('Start trimming polyA sequence from reads...')
    click.echo("  Input FASTQ file: {0}".format(ifastq))
    click.echo("  Output prefix: {0}".format(output))
    if not log:
        log = "{0}.log".format(ifastq)
    click.echo("  Log prefix: {0}".format(log))
    click.echo("  Nucleotide trimming list: {0}".format(",".join(nucleotide_trimming_list)))
    click.echo("  Nucleotide trimming direction: {0}".format(nucleotide_trimming_direction))
    click.echo("  Minimum trimmed read length: {0}".format(minimum_length))
    click.echo("  Barcode sequence length: {0}".format(barcode_seq_length))
    click.echo("  Barcode direction: {0}".format(barcode_direction))

    # Trimming polyA/T and Barcode sequence from FASTQ reads
    fastq_file = open(ifastq, 'r')
    log_file = open(log, 'w')
    output_file = open(output, 'w')
    trim_polya_run(fastq_file, nucleotide_trimming_list, nucleotide_trimming_direction,
                   minimum_length, output_file, log_file,
                   barcode_seq_length, barcode_direction)

    # logging - Elapsed time
    end_time = time.time() - start_time
    finish_time_list = time_checker(end_time)
    now_time("PASScan - trim_polya was successfully finished: {0}:{1}:{2} elapsed".format(finish_time_list[0], finish_time_list[1], finish_time_list[2]))

@cmd.command()
@click.option('-i', '--ibam', 'ibam', callback=validate_ibam, help='Input file in BAM file. [required]')
@click.option('-g', '--reference-genome-fasta', 'gfasta', callback=validate_gfasta, help='Input file in Reference FASTA file. [required]')
@click.option('-o', '--output-prefix', 'output', callback=validate_output, help='Prefix of output file. [required]')
@click.option('-l', '--log-prefix', 'log', callback=validate_log, help='Prefix of log file. Default: None')
@click.option('-m', '--minimum-polya-length', 'minimum_polya', type=int, default=6, help='Discard polyA site-supporting reads that are shorter than polyA LENGTH [INT]. Default: 6')
@click.option('-d', '--polya-direction', 'polya_direction', type=click.Choice(['three', 'five']), default='three', help='PolyA direction on reads. Default: three')
def get_pass_read(ibam, gfasta, output, log, minimum_polya, polya_direction):
    start_time = time.time()
    now_time("Beginning PASScan run (v0.1.0)")
    print("-"*50)
    now_time('Start extracting a poly(A) site-supporting (PASS) reads from BAM file...')
    click.echo("  Input BAM file: {0}".format(ibam))
    click.echo("  Reference FASTA file: {0}".format(gfasta))
    click.echo("  Output prefix: {0}".format(output))
    if not log:
        log = "{0}.log".format(os.path.splitext(output)[0])
    click.echo("  Log prefix: {0}".format(os.path.splitext(output)[0]))
    click.echo("  Minimum polyA length: {0}".format(minimum_polya))
    click.echo("  PolyA Direction on reads: {0}".format(polya_direction))

    # Extracting PASS reads from BAM file
    output_file = open(output, 'w')
    internal_polyA_file = open("{0}_internal_polyA.bed".format(os.path.splitext(output)[0]), 'w')
    log_file = open(log, 'w')
    get_pass_read_run(ibam, gfasta, output_file, log_file, internal_polyA_file, minimum_polya, polya_direction)

    # logging - Elapsed time
    end_time = time.time() - start_time
    finish_time_list = time_checker(end_time)
    now_time("PASScan - get_pass_read was successfully finished: {0}:{1}:{2} elapsed".format(finish_time_list[0], finish_time_list[1], finish_time_list[2]))


def main():
    cmd()

if __name__ == '__main__':
    main()
