#!/usr/bin/env python3

"""
Given a reference FASTA and a queries FASTA, mark queries that are too short.

Usage:
    fasta_compare_seq_len.py <REFERENCE> <QUERIES> [--verbose] [--help]

Arguments:
    <REFERENCE>                  A FASTA file containing the reference sequences.
    <QUERIES>                    A FASTA file containing the query sequences.
    --verbose                    Print the progressions of the program to the terminal (Standard Error).
    --help                       Show this help message and exit.
"""

#### LIBS ####
from Bio import SeqIO           # to work with sequence files
from docopt import docopt       # to create the argument parser
import numpy as np              # to work with ndarrays and do maths


if __name__ == '__main__':
    #### PARSE ARGS ####
    args = docopt(__doc__)
    reference_file = args['<REFERENCE>']
    queries_file = args['<QUERIES>']
    out_prefix = queries_file.split('.')[0].split('/')[-1]

    #### MAIN ####
    # Get the length of all reference sequences
    print("Reading reference file %s..." % reference_file)
    len_list = []
    for seq_record in SeqIO.parse(reference_file, "fasta"):
        len_list.append(len(seq_record))
    print("Done reading reference file %s." % reference_file)

    # Compute the metric on the reference lengths
    ref_mean = np.mean(len_list)
    ref_median = np.median(len_list)
    ref_max = np.max(len_list)

    # Iterate over the query sequences and mark those that are too short
    print("Reading query file %s..." % queries_file)
    out_string = '1_seqid\t2_rel_len_mean\t3_rel_len_median\t4_rel_len_max\n'
    with open(out_prefix + '.rel_len.tsv', 'w') as wh:
        for seq_record in SeqIO.parse(queries_file, "fasta"):
            seq_len = len(seq_record)
            seqid = seq_record.id
            rel_len_mean = seq_len / ref_mean
            rel_len_median = seq_len / ref_median
            rel_len_max = seq_len / ref_max
            out_string += '%s\t%.2f\t%.2f\t%.2f\n' % (seqid, rel_len_mean, rel_len_median, rel_len_max)
        wh.write(out_string)
    print("Done reading query file %s." % queries_file)