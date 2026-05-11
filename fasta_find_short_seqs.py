#!/usr/bin/env python3

"""
Given a reference FASTA and a queries FASTA, mark queries that are too short.

Usage:
    fasta_find_short_seqs.py <REFERENCE> <QUERIES> [--metric=STR] [--threshold=FLOAT] [--verbose] [--help]

Arguments:
    <REFERENCE>                  A FASTA file containing the reference sequences.
    <QUERIES>                    A FASTA file containing the query sequences.
    --metric STR                 The metric to use for comparison. One of 'mean', 'median', 'max' [default: max].
    --threshold FLOAT            The threshold for marking short sequences [default: 0.7].
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
    metric = args['--metric']
    threshold = float(args['--threshold'])
    out_prefix = queries_file.split('.')[0].split('/')[-1]

    if metric not in ['mean', 'median', 'max']:
        print("Error: metric must be one of 'mean', 'median', 'max'.")
        exit(1)

    #### MAIN ####
    # Get the length of all reference sequences
    print("Reading reference file %s..." % reference_file)
    len_list = []
    for seq_record in SeqIO.parse(reference_file, "fasta"):
        len_list.append(len(seq_record))
    print("Done reading reference file %s." % reference_file)

    # Compute the metric on the reference lengths
    if metric == 'mean':
        ref_metric = np.mean(len_list)
    elif metric == 'median':
        ref_metric = np.median(len_list)
    elif metric == 'max':
        ref_metric = np.max(len_list)
    
    # Iterate over the query sequences and mark those that are too short
    print("Reading query file %s..." % queries_file)
    out_string = ''
    with open(out_prefix + '_short_seqs.txt', 'w') as out_short:
        for seq_record in SeqIO.parse(queries_file, "fasta"):
            seq_len = len(seq_record)
            seqid = seq_record.id
            if seq_len < threshold * ref_metric:
                out_string += '%s\n' % seqid
        out_short.write(out_string)
    print("Done reading query file %s." % queries_file)