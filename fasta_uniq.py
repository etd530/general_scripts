#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Given an input FASTA file, extract unique sequences and save them to a new FASTA file.
Sequences are prioritized based on their order of appearance.

By default, sequence comparison is case-insensitive. Use --case-sensitive to
consider sequences with distinct case as different.

A second output file contains the IDs associated with each unique sequence.
Each line contains the IDs of identical sequences, separated by commas.

Usage:
    fasta_uniq.py [--case-sensitive] <input_file> <output_file> <id_output_file>
    fasta_uniq.py -h | --help

Arguments:
    <input_file>                Input FASTA file
    <output_file>               Output FASTA file with unique sequences
    <id_output_file>            Output file containing IDs of identical sequences

Options:
    --case-sensitive            Treat uppercase and lowercase sequences as different
    -h, --help                  Show this help message and exit
"""

#### LIBS ####
from Bio import SeqIO
from docopt import docopt


if __name__ == "__main__":
    #### VARS ####
    arguments = docopt(__doc__, version='fasta_uniq 1.2')

    input_file = arguments['<input_file>']
    output_file = arguments['<output_file>']
    id_output_file = arguments['<id_output_file>']
    case_sensitive = arguments['--case-sensitive']

    #### MAIN ####
    # Dictionary to store unique sequences.
    # The first record encountered for each sequence is retained.
    unique_sequences = {}

    # Dictionary to store all IDs associated with each sequence.
    # Lists preserve the order in which IDs were encountered.
    sequence_ids = {}

    # Read the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq)

        # Use an uppercase sequence as the comparison key by default.
        # With --case-sensitive, use the original sequence.
        if case_sensitive:
            seq_key = seq_str
        else:
            seq_key = seq_str.upper()

        # Store the first occurrence of each unique sequence
        if seq_key not in unique_sequences:
            unique_sequences[seq_key] = record
            sequence_ids[seq_key] = []

        # Store the ID associated with this sequence.
        # Avoid adding the same ID more than once.
        if record.id not in sequence_ids[seq_key]:
            sequence_ids[seq_key].append(record.id)

    # Write the unique sequences to the output FASTA file
    with open(output_file, "w") as output_handle:
        for record in unique_sequences.values():
            SeqIO.write(record, output_handle, "fasta")

    # Write the IDs associated with each unique sequence.
    # The order follows the order of first occurrence of each sequence.
    with open(id_output_file, "w") as id_handle:
        for seq_key in unique_sequences:
            id_handle.write(",".join(sequence_ids[seq_key]) + "\n")
