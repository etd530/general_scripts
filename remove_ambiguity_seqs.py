#!/usr/bin/env python3

"""
Remove sequences containing ambiguous bases/aminoacids from a FASTA file.

Usage:
	remove_ambiguity_seqs.py --fasta=FASTA --type=STR [--keep_stops] [-h, --help]

Options:
	--fasta=FASTA                          The FASTA file with the sequences to revise.
	--type=STR                             Either 'n' (nucleotides) or 'p' (protein).
	--keep_stops                           Do NOT count stop codons as invalid characters.
	-h, --help                             Print this message help and exit.
"""

#### LIBS ####
import sys
from docopt import docopt
from Bio.SeqIO import FastaIO

#### MAIN ####
if __name__ == "__main__":
	#### VARS ####
	arguments = docopt(__doc__, version='remove_ambiguity_seqs 1.0')
	print(arguments)
	fasta = arguments['--fasta']
	seq_type = arguments['--type']
	keep_stops = arguments['--keep_stops']

	prefix = fasta.replace('.fasta', '')
	outfile = prefix + '.no_ambiguities.fasta'

	#### MAIN ####
	# Chose either nucleotide or protein alphabet
	if seq_type == 'n':
		alphabet = {'A', 'T', 'C', 'G'}
	elif seq_type == 'p':
		alphabet = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
	else:
		exit("ERROR: --type must be either 'n' for nucleotide sequences or 'p' for protein sequences.")

	# If option to keep stop codons is set, add the asterisk to the valid alphabet
	if keep_stops:
		alphabet.add('*')

	# Iterate over FASTA file
	with open(outfile, 'w') as wh:
		with open(fasta, 'r') as fh:
				for seq in FastaIO.FastaIterator(fh):
					keep_seq = True
					for char in seq.seq:
						if char not in alphabet:
							print("Invalid character '%s' found in sequence %s" % (char, seq.id))
							keep_seq = False
							break
					if keep_seq:
						wh.write('>%s\n%s\n' % (seq.id, seq.seq))
