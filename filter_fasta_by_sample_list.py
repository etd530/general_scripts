#!/usr/bin/env python3

"""
Filter samples in a FASTA file based on a list of samples to keep.

Usage:
	filter_fasta_by_sample_list.py <FASTA> <SAMPLES_LIST>
"""

#### LIBS ####
from docopt import docopt      # to create the argument parser
from Bio.SeqIO import FastaIO  # to work with fasta sequences

if __name__ == '__main__':
	#### ARGS ####
	args = docopt(__doc__)
	print(args)

	fasta = args['<FASTA>']
	samples2keep = args['<SAMPLES_LIST>']

	#### MAIN ####
	samples_list = []
	with open(samples2keep, 'r') as fh:
		for sample in fh:
			samples_list.append(sample.strip())
	with open(fasta, 'r') as fh:
		with open(fasta + '.samples_filtered.fasta', 'w') as wh:
			for sequence in FastaIO.FastaIterator(fh):
				if sequence.id in samples_list:
					wh.write('>%s\n%s\n' % (sequence.id, sequence.seq))