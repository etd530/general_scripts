#!/usr/bin/env python3

"""
Given a FASTA file of unaligned sequences, compute basic statistics of the sequences contained in it.

Usage:
	fasta_seqstats.py <FASTA> [--verbose] [--help]

Arguments:
	<FASTA>                      A FASTA file containing unaligned sequences.
	--verbose                    Print the progressions of the program to the terminal (Standard Error).
	--help                       Show this help message and exit.
"""

#### LIBS ####
from Bio import SeqIO           # to work with sequence files7
from docopt import docopt       # to create the argument parser
import matplotlib.pyplot as plt # basic plotting
import seaborn                  # nice plotting
import sys                      # exit with OK or error statuses
import numpy as np              # to work with ndarrays and do maths

#### FUNS ####

if __name__ == '__main__':
	#### PARSE ARGS ####
	args = docopt(__doc__)
	fasta_file = args['<FASTA>']

	fasta_prefix = fasta_file.split('.')[0].split('/')[-1]

	#### GRAPHICAL PARAMS ####
	# seaborn.set(font_scale = 1.25)

	#### MAIN ####
	# Iterate over the fasta file
	print("Reading file %s..." % fasta_file)
	
	len_seqs = []
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		len_seqs.append(len(seq_record))
	
	print("Done reading file %s. Computing statistics..." % fasta_file)
	# Count the number of sequences in the file
	nb_seq = len(len_seqs)
	if args['--verbose']:
		print("Number of sequences: %s" % nb_seq)
	
	# Get the mean length of the sequences
	mean_len_seqs = np.mean(len_seqs)
	if args['--verbose']:
		print("Mean length of sequences: %s" % mean_len_seqs)
	# Get the standard deviation of the length of the sequences
	std_len_seqs = np.std(len_seqs)
	if args['--verbose']:
		print("Standard deviation of the length of sequences: %s" % std_len_seqs)
	# Get the median length of the sequences
	median_len_seqs = np.median(len_seqs)
	if args['--verbose']:
		print("Median length of sequences: %s" % median_len_seqs)
	# Get the minimum and maximum length of the sequences
	min_len_seqs = min(len_seqs)
	max_len_seqs = max(len_seqs)
	if args['--verbose']:
		print("Minimum length of sequences: %s" % min_len_seqs)
		print("Maximum length of sequences: %s" % max_len_seqs)
	# Get the 5% and 95% quantiles of the length of the sequences
	q5_len_seqs = np.percentile(len_seqs, 5)
	q95_len_seqs = np.percentile(len_seqs, 95)
	if args['--verbose']:
		print("5% quantile of the length of sequences: %s" % q5_len_seqs)
		print("95% quantile of the length of sequences: %s" % q95_len_seqs)
	# Generate log file with the statistics of the sequences
	with open(fasta_prefix + '_seqstats.log', 'w') as log_file:
		write_string = "Statistics of the sequences contained in the file %s:\n" % fasta_file
		write_string += "---------------------------------------------\n"
		write_string += "Number of sequences: %s\n" % nb_seq
		write_string += "Mean length of sequences: %s\n" % mean_len_seqs
		write_string += "Standard deviation of the length of sequences: %s\n" % std_len_seqs
		write_string += "Median length of sequences: %s\n" % median_len_seqs
		write_string += "Minimum length of sequences: %s\n" % min_len_seqs
		write_string += "Maximum length of sequences: %s\n" % max_len_seqs
		# write_string += "5% quantile of the length of sequences: %s\n" % q5_len_seqs
		# write_string += "95% quantile of the length of sequences: %s\n" % q95_len_seqs
		log_file.write(write_string)
	# Generate a histogram of the length of the sequences
	plt.figure(figsize = (10, 6))
	seaborn.histplot(len_seqs, bins = 30, kde = True)
	plt.title("Distribution of the length of sequences in the file %s" % fasta_file)
	plt.xlabel("Length")
	plt.ylabel("Frequency")
	plt.savefig(fasta_prefix + '_seqstats.png', dpi = 300, bbox_inches = 'tight')

	if args['--verbose']:
		print('Execution finished.')
	sys.exit()
