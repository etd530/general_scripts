#!/usr/bin/env python3

"""
Compute pariwise distances (uncorrected p-distance) for a given set of sequences in one or more multiple sequence alignments containing the same samples
and return the distance matrix for each locus as well as the mean for the concatenation of all loci.

Samples do not need to be present in all files; if a loci is missing from a sample the total length of alignment is just shorter.

WARNING: Differences are counted base on different characters, i.e. an ambiguity code (e.g. R) is considered different from A, T, C and G (SHOULD UPDATE).

Usage:
	pairwise_pdist_matrix.py <FASTA>... [--no_gaps] [--verbose] [--help]

Arguments:
	<FASTA>                      One or more multiple sequence alignments in FASTA format, containing the same samples.
	--no_gaps                    Do not count gap characters ('-') as a difference for p-distance calculation.
	--verbose                    Print the progressions of the program to the terminal (Standard Error).
	--help                       Show this help message and exit.
"""

#### LIBS ####
from Bio import AlignIO         # to work with Seq alignments
from docopt import docopt       # to create the argument parser
import pandas as pd             # to work with data frames
import matplotlib.pyplot as plt # basic plotting
import seaborn                  # nice plotting
import sys                      # exit with OK or error statuses
import numpy as np              # to work with ndarrays and do maths

#### FUNS ####
def number_of_substitutions(seq1, seq2, no_gaps = False):
    assert len(seq1) == len(seq2)
    if no_gaps:
    	return sum(1 if nuc1 != nuc2 and nuc1 != '-' and nuc2 != '-' else 0 for(nuc1, nuc2) in zip(seq1, seq2))
    else:
    	return sum(0 if nuc1 == nuc2 else 1 for (nuc1, nuc2) in zip(seq1, seq2))

def ali_len_msa2pairwise(seq1, seq2, no_gaps = False):
	assert len(seq1) == len(seq2)
	if no_gaps:
		return sum(0 if nuc1 == '-' or nuc2 == '-' else 1 for (nuc1, nuc2) in zip(seq1, seq2))
	else:
		return sum(0 if nuc1 == '-' and nuc2 == '-' else 1 for (nuc1, nuc2) in zip(seq1, seq2))

def plot_heatmap(dataframe, title, outfile, outformat = "svg", log_scale = False):
    nrow = len(dataframe.index)
    ncol = len(dataframe.columns)
    figuresize = [30+ncol, 20+nrow]
    plt.figure(figsize = figuresize)
    plt.title(title, fontsize = 120)
    # cmap = seaborn.cm.viridis
    cmap = "viridis"
    if log_scale:
        ax = heatmap = seaborn.heatmap(dataframe, square=True, fmt = '.2g', annot=True, norm = LogNorm(clip = True), cmap = cmap)
    else:
        ax = heatmap = seaborn.heatmap(dataframe, square=True, fmt = '.2g', annot=True, cmap = cmap) # vmin = 0, vmax = 1)
    plt.savefig(outfile, format = outformat)
    plt.close()

if __name__ == '__main__':
	#### PARSE ARGS ####
	args = docopt(__doc__)
	no_gaps = args['--no_gaps']
	if args['--verbose']:
		if no_gaps:
			print('Computing distances without considering gaps...')
		else:
			print('Computing distances considering gaps...')

	#### GRAPHICAL PARAMS ####
	seaborn.set(font_scale = 1.25)

	#### MAIN ####
	sequences_set = []

	for file in args['<FASTA>']:
		ali = AlignIO.read(file, "fasta")
		sequence_ids = [seq.id for seq in ali]
		sequences_set = sequences_set + sequence_ids

	sequences_set = list(set(sequences_set))
	sequences_set.sort()

	subs_num_per_locus = []
	ali_len_per_locus = []
	# concatenated_len = 0

	for file in args['<FASTA>']:
		print(file)
		ali = AlignIO.read(file, "fasta")
		# concatenated_len += ali.get_alignment_length()
		pairwise_subs_num = pd.DataFrame(index = sequences_set, columns = sequences_set, dtype = 'float')
		# pairwise_subs_num = pd.DataFrame(index = sequences_set[1:], columns = sequences_set[:-1], dtype = 'float')
		# pairwise_ali_len = pd.DataFrame(index = sequences_set[1:], columns = sequences_set[:-1], dtype = 'float')
		pairwise_ali_len = pd.DataFrame(index = sequences_set, columns = sequences_set, dtype = 'float')
		nb_seq = len(sequences_set)
		for i in range(nb_seq):
			# print('Looking for sequence: %s' % sequences_set[i])
			seq1_found = False
			for seq in ali:
				if seq.id == sequences_set[i]:
					seq1_found = True
					# print('Found sequence %s' % seq.id)
					rec1 = seq
					break
			if seq1_found:
				j = i + 1
				while j < nb_seq:
					# print('Looking for sequence: %s' % sequences_set[j])
					seq2_found = False
					for seq in ali:
						if seq.id == sequences_set[j]:
							seq2_found = True
							# print('Found sequence %s' % seq.id)
							rec2 = seq
							break
					if seq2_found:
						# try:
						# print(rec1.id)
						# print(len(rec1))
						# print(rec2.id)
						# print(len(rec2))
						subs_num = number_of_substitutions(rec1.seq, rec2.seq, no_gaps = no_gaps)
						pairwise_subs_num[rec1.id][rec2.id] = subs_num
						pairwise_subs_num[rec2.id][rec1.id] = subs_num
						# print('Number of substitutions between %s and %s in loci %s: %s' % (rec1.id, rec2.id, file, str(number_of_substitutions(rec1.seq, rec2.seq))))
						# print('Length of alignment between %s and %s in loci %s: %s' % (rec1.id, rec2.id, file, str(len(rec1.seq))))
						ali_len = ali_len_msa2pairwise(rec1.seq, rec2.seq, no_gaps = no_gaps)
						pairwise_ali_len[rec1.id][rec2.id] = ali_len
						pairwise_ali_len[rec2.id][rec1.id] = ali_len
						# except KeyError:
						# 	pairwise_subs_num[rec2.id][rec1.id] = number_of_substitutions(rec1.seq, rec2.seq)
					j += 1
		pairwise_dist_matrix = pairwise_subs_num/len(ali[0])
		if no_gaps:
			pairwise_dist_matrix.to_csv(file.replace('.fasta', '') + '.pdist_matrix.no_gaps.tsv', sep = '\t')
		else:
			pairwise_dist_matrix.to_csv(file.replace('.fasta', '') + '.pdist_matrix.tsv', sep = '\t')
		subs_num_per_locus.append(pairwise_subs_num.to_numpy())
		ali_len_per_locus.append(pairwise_ali_len.to_numpy())
		if no_gaps:
			plot_heatmap(pairwise_dist_matrix, title = 'Pairwise p-distances', outfile = file.replace('.fasta', '') + '.pdist_matrix.no_gaps.png', outformat = 'png')
			plot_heatmap(pairwise_dist_matrix, title = 'Pairwise p-distances', outfile = file.replace('.fasta', '') + '.pdist_matrix.no_gaps.pdf', outformat = 'pdf')
		else:
			plot_heatmap(pairwise_dist_matrix, title = 'Pairwise p-distances', outfile = file.replace('.fasta', '') + '.pdist_matrix.png', outformat = 'png')
			plot_heatmap(pairwise_dist_matrix, title = 'Pairwise p-distances', outfile = file.replace('.fasta', '') + '.pdist_matrix.pdf', outformat = 'pdf')

	if args['--verbose']:
		print('Finished copmutation for each locus, computing mean for concatenation of all loci...')

	subs_num_per_locus = np.array(subs_num_per_locus)
	ali_len_per_locus = np.array(ali_len_per_locus)
	subs_num_total = np.nansum(subs_num_per_locus, axis = 0)
	ali_len_total = np.nansum(ali_len_per_locus, axis = 0)
	ali_len_total[ali_len_total == 0] = np.nan
	# print(pd.DataFrame(data = subs_num_total, index = pairwise_subs_num.index, columns = pairwise_subs_num.columns))
	# print(pd.DataFrame(data = ali_len_total, index = pairwise_subs_num.index, columns = pairwise_subs_num.columns))
	pairwise_dist_total_matrix = subs_num_total/ali_len_total
	# pairwise_dist_total_matrix[np.triu_indices(pairwise_dist_total_matrix.shape[0], 1)] = np.nan
	pairwise_dist_total_matrix = pd.DataFrame(data = pairwise_dist_total_matrix, index = pairwise_subs_num.index, columns = pairwise_subs_num.columns)
	# print(pairwise_dist_total_matrix)
	pairwise_dist_total_matrix.to_csv('concatenated_loci.pdist_matrix.tsv', sep = '\t')
	if no_gaps:
		plot_heatmap(pairwise_dist_total_matrix, title = 'Pairwise p-distances', outfile = 'concatenated_loci.pdist_matrix.no_gaps.png', outformat = 'png')
		plot_heatmap(pairwise_dist_total_matrix, title = 'Pairwise p-distances', outfile = 'concatenated_loci.pdist_matrix.no_gaps.pdf', outformat = 'pdf')
	else:
		plot_heatmap(pairwise_dist_total_matrix, title = 'Pairwise p-distances', outfile = 'concatenated_loci.pdist_matrix.png', outformat = 'png')
		plot_heatmap(pairwise_dist_total_matrix, title = 'Pairwise p-distances', outfile = 'concatenated_loci.pdist_matrix.pdf', outformat = 'pdf')
	if args['--verbose']:
		print('Execution finished.')
	sys.exit()
