#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given a set of proteomes in FASTA format and a set of genes, output one file containing all homologs of each gene from each proteome.

Usage:
	proteomes2homologs.py <genes> <proteomes> <proteomes> ...
	proteomes2homologs.py -h | --help

Arguments:
	<genes>                               List of gene names separated by commas
	<proteomes> <proteomes> ...           Two or more proteome files in FASTA format

Options:
	-h, --help                            Show this help message and exit
"""

#### LIBS ####
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from docopt import docopt

if __name__ == "__main__":
	#### VARS ####
	arguments = docopt(__doc__, version='proteomes2homologs 1.0')
	proteome_files = arguments['<proteomes>']
	gene_list = arguments['<genes>'].split(',')

	#### MAIN ####
	# Make a dictionary where keys will be the target genes
	genes_dict = {}
	
	# Process each proteome file and sort all genes by their gene name
	for proteome_file in proteome_files:
		proteome_record = SeqIO.parse(open(proteome_file, "r"), "fasta")
		for seq_record in proteome_record:
			gene_name = seq_record.id.split('.')[-1]
			if gene_name in genes_dict.keys():
				genes_dict[gene_name].append(seq_record)
			else:
				genes_dict[gene_name] = [seq_record]

	# Write one file for each gene
	for gene in genes_dict.keys():
		if gene in gene_list:
			with open(f"{gene}.faa", "w") as output_handle:
				for seq_record in genes_dict[gene]:
					output_handle.write(f'>{seq_record.id}\n{str(seq_record.seq)}\n')
