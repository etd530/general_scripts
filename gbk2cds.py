#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract CDS nucleotide sequences from GenBank files and the reference genome and convert them to FASTA format (nucleotide).
Modified from original script in Stack Exchange by user pippo1980.

Usage:
	gbk2cds.py <gbk_file> <ref_genome>
	gbk2cds.py -h | --help
	 
Arguments:
	<gbk_file>                 GenBank file to be converted
	<ref_genome>               Reference genome in FASTA format
	 
Options:
	-h, --help                 Show this help message and exit
"""

#### LIBS ####
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from docopt import docopt


if __name__ == "__main__":
	#### VARS ####
	arguments = docopt(__doc__, version='gbk2faa 1.0')
	file_name = arguments['<gbk_file>']

	#### MAIN ####
	# List to store all the CDS entries
	all_entries = []

	# Read the reference genome
	ref_genome = SeqIO.read(open(arguments['<ref_genome>'], "r"), "fasta")
	
	# Read the GenBank file
	gb_record = SeqIO.read(open(file_name,"r"), "genbank")
	gb_feature = gb_record.features

	# Extract the accession number of the GenBank file
	accession = gb_record.annotations['accessions'][0]

	# Process the features in the GenBank file and extract the coordinates of the CDS
	for seq_feature in gb_feature:
		if seq_feature.type=="CDS":
				if 'translation' in seq_feature.qualifiers:
					if seq_feature.qualifiers['translation'][0] != '':
						# extract the gene name
						if 'gene' in seq_feature.qualifiers:
							gene = seq_feature.qualifiers['gene'][0]
					
						else:
							gene = 'not defined'
						
						# extract the product name
						if 'product' in seq_feature.qualifiers:
							product = seq_feature.qualifiers['product'][0]
						else:
							product = 'not defined'
						
						# extract the protein id
						if 'protein_id' in seq_feature.qualifiers:
							protein_id = seq_feature.qualifiers['protein_id'][0]
						else:
							protein_id = 'not defined'
							
						# extract the old locus tag
						if 'old_locus_tag' in seq_feature.qualifiers:
							old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
						else:
							old_locus_tag = 'not defined'
						
						# extract the locus tag
						if 'locus_tag' in seq_feature.qualifiers:
							locus_tag = seq_feature.qualifiers['locus_tag'][0]
						else:
							locus_tag = 'not defined'
				
						# extract the strand
						if seq_feature.location.strand == 1:
							complement = 'no complement'

						if seq_feature.location.strand == -1:
							complement = str(seq_feature.location).strip('[]').split('(-)')[0]

						# extract the start and end coordinates of the CDS feature
						start = seq_feature.location.start -1 # GenBank coordinates are 1-based, so we subtract 1 because Python uses 0-based indexing
						end = seq_feature.location.end

						# Extract the nucleotide sequence from the reference genome
						if seq_feature.location.strand == 1:
							nucleotide_seq = ref_genome.seq[start:end]
						else:
							nucleotide_seq = ref_genome.seq[start:end].reverse_complement()
						
						pippo = SeqRecord(
									Seq(nucleotide_seq),
									id = accession + '.' + gene,
									description = ('')
									)
						all_entries.append(pippo)

	# write file with single-line FASTA sequences
	with open(f'{file_name[:-3]}.cds.fasta', 'w') as fasta_file:
		for entry in all_entries:
			fasta_file.write(f'>{entry.id}\n{str(entry.seq)}\n')
