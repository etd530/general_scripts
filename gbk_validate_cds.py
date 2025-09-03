#!/usr/bin/env python3

"""
Given a GenBank file with annotated CDS and reference FASTA, confirm that 
the coordinates of the CDS give the expcted coding sequences.

Usage:
	gbk_validate_cds.py --gb=GBK --fasta=FASTA [-h, --help]

Options:
	--gb GBK                    GenBank flat file to revise.
	--fasta FASTA               FASTA file containing the DNA sequence from which the CDS originate.
	-h, --help                  Print this message help and exit.
"""
#### LIBS ####
from docopt import docopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

if __name__ == "__main__":
	#### VARS ####
	arguments = docopt(__doc__, version='gbk_validate_cds 1.0')
	gbk = arguments['--gb']
	fasta = arguments['--fasta']

	#### MAIN ####
	# Read the FASTA file
	nt_seq = SeqIO.read(fasta, "fasta")
	
	# Read the GenBank file
	gb_record = SeqIO.read(open(gbk,"r"), "genbank")
	gb_feature = gb_record.features

	# Iterate over features and find the CDS to evaluate
	for seq_feature in gb_feature:
		type = seq_feature.type
		if type == 'CDS':
			# get coordinates, strand, and expected AA sequence
			strand = seq_feature.location.strand
			start = seq_feature.location.start
			end = seq_feature.location.end
			expected_aa_seq = seq_feature.qualifiers['translation'][0]
			# print(expected_aa_seq)
			# if an internal stop codon is found, print a warning
			if re.search(r'\*', expected_aa_seq[:-1]):
				print("WARNING: internal stop codon found in CDS %s. Skipping this feature." % seq_feature)
			
			# extact the NT sequence from the fasta
			cds = nt_seq[start:end].seq
			if strand == -1:
				cds = cds.reverse_complement()
			if len(cds) % 3 != 0:
				print("WARNING: CDS sequence for %s is not a multiple of three. Translating in frame = 1\n         If the translation is fine, this may be an incomplete stop codon and you can ignore this warning. If unsure, leave as it is." % seq_feature.qualifiers['gene'][0])
			
			actual_aa_seq = cds.translate(table = 5)
			
			if actual_aa_seq != expected_aa_seq:
				# make sure they are equal when accounting for stop codons being in the annotated CDS but not in the translation
				if actual_aa_seq[-1] == '*' and len(actual_aa_seq) == len(expected_aa_seq) +1:
					if actual_aa_seq[:-1] == expected_aa_seq:
						continue
				# check they differ only for the first AA due to non-canonical start codons
				if actual_aa_seq[1:-1] == expected_aa_seq[1:] and actual_aa_seq[0] == 'I' and expected_aa_seq[0] == 'M':
					print("WARNING: sequences differ only in first residue (Ile vs Met)")
				else:
					print("WARNING: expected and obtained AA sequences do not match!!!")
					print("Gene name: %s" % seq_feature.qualifiers['gene'][0])
					print("Expected sequence: %s" % expected_aa_seq)
					print("Obtained sequence: %s" % actual_aa_seq)
