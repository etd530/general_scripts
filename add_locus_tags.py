#!/usr/bin/env python3

"""
This script adds locus tags to a GenBank file based on a provided prefix.

Usage:
	add_locus_tags.py --input=GBK --output=GBK --prefix=STR

Options:
	-i, --input	GBK	        Input GenBank file.
	-o, --output GBK		Output GenBank file with locus tags.
	-p, --prefix STR	    Prefix for locus tags.
	-h, --help              Show this help message and exit.
"""

#### LIBS ####
from docopt import docopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys

if __name__ == "__main__":
	#### VARS ####
	arguments = docopt(__doc__, version='add_locus_tags 1.0')
	infile = arguments['--input']
	outfile = arguments['--output']
	prefix = arguments['--prefix']

	#### MAIN ####
	# Read the input GenBank file
	gb_record = SeqIO.read(open(infile, "r"), "genbank")
	gb_features = gb_record.features

	# Process the genbank file, for each record, if no locus tag is available, add one based on the prefix
	feature_count = 0 # to add numbers for successive features
	for seq_feature in gb_features:
		# print(seq_feature)
		if seq_feature.type=="gene":
			feature_count += 1 # add one if we have found a new gene feature
		if seq_feature.type in ["gene", "tRNA", "CDS", "rRNA"]:
			if 'locus_tag' not in seq_feature.qualifiers:
				new_tag = prefix + "%03d" % feature_count
				# print(new_tag)
				seq_feature.qualifiers['locus_tag'] = new_tag
				# print(seq_feature.qualifiers)

			elif 'locus_tag' in seq_feature.qualifiers:
				sys.exit("ERROR: Some features already contain locus tags!")

	# write output GBK file
		with open(outfile, 'w') as fh:
			SeqIO.write(gb_record, fh, "genbank")
