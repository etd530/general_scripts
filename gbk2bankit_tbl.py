#!/usr/bin/env python3

"""
Convert a GenBank file to the table format required for submission to BankIt.

Usage:
	gff2bankit_tbl.py <GBK> [-h, --help]

Arguments:
	<GBK>         An input GBK file.

Options:
	-h, --help    Print this message help and exit
"""
#### LIBS ####
from docopt import docopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

if __name__ == "__main__":
	#### VARS ####
	arguments = docopt(__doc__, version='gbk2bankit_tbl 1.0')
	infile = arguments['<GBK>']
	file_prefix = infile.replace('.gb', '')
	outfile =  file_prefix + '.tbl'

	#### MAIN ####
	# Read the GenBank file
	gb_record = SeqIO.read(open(infile,"r"), "genbank")
	gb_feature = gb_record.features

	# Extract the locus name from the GenBank file
	locus = gb_record.name

	with open(outfile, 'w') as fh:
		fh.write('>Feature %s\n' % locus)
		for seq_feature in gb_feature:
			type = seq_feature.type
			strand = seq_feature.location.strand
			if strand == 1:
				start = seq_feature.location.start + 1 # we add one because bankit files are 1-based
				end = seq_feature.location.end
			else:
				start = seq_feature.location.end
				end = seq_feature.location.start + 1 # we add one because bankit files are 1-based
			# Check if it's a CDS and if so, if it has internal stop codons
			if type == 'CDS':
				aa_seq = seq_feature.qualifiers['translation'][0]
				# if an internal stop codon is found, continue to the next feature
				if re.search(r'\*', aa_seq[:-1]):
					print("Warning: internal stop codon found in CDS %s. Skipping this feature." % seq_feature)
					continue

			fh.write("%s\t%s\t%s\n" % (start, end, type))
			for qualifier in seq_feature.qualifiers:
				name = qualifier
				value = seq_feature.qualifiers[name][0]
				# we change the qualifier "gene" to "product" for the CDS to comply with the format
				if name == 'gene' and seq_feature.type == 'CDS':
					name = 'product'
				# check that tRNA product qualifiers don't end in a number, else fix it
				if name == 'product' and seq_feature.type in ['tRNA', 'gene']:
					if re.search(r'\d$', value):
						print("Warning: tRNA/gene product qualifier %s ends in a number. Fixing it by removing the number." % value)
						value = re.sub(r'\d$', '', value)
				# write the qualifiers except translation, which will be processed by GenBank
				if name != 'translation':
					fh.write("\t\t\t%s\t%s\n" % (name, value))
