#!/usr/bin/env python3

"""
Transform NUCmer alignment output in TSV format to GFF3 format.

Usage:
	nucmer2gff.py <nucmer_tsv>
	nucmer2gff.py (-h | --help)

Arguments:
	nucmer_tsv    Path to the NUCmer TSV output file.

Options:
	-h --help     Show this help message and exit.
"""

#### LIBS ####
from docopt import docopt

if __name__ == '__main__':
	
	#### VARS ####
	args = docopt(__doc__)

	nucmer_tsv = args['<nucmer_tsv>']
	nucmer_gff = nucmer_tsv.replace('.tsv', '.gff')

	#### MAIN ###
	gff_file = '##gff-version 3\n'
	with open(nucmer_tsv) as fh:
		i = 0
		for line in fh:
			if not line.startswith('S1'):
				i+=1
				line_list = line.strip().split()
				chromosome = line_list[7]
				start = line_list[0]
				end = line_list[1]
				if start < end:
					strand = '+'
				else:
					strand = '-'
				attributes = 'nuwt_%s' % i
				gff_line='\t'.join([chromosome, '.', 'NUCmer_NUWT', str(start), str(end), '.', strand, '.', attributes])
				gff_file=gff_file + gff_line + '\n'

	with open(nucmer_gff, 'w') as fh:
		fh.write(gff_file)
