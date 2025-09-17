#!/usr/bin/env python3

"""
Given a txt file with species names, get the corresponding taxids.

Usage:
	get_taxids.py <TXT> <TAXDB> [-h, --help]

Arguments:
	<TXT>                     A plain text file containing species names.
	<TAXDB>                   Path to the location of the taxdb files (nodes, names, and merged)

Options:
	-h, --help                Show this help message and exit.
"""

#### LIBS ####
import taxopy                  # to work with NCBI Taxonomy database
from docopt import docopt


if __name__ == '__main__':
	#### ARGS #####
	arguments = docopt(__doc__, version='get_taxids 1.0')
	file = arguments['<TXT>']
	taxdb_path = arguments['<TAXDB>']

	#### MAIN ####
	# Create taxdb object
	taxdb_nodes = taxdb_path + '/nodes.dmp'
	taxdb_names = taxdb_path + '/names.dmp'
	taxdb_merged = taxdb_path + '/merged.dmp'
	
	print('Reading taxonomy database, please wait...')
	taxdb = taxopy.TaxDb(nodes_dmp = taxdb_nodes, names_dmp = taxdb_names, merged_dmp = taxdb_merged)

	# Iterate through the file and write the corresponding taxids
	taxid_string = ''
	with open(file, 'r') as fh:
		with open(file.replace('.txt', '') + '.taxids.txt', 'w') as wh:
			for sp in fh:
				sp = sp.strip('\n')
				taxid = taxopy.taxid_from_name(sp, taxdb)
				if len(taxid) > 1:
					print("WARNING: more than one taxid found for species %s. Printing all of them." % sp)
					taxid = ','.join(taxid)
				elif len(taxid) == 1:
					taxid = str(taxid[0])
				else:
					print('WARNING: No taxid found for species %s!!! Printing NA.\nYou should revise the NCBI taxonomy database to ensure this name is considered valid.' % sp)
					taxid = 'NA'
				taxid_string = taxid_string + taxid + '\n'
			wh.write(taxid_string)
