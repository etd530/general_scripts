#!/usr/bin/env python3

"""
Given a TSV file output from BLAST, get the MRCA of the taxids of the hits. Optionally, provide a taxid to restrict the search.

Usage:
    blast_get_mrca.py <BLAST_TSV> <TAXDB> --colnum=<INT> [--taxid=<TAXID>] [-h, --help]

Arguments:
    <BLAST_TSV>               A TSV file output from BLAST with taxids. Query ID is expected in column 0.
    <TAXDB>                   Path to the location of the taxdb files (nodes, names, and merged)

Options:
    --colnum=<INT>            Column number (0-based) in the BLAST TSV file where taxids are located.
    --taxid=<TAXID>           Taxid to restrict the MRCA search to a specific clade. [default: 1]
    -h, --help                Show this help message and exit.
"""

#### LIBS ####
import taxopy                  # to work with NCBI Taxonomy database
import sys
import pandas as pd
from docopt import docopt

if __name__ == '__main__':
	#### ARGS #####
	arguments = docopt(__doc__, version='blast_get_mrca 1.0')
	blast_tsv = arguments['<BLAST_TSV>']
	taxdb_path = arguments['<TAXDB>']
	colnum = int(arguments['--colnum'])
	restrict_taxid = int(arguments['--taxid'])

	#### MAIN ####
	# Create taxdb object
	taxdb_nodes = taxdb_path + '/nodes.dmp'
	taxdb_names = taxdb_path + '/names.dmp'
	taxdb_merged = taxdb_path + '/merged.dmp'
	
	print('Reading taxonomy database, please wait...')
	taxdb = taxopy.TaxDb(nodes_dmp = taxdb_nodes, names_dmp = taxdb_names, merged_dmp = taxdb_merged)

	# Read BLAST TSV as Pandas dataframe
	df = pd.read_csv(blast_tsv, sep='\t', header=None)
	
	# Get list of query IDs
	query_ids = df[0].unique().tolist()

	# For each query, find MRCA of taxids in specified column, restricted to given taxid
	out_string = "query\tMRCA_name\n"
	for query_id in query_ids:
		# print(f'Processing query ID: {query_id}')
		query_df = df[df[0] == query_id]

		# Collect taxids from the specified column from the subset dataframe
		taxa_set = set()
		for line in query_df.itertuples(index=False):
			taxid_list = str(line[colnum]).split(';')  # in case of multiple taxids separated by semicolon
			for taxid_str in taxid_list:
				try:
					taxid = int(taxid_str)
					taxon = taxopy.Taxon(taxid, taxdb)
				except ValueError:
					sys.exit(f"ERROR: Invalid taxid '{taxid_str}' in column {colnum}. Exiting.")
				except taxopy.exceptions.TaxidError:
					print("The input integer is not a valid NCBI taxonomic identifier: %s. Please try update the taxdump files." % taxid_str)
				if restrict_taxid in taxon.taxid_lineage:
					taxa_set.add(taxon)

		# Get MRCA
		if taxa_set:
			if len(taxa_set) > 1: # function to find LCA requires at least 2 taxa
				mrca_taxon = list(taxa_set)[0]
				mrca_taxon = taxopy.find_lca(list(taxa_set), taxdb)
				mrca_name = mrca_taxon.name
			elif len(taxa_set) == 1: # if there is only one taxon, that is the MRCA
				mrca_taxon = list(taxa_set)[0]
				mrca_name = mrca_taxon.name
			out_string = out_string + f'{query_id}\t{mrca_name}\n'
		else:
			print('No valid taxids found in the BLAST TSV file, or no hits belonging to the restricted taxonomic search.')
		
		# Write output to file
		with open(f'{blast_tsv}.mrcas.tsv', 'w') as out_fh:
			out_fh.write(out_string)
