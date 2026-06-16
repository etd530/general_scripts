#!/usr/bin/env python3

"""
Given a TSV file output from BLAST, get summary statistics of the taxonomy of the hits.
You can get the MRCA of the taxids of the hits. Optionally, provide a taxid to restrict the search.
You can also compute an alien index score for belonging to a given taxid (log10(best belonging eval) - log10(best nonbelonging eval)).

Usage:
    blast_taxonomy.py <BLAST_TSV> <TAXDB> --eval_colnum=<INT> --taxids_colnum=<INT> [--taxid=<TAXID> --alien-index=<TAXID>] [-h, --help]

Arguments:
    <BLAST_TSV>               A TSV file output from BLAST with taxids. Query ID is expected in column 0.
    <TAXDB>                   Path to the location of the taxdb files (nodes, names, and merged)

Options:
    --taxids_colnum=<INT>     Column number (0-based) in the BLAST TSV file where taxids are located.
    --taxid=<TAXID>           Taxid to restrict the MRCA search to a specific clade. [default: 1]
    --eval_colnum=<INT>       Column number (0-based) in the BLAST TSV file where e-values are located, for alien index computation.
    --alien-index=<TAXID>     Taxid to use to compute an alien index score for belonging to that taxon. [default: 1]
    -h, --help                Show this help message and exit.
"""

#### LIBS ####
import taxopy                  # to work with NCBI Taxonomy database
import sys
import pandas as pd
from docopt import docopt
import numpy as np

#### FUNS ####
def alien_index(best_belonging_eval, best_nonbelonging_eval):
	if best_belonging_eval == float('inf') and best_nonbelonging_eval == float('inf'):
		return np.nan  # if no hits at all, return NA
	elif best_belonging_eval == float('inf'):
		return -np.inf  # if no belonging hits, AI is negative infinity
	elif best_nonbelonging_eval == float('inf'):
		return np.inf  # if no non-belonging hits, AI is positive infinity
	else:
		return np.log10(best_nonbelonging_eval) - np.log10(best_belonging_eval) # else we do the actual number

if __name__ == '__main__':
	#### ARGS #####
	arguments = docopt(__doc__, version='blast_taxonomy 1.0')
	blast_tsv = arguments['<BLAST_TSV>']
	taxdb_path = arguments['<TAXDB>']
	taxids_colnum = int(arguments['--taxids_colnum'])
	restrict_taxid = int(arguments['--taxid'])
	alien_index_taxid = int(arguments['--alien-index'])
	eval_colnum = int(arguments['--eval_colnum'])
	
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

	# Prepare header line for output file
	if alien_index_taxid != 1:
		out_string = "query\tMRCA_name\tMRCA_lineage\tAlien_index_score\n"
	else:
		out_string = "query\tMRCA_name\tMRCA_lineage\n"

	# Iterate over each query to find the MRCA of taxids in specified column, restricted to given taxid
	# Also compute alien index if requested, by tracking the best e-value of belonging and non-belonging hits to the specified taxid
	for query_id in query_ids:         # for each query sequence in the input file
		# print(f'Processing query ID: {query_id}')
		query_df = df[df[0] == query_id]
		min_belonging_eval = float('inf')
		min_nonbelonging_eval = float('inf')
		
		# Collect taxids from the specified column from the subset dataframe
		taxa_set = set()
		for line in query_df.itertuples(index=False):
			taxid_list = str(line[taxids_colnum]).split(';')  # in case of multiple taxids separated by semicolon
			eval = line[eval_colnum]
			for taxid_str in taxid_list:
				# get the taxid and add it to the set of taxa if it belongs to the restricted taxid
				try:
					taxid = int(taxid_str)
					# print(f'Processing taxid: {taxid} with e-value: {eval}')
					taxon = taxopy.Taxon(taxid, taxdb)
					# print(f'Taxon name: {taxon.name}')
				except ValueError:
					print(f"WARNING: Invalid taxid '{taxid_str}' for query {query_id}. Please make sure you provided the correct column for the taxids!")
				except taxopy.exceptions.TaxidError:
					print("The input integer is not a valid NCBI taxonomic identifier: %s. Please try update the taxdump files." % taxid_str)
				if restrict_taxid in taxon.taxid_lineage:
					# print("Adding taxid %s to the set of taxa for query %s" % (taxid, query_id))
					taxa_set.add(taxon)
					# print(taxa_set)
				# if alien index is requested, update the best belonging and non-belonging e-values
				if alien_index_taxid != 1:
					if alien_index_taxid in taxon.taxid_lineage:
						if eval < min_belonging_eval:
							min_belonging_eval = eval
					else:
						if eval < min_nonbelonging_eval:
							min_nonbelonging_eval = eval

		# Compute alien index if requested
		if alien_index_taxid != 1:
			ai_score = alien_index(min_belonging_eval, min_nonbelonging_eval)
			# print(ai_score)

		# print(taxa_set)
		# Get MRCA
		if taxa_set:
			if len(taxa_set) > 1: # function to find LCA requires at least 2 taxa
				mrca_taxon = list(taxa_set)[0]
				mrca_taxon = taxopy.find_lca(list(taxa_set), taxdb)
				mrca_name = mrca_taxon.name
				mrca_lineage = mrca_taxon.name_lineage
				print(f'MRCA of taxids for query {query_id} is {mrca_name} (taxid {mrca_taxon.taxid})')
			elif len(taxa_set) == 1: # if there is only one taxon, that is the MRCA
				mrca_taxon = list(taxa_set)[0]
				mrca_name = mrca_taxon.name
				mrca_lineage = mrca_taxon.name_lineage
				print(f'Only one taxid found for query {query_id}, so MRCA is {mrca_name} (taxid {mrca_taxon.taxid})')
			if alien_index_taxid != 1:
				out_string = out_string + f'{query_id}\t{mrca_name}\t{mrca_lineage}\t{ai_score}\n'
			else:
				out_string = out_string + f'{query_id}\t{mrca_name}\t{mrca_lineage}\n'
		else:
			print('No valid taxids found in the BLAST TSV file, or no hits belonging to the restricted taxonomic search.')
			if alien_index_taxid != 1:
				out_string = out_string + f'{query_id}\tNo_valid_taxids_found\tNo_valid_taxids_found\t{ai_score}\n'
			else:
				out_string = out_string + f'{query_id}\tNo_valid_taxids_found\tNo_valid_taxids_found\n'
		
	# Write output to file
	with open(f'{blast_tsv}.taxonomy.tsv', 'w') as out_fh:
		out_fh.write(out_string)
