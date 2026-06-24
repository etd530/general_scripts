#!/usr/bin/env python3

"""
Estimate a viral genome's completeness based on its length relative to that of a set of reference genomes of known affinity.

Usage:
    checkv_completeness.py <queries.fasta> <diamond.tsv> <db_lens.tsv> [-h, --help]

Arguments:
    <queries.fasta>            Input FASTA with a series of query contigs.
    <diamond.tsv>              Input TSV with diamond tabular output.
    <db_lens.tsv>              Input TSV with reference genomes in the first column and their length (in base pairs) in the second column.

Options:
    -h, --help                Show this help message and exit.
"""

#### LIBS ####
import docopt
import time
import pandas as pd
import math
from Bio import SeqIO

#### FUNS ####
def expected_len(weights, lengths):
    expected_length = sum(l * w for l, w in zip(lengths, weights)) / sum(weights)
    return expected_length

if __name__ == "__main__":
    
    #### ARGS ####
    args = docopt.docopt(__doc__, version="CheckV v0.9.0")
    query_file = args["<queries.fasta>"]
    diamond_file = args["<diamond.tsv>"]
    db_file = args["<db_lens.tsv>"]

    #### MAIN ####
    # Read diamond output as pandas df
    print("Reading diamond output...")
    diamond_df = pd.read_csv(diamond_file, sep='\t', header=None)
    diamond_df.rename(columns = {0: 'qseqid', 1: 'qlen', 2: 'sallseqid', 3: 'slen', 4: 'pident', 5: 'length', 6: 'mismatch', 7: 'gapopen', 8: 'qstart', 9: 'qend', 10: 'sstart', 11: 'send', 12: 'evalue', 13: 'bitscore', 14: 'stitle'}, inplace = True)
    # print(diamond_df["qseqid"][0])

    # Remove protein index from qseqid since we will use the full contigs
    diamond_df['qseqid'] = diamond_df['qseqid'].apply(lambda x: '_'.join(x.split('_')[:-1]))

    # Turn subject IDs to lowercase
    diamond_df['sallseqid'] = diamond_df['sallseqid'].apply(lambda x: x.lower())
    print(diamond_df['sallseqid'])

    # Read reference sequences length dataframe
    print("Reading TSV of reference sequence lengths...")
    db_lens_df = pd.read_csv(db_file, sep='\t', header=None)
    db_lens_df.rename(columns = {0: 'virus', 1: 'length'}, inplace = True)

    # Iterate over queries and make a dictionary of their lengths
    print("Extracting the length of each query sequence...")
    queries_dict = {}
    for seq_record in SeqIO.parse(query_file, "fasta"):
        queries_dict[seq_record.id] = len(seq_record)
    
    # Iterate over queries dict
    outstring="query\tlength\tlength_estimate\tcompletness\n"
    print("Iterating over the Diamond results for each query...")
    for query, query_len in queries_dict.items():
        # Get subset of diamond results from that query
        subset_df = diamond_df[diamond_df['qseqid'] == query]
        
        # If the subset is not empty, then iterate over the entries and get the hit names and bitscores
        if len(subset_df.index) != 0:
            bitscores_dict = {}
            bitscore_sum = 0
            for index, row in subset_df.iterrows():
                subject = row['sallseqid']
                bitscore = float(row['bitscore'])
                
                # If the key already exists but the new value is higher, add it to the total, substract to old one, and replace it in the dictionary
                if subject in bitscores_dict.keys() and bitscore > bitscores_dict[subject]:
                    bitscore_sum += bitscore
                    bitscore_sum = bitscore_sum - bitscores_dict[subject]
                    bitscores_dict[subject] = bitscore
                # else if it does not exist, simply add (in the case it exists but it is better than the new one, we do nothing so we only keept the best hit for each ref sequence)
                elif subject not in bitscores_dict.keys():
                    bitscores_dict[subject] = bitscore
                    bitscore_sum += bitscore

            # Get the genome lengths of the query hits into a list
            # At the same time make list of bitscores, keeping only those of viruses present in the lenght df
            bitscores_list = []
            virlens_list = []
            for key in bitscores_dict.keys():
                try:
                    virlen = db_lens_df['length'][db_lens_df['virus'] == key].values[0]
                    virlens_list.append(virlen)

                    bitscore = bitscores_dict[key]
                    bitscores_list.append(bitscore)

                except IndexError:
                    print("Error when getting length for virus %s, it will be ignored for completeness calculation." % key)
            
            # Transform the bitscores to weights by dividing each by the total
            weights_list = [x/sum(bitscores_list) for x in bitscores_list]

            if len(weights_list) == 0 and len(virlens_list) == 0:
                print("WARNING: No virus hits with available lengths for query %s; no completeness estimate can be obtained for it!" % query)
                outstring = outstring + '%s\t%s\tNA\tNA\n' % (query, query_len)
                continue

            try:
                assert math.isclose(1, sum(weights_list), rel_tol=1e-6) # need to use math.isclose since the weights are low and suffer from float point error
            except AssertionError:
                print(weights_list)
                print(sum(weights_list))
                print(bitscore_sum)
                print(bitscores_dict.values())
                print(sum(bitscores_dict.values()))
                exit()
            
            # Compute the expected length and the completness estimate of the query
            len_estimate = expected_len(weights = weights_list, lengths = virlens_list)
            completeness = query_len/len_estimate
            outstring = outstring + '%s\t%s\t%s\t%s\n' % (query, query_len, len_estimate, completeness)

    # Write the output to a file
    with open('completness_estimates.tsv', 'w') as fh:
        fh.write(outstring)
    
    print("Execution done.")