#!/usr/bin/env python3
"""
Subset FASTA headers to keep only some fields separated by a given delimiter.

Usage:
	subset_fasta_headers.py --delim=STR --fields=STR [--threads=INT] [--verbose] [--help] <FASTA>...

Arguments:
	<FASTA>                      One or more FASTA files to subset headers from.
	--delim=STR                  Delimiter used to separate fields in the header (default: '|').
	--fields=STR                 Comma-separated list of fields to keep in the header (default: '1').
	--threads=INT                Number of threads to use (default: 1).
	--verbose                    Print the progressions of the program to the terminal (Standard Error).
	--help                       Show this help message and exit.
"""

#### LIBS ####
from Bio import SeqIO
from docopt import docopt
from concurrent.futures import ThreadPoolExecutor

#### FUNS ####
def subset_fasta_headers(fasta_file, delimiter, fields, threads=1, verbose=False):
	"""
	Subset FASTA headers to keep only some fields separated by a given delimiter.

	Arguments:
		fasta_files: List of FASTA files to subset headers from.
		delimiter: Delimiter used to separate fields in the header.
		fields: Comma-separated list of fields to keep in the header.
		threads: Number of threads to use (default: 1).
		verbose: Print the progressions of the program to the terminal (Standard Error).
	"""
	prefix = fasta_file.split('/')[-1].rstrip('.fasta')
	with open(fasta_file, 'r') as input_handle, open(f"{prefix}.renamed.fasta", 'w') as output_handle:
		for record in SeqIO.parse(input_handle, "fasta"):
			header_fields = record.id.split(delimiter)
			try:
				subset_header = delimiter.join([header_fields[i] for i in fields])
			except IndexError:
				human_readable_fields = [str(i+1) for i in fields]
				print(f"ERROR: One or more fields {human_readable_fields} are out of range for header '{record.id}' in file '{fasta_file}'.")
				exit(1)
			record.id = subset_header
			record.description = ""
			SeqIO.write(record, output_handle, "fasta")

	if verbose:
		print(f"Subsetted headers written to {fasta_file}.renamed.fasta")

#### MAIN ####
if __name__ == "__main__":

	# Parse command line arguments
	args = docopt(__doc__)
	fasta_files = args['<FASTA>']
	delimiter = args['--delim']
	fields = args['--fields'].split(',')
	fields = [int(field)-1 for field in fields]  # Convert fields to integers; we subtract 1 to convert to 0-based indexing
	threads = int(args['--threads']) if args['--threads'] else 1
	verbose = args['--verbose']

	# Check that mandatory arguments are provided
	if not delimiter:
		print("Error: Delimiter is required.")
		sys.exit(1)
	if not fields:
		print("Error: Fields are required.")
		sys.exit(1)
	if not fasta_files:
		print("Error: At least one FASTA file is required.")
		sys.exit(1)
	# Check that fields are valid
	for field in fields:
		if field < 0:
			print(f"Error: Invalid field index {field}. Must be a non-negative integer.")
			sys.exit(1)
	# Check that threads is a positive integer
	if threads <= 0:
		print(f"Error: Invalid number of threads {threads}. Must be a positive integer.")
		sys.exit(1)	

	# Process FASTA files in parallel
	with ThreadPoolExecutor(max_workers=threads) as executor:
		futures = [
			executor.submit(subset_fasta_headers, fasta_file, delimiter, fields, threads, verbose)
			for fasta_file in fasta_files
		]
		for future in futures:
			future.result()  # Wait for all to finish
