#!/usr/bin/env bash

# Convert GenBank format to EMBL format using Biopython.
# Usage: gbk2embl.sh input.gbk output.embl

infile="$1"
outfile="$2"
python -c "import sys; from Bio import SeqIO; SeqIO.convert(sys.stdin, 'genbank', sys.stdout, 'embl');" < $infile > $outfile
