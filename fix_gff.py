#!/usr/bin/env python3

"""
Correct GFF Coordinates
Usage:
  correct_gff_coordinates.py <fasta_file> <gff_file> <output_file>
  correct_gff_coordinates.py (-h | --help)

Arguments:
  <fasta_file>   Path to the genome FASTA file.
  <gff_file>     Path to the GFF file with genome-wide coordinates.
  <output_file>  Path to save the corrected GFF file.

Options:
  -h --help      Show this help message.
"""

#### LIBS ####
from Bio import SeqIO
from docopt import docopt

#### FUNS ####
def parse_fasta_sizes(fasta_file):
    """
    Parse the FASTA file and calculate the sizes of each contig.
    
    :param fasta_file: Path to the FASTA file.
    :return: A dictionary mapping contig names to their sizes.
    """
    contig_sizes = {}
    with open(fasta_file, 'r') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            contig_sizes[record.id] = len(record.seq)
    return contig_sizes

def adjust_gff_coordinates(gff_file, contig_sizes, output_file):
    """
    Adjust GFF coordinates based on contig sizes.
    
    :param gff_file: Path to the GFF file with genome-wide coordinates.
    :param contig_sizes: Dictionary of contig sizes.
    :param output_file: Path to save the corrected GFF file.
    """
    with open(gff_file, 'r') as gff, open(output_file, 'w') as out:
        current_contig = None
        current_offset = 0
        contig_boundaries = {}
        
        # Compute the cumulative offsets for contigs
        for contig, size in contig_sizes.items():
            contig_boundaries[contig] = current_offset
            current_offset += size
        
        for line in gff:
            # Skip comments and header lines
            if line.startswith('#'):
                out.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                out.write(line)
                continue  # Malformed line, just write it unchanged
            
            start = int(fields[3])
            end = int(fields[4])
            
            # Identify the corresponding contig
            for contig, offset in contig_boundaries.items():
                if start > offset and start <= offset + contig_sizes[contig]:
                    current_contig = contig
                    start -= offset
                    end -= offset
                    break
            
            # Adjust the fields
            fields[0] = current_contig
            fields[3] = str(start)
            fields[4] = str(end)
            
            # Write the updated line to the output file
            out.write('\t'.join(fields) + '\n')

#### MAIN ####
if __name__ == "__main__":
    # Parse command-line arguments
    args = docopt(__doc__)
    
    # Input files
    fasta_file = args['<fasta_file>']
    gff_file = args['<gff_file>']
    output_file = args['<output_file>']
    
    # Parse contig sizes from the FASTA file
    contig_sizes = parse_fasta_sizes(fasta_file)
    
    # Adjust GFF coordinates
    adjust_gff_coordinates(gff_file, contig_sizes, output_file)
    
    print(f"Corrected GFF saved to {output_file}")
