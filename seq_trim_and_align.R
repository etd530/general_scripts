#!/usr/bin/env Rscript

#### LIBS ####
library(Biostrings)
library(msa)
library(DECIPHER)

#### VARS ####
fasta_file <- "DataBarcodingLibrary/fasta_genera_not_aligned_ONLYNEW/Formica_not_aligned.fasta"  # Replace with your file path
sequences <- readDNAStringSet(fasta_file)


#### MAIN ####
# Create a new sequence
new_sequence <- DNAStringSet("TATCTTATACTTCCTTCTTGCTATCTGAGCTAGTCTAGTAGGATCTTCAATAAGAATAATTATCCGATTAGAACTAGGTACTTGCAACTCATTAATTAATAATGATCAAATTTATAATACCGTAATTACTAGCCACGCTTTTATTATAATTTTCTTTATAGTTATACCCTTTATAATTGGTGGGTTTGGTAATTTTCTTGTCCCCTTGATACTAGGCTCCCCTGATATAGCCTACCCTCATTTAAATAATATAAGTTTCTGATTATTGCCCCCCTCAATGGCCCTACTTTTATTAAGCAATTTTATTAATGTCGGAGTAGGAACTGGATGAACTATTTACCCTCCCCTCGCCTCTAATATCTTTTACTCCGGCCCCTCAATTGATTTATCAATTTTTTCTCTTCATATCGCCGGAATATCATCAATTCTTGGCGCTATTAATTTTGTTTCTACAATTTTAAACATACATCATAAAAATTTTTCTATAGAAAAAATTCCATTACTAGTTTGATCAATTATAATTACAGCAGTTCTACTACTACTTTCCTTATCAGTATTAGCAGGAGCTATTACTATGCTTTTAACTGATCGTAACTTAAATACATCATTTTTTGACCCATCTGGGGGAGGAGACCCTATTTTATACCAACATTTATTT")  # Replace with your sequence
names(new_sequence) <- "JQ742645"

# Combine the new sequence with the existing ones
combined_sequences <- append(sequences, new_sequence)

# Perform multiple sequence alignment
alignment <- msa(combined_sequences, method = "ClustalW")

# Convert the alignment to a DNAStringSet object for further manipulation
aligned_sequences <- as(alignment, "DNAStringSet")
BrowseSeqs(aligned_sequences)

# Step 3: Get the reference sequence
reference_sequence <- aligned_sequences[length(aligned_sequences)]  # Specify the reference sequence
library(stringr)
reference_sequence <- str_split(as.character(reference_sequence), "")[[1]]

# Step 4: Trim all sequences to the start and end of the reference sequence
trimmed_sequences <- subseq(aligned_sequences, start = min(which(reference_sequence != "-")),
                                               end = max(which(reference_sequence != "-")))
BrowseSeqs(trimmed_sequences)

# Remove reference by name
trimmed_sequences <- trimmed_sequences[!names(trimmed_sequences) %in% "sequence_to_remove"]

# Step 6: Save the trimmed sequences to a new FASTA file
output_fasta <- "DataBarcodingLibrary/fasta_genera_not_aligned_ONLYNEW/Formica_not_aligned_OK.fasta"
writeXStringSet(trimmed_sequences, filepath = output_fasta)

cat("Trimmed sequences saved to:", output_fasta)
