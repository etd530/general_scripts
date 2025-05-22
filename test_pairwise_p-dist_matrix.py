#!/usr/bin/env python3

#### LIBS ####
import pairwise_pdist_matrix as p_dist

#### Tests for number_of_substitutions() ####
def test_number_of_substitutions():
	assert p_dist.number_of_substitutions("ACGT", "ACGT") == 0
	assert p_dist.number_of_substitutions("ACGT", "ACGA") == 1
	assert p_dist.number_of_substitutions("ACGT", "TGCA") == 4
	assert p_dist.number_of_substitutions("ACGT", "ACG-") == 1
	assert p_dist.number_of_substitutions("ACGT", "ACG-", no_gaps = True) == 0
	assert p_dist.number_of_substitutions("ACGT--", "--GTAA") == 4
	assert p_dist.number_of_substitutions("ACGT--", "--GTAA", no_gaps = True) == 0
	assert p_dist.number_of_substitutions("A-G-", "A-G-") == 0
	assert p_dist.number_of_substitutions("A-G-", "A-G-", no_gaps = True) == 0
def test_ali_len_msa2pairwise():
	assert p_dist.ali_len_msa2pairwise("ACGT", "ACGT") == 4
	assert p_dist.ali_len_msa2pairwise("ACGT", "ACGA") == 4
	assert p_dist.ali_len_msa2pairwise("ACGT", "TGCA") == 4
	assert p_dist.ali_len_msa2pairwise("ACGT", "ACG-") == 4
	assert p_dist.ali_len_msa2pairwise("ACGT", "ACG-", no_gaps = True) == 3
	assert p_dist.ali_len_msa2pairwise("ACGT--", "--GTAA") == 6
	assert p_dist.ali_len_msa2pairwise("ACGT--", "--GTAA", no_gaps = True) == 2
	assert p_dist.ali_len_msa2pairwise("A-G-", "A-G-") == 2
	assert p_dist.ali_len_msa2pairwise("A-G-", "A-G-", no_gaps = True) == 2
