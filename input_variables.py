#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script saves some input variables needed by the toehold_master.py script

# The path to the sequence for which toeholds will be designed (FASTA-formatted)
input_seq = "<>"
# The length of the unpaired region of the recognition sequence (upstream of the hairpin) 
length_unpaired = <>

# The length of the paired region of the recognition sequence (in the hairpin) --> has to be a multiple of 3 (15 for example)
length_paired = <>

# The path to the output folder
output_folder = "<>"

# The reporter sequence to be simulated at the end of each toehold switch
reporter = ""

# The molecule type that is provided as input (RNA or DNA)
mol_type = ""

# The path to the list of reference genomes to which candidate toeholds should be mapped
# The list should be tab-delimited, with two columns:
# 1.- Path to the genome file
# 2.- Name tag for this genome
reference_list = "<>"

# The lower threshold of sequence identity for the hits to be retained
pct_ident = 90

# The upper threshold for the evalue of hits to be retained
evalue = 1e-6

# The minimum number of unpaired residues in the secondary structure of the target RNA for a candidate trigger to be considered
min_unpaired = <>
