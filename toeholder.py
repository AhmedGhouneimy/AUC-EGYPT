#!/usr/bin/env python
# coding: utf-8

# # Designing toeholds from a target sequence or genome
# 
# This script will receive a target sequence or genome and look for candidates of sites for which we could design toeholds. Candidates will be selected as follows:
# - Regions having 2 weak pairs and 1 strong at the base of the candidate
# - User defines lengths a (unpaired part of the recognition sequence) and b (paired part)


# Load libraries
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from Bio import SeqIO
import os
import csv
from collections import OrderedDict
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shutil import copyfile
import glob
from toeholder_helper_functions import *
from input_variables import *

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Make a copy of the input parameter file in the new directory
copyfile('input_variables.py', os.path.join(output_folder, 'input_variables.py'))

selected = []

# ## 1.- Use a loop to find candidates
# Parse the input sequence
input_seq_record = next(SeqIO.parse(input_seq, 'fasta'))

# If input sequence is DNA, parse and transcribe
if mol_type == 'DNA':   
    full_input_seq = str(input_seq_record.seq).upper().replace('T', 'U')   
# If it is RNA, I should prepare a NUPACK run to look at its secondary structure
elif mol_type == 'RNA':
    full_input_seq = str(input_seq_record.seq).upper()
#print(full_input_seq)    
# Calculate the secondary structure with NUPACK

output_file = os.path.join(output_folder, 'mRNA.in')#output file is the path to the mrna input to nupack
"""
prepare_nupack_input(full_input_seq, output_file)
nupack_mfe(output_file[0:-3])#Now the output file conatins the nupack outout

# Parse with NUPACK
mRNA_dict = parse_nupack(os.path.join(output_folder, 'mRNA.mfe'))
mRNA_structure = mRNA_dict['structure']
mRNA_sequence = mRNA_dict['sequence']
"""
#AUC-Egypt Was Here
def get_trigger_region(pos, seq, len_trigger, len_sensed):
    #critical_pt = (len_trigger - len_sensed)//2
    critical_pt = len_trigger //2 #half point in the trigger
    trigger = ""
    if pos >= critical_pt and pos <= len(seq)-critical_pt: #if it's far away from both ends --> trigger seq is in the middle
        start = pos-critical_pt#- int(0.5*len_sensed)
        end = pos+critical_pt#+ int(0.5*len_sensed)
        trigger = seq[start :  end]
        sp = critical_pt-int(0.5*len_sensed) #sp start pos in the short trigger seq
        ep =critical_pt+int(0.5*len_sensed) #ep end pos in the short trigger seq
        #print (ep == sp+len_sensed)
    elif pos <= critical_pt:
        start = 0
        end = len_trigger
        trigger = seq[start:end]
        sp = pos
        ep = pos + len_sensed
    elif pos >= len(seq) - critical_pt:
        start = len(seq)-len_trigger
        end = len(seq)-1
        trigger = seq[start:]
        sp = len_trigger-(pos % critical_pt)
        ep = sp + len_sensed
    return trigger, sp, ep
        
for hairpin_start_pos in range(length_unpaired, len(full_input_seq) - length_paired):
    # Check the base of the hairpin
    #hairpin_base = mRNA_sequence[hairpin_start_pos:hairpin_start_pos + 3]
    hairpin_base = full_input_seq[hairpin_start_pos:hairpin_start_pos + 3]

    # Discard cases for which the base of the hairpin does not have two weak pairs and a strong one
    if hairpin_base.count('A') + hairpin_base.count('U') != 2:
        continue
    pos = hairpin_start_pos - length_unpaired
    seq = full_input_seq
    len_trigger = 200 #make it user defined
    len_sensed = 30 #make it user defined
    trigger, start_pos, end_pos = get_trigger_region(pos, seq, len_trigger, len_sensed)
    #end_pos = hairpin_start_pos + length_paired
    #start_pos = hairpin_start_pos - length_unpaired
    #sub_seq = mRNA_sequence[start_pos:end_pos]
    #NUPACK THE TRIGGER
    prepare_nupack_input(trigger, output_file)
    nupack_mfe(output_file[0:-3])#Now the output file conatins the nupack outout

# Parse with NUPACK
    mRNA_dict = parse_nupack(os.path.join(output_folder, 'mRNA.mfe'))
    mRNA_structure = mRNA_dict['structure']
    mRNA_sequence = mRNA_dict['sequence']
    
    mRNA_binding_energy = mRNA_dict['energy']

    sub_seq = mRNA_sequence[start_pos:end_pos]

    sub_struc = mRNA_structure[start_pos:end_pos]
    
    # Count the number of unpaired bases
    unpaired = sub_struc.count('.')
    if unpaired >= min_unpaired:
        selected.append([unpaired, sub_struc, sub_seq, start_pos + 1, end_pos, length_unpaired, length_paired, mRNA_sequence, mRNA_structure, mRNA_binding_energy])
    
tmp = pd.DataFrame(np.array(selected), columns=['Non paired count', 'Structure', 'Sequence', 'Start position in trigger ', 'End position in trigger', 'Length unpaired trigger', 'Length paired trigger','Full trigger seq','full trigger struct','full_trigger_binding_energy'])
df = tmp.sort_values(by = 'Non paired count', ascending = False)

df.to_csv(path_or_buf=os.path.join(output_folder, 'toehold_candidates.txt'), sep = '\t', index = False)
header_df = list(df.columns)

# Loop through the folders and gather the information in a table
results = []

# Retrieve the binding energy of the mRNA on its own
#mRNA_binding_energy = mRNA_dict['energy']

# Prepare toeholds for each of the sequences
for index, row in df.iterrows():
    sequence = row[2] #sensed region
    full_trigger= row[7]#full trigger
    trigger_energy = row[9]
    

    # ## 2.- Generate toeholds for each of the candidates
    toehold_folder = os.path.join(output_folder, str(index))#1,2,3,4,....
    generate_toehold(sequence, mol_type, reporter, toehold_folder, length_unpaired) #sequence is the portion sensed

    # ## 3.- Check how well the toeholds bind to the target in the mRNA
    if not os.path.exists(os.path.join(toehold_folder, 'switch1_python.mfe')):
        # This is the case in which my generator function found a stop codon and stopped
        continue
    # For all other cases, I need to write down the file that will work as the input for NUPACK
    # Read the file to extract the toehold's sequence
    toehold_dict = parse_nupack(os.path.join(toehold_folder, 'switch1_python.mfe'))
    toehold = toehold_dict['sequence']
    toehold_struct = toehold_dict['structure']
    toehold_binding_energy = toehold_dict['energy']
    
    output_file = os.path.join(toehold_folder, 'toehold_mRNA.in')
    prepare_nupack_input_two(toehold, full_trigger, output_file) #potentially pass a small trigger
    
    # Run NUPACK
    nupack_mfe(output_file[0:-3])

    # Use the last toehold we prepared to save the length of the ones I am generating
    len_toehold = len(toehold)

    # Retrieve the information from the table about this toehold
    # toehold_num = int(toehold_folder.split('/')[-2])
    toehold_num = int(index)
    toehold_data = [entry for entry in row]#??
    #print ("toehold_data")
    #print (toehold_data) -->#['16', '(.((....(((.((.((((((.........', 'GUUGAUUUUUGUGGAAAGGGCUAUCAUCUU', '86', '115', '15', '15', 'AUUAAUUAGAGCUGCAGAAAUCAGAGCUUCUGCUAAUCUUGCUGCUACUAAAAUGUCAGAGUGUGUACUUGGACAAUCAAAAAGAGUUGAUUUUUGUGGAAAGGGCUAUCAUCUUAUGUCCUUCCCUCAGUCAGCACCUCAUGGUGUAGUCUUCUUGCAUGUGACUUAUGUCCCUGCACAAGAAAAGAACUUCACAACUG', '.....(((((((.((((....((((...))))......)))).))).))))..((((.((((....)))).))))........(((.((....(((.((.((((((...........)))))).)).)))....)).)))...(((.((((((((((((.(.(((....)))).))))..)))..)).))).))).....', '-45.3']
    #toehold_data is the trigger info and it will append other values
    # toehold_data = [entry for entry in df.iloc[toehold_num]]

    # Check if we have a result
    result_path = os.path.join(toehold_folder, 'toehold_mRNA.mfe')

    # start_pos = int(toehold_data[3]) + len_toehold
    # end_pos = int(toehold_data[4]) + len_toehold
    start_pos = int(toehold_data[3]) + len_toehold #start pos in full trigger
    end_pos = int(toehold_data[4]) + len_toehold # end_pos in full trigger

    toehold_data.append(toehold_num)
    #   RBSLINKER
    RBS_LINKER_dict = parse_nupack(os.path.join(toehold_folder, 'RBS_LINKER.mfe'))
    deltaG_RBS_LINKER = RBS_LINKER_dict['energy']
    
    if os.path.exists(result_path): # this if appends switch-mrna, percent maches and siwtch energy
        # Parse the output file
        toehold_mRNA_dict = parse_nupack(result_path)
        
        # Add the binding energy mrna+switch
        binding_energy = toehold_mRNA_dict['energy']
        toehold_data.append(binding_energy)
        #binding_energy = trigger_energy #trigger energy is already there
        #toehold_data.append(binding_energy)
        
        # Add to the current row the percentage of residues that are successfully bound
        pair_list = toehold_mRNA_dict['paired_bases']
        percentage_correct = check_matches(start_pos, end_pos, pair_list) #check the start and end and print pair_list
        toehold_data.append(percentage_correct)
        
        # Get the binding energy for the toehold on its own
        '''
        toehold_alone_path = os.path.join(toehold_folder, 'switch1_python.mfe')
        toehold_alone_dict = parse_nupack(toehold_alone_path)
        toehold_alone_binding_energy = toehold_alone_dict['energy']
        toehold_data.append(toehold_alone_binding_energy)
        #toehold_data.append(deltaG_RBS_LINKER)
        '''
        toehold_data.append(toehold_binding_energy)
    else:
        # There are no values for binding energy or percentage of correct bases
        toehold_data.append('NA')
        toehold_data.append('NA')
        toehold_data.append('NA')
    
    # Add the binding energy of the mRNA on its own
    #toehold_data.append(trigger_energy) already appended
    
    # Look at the GC content
    sequence = toehold_data[2]
    gc_content = round(float(sequence.count('G') + sequence.count('C'))*100/len(sequence), 2)
    toehold_data.append(gc_content)
    toehold_data.append(deltaG_RBS_LINKER)
    #deltaG_binging = binding energy - (toehold_alone + rna_alone)

    MFE_Difference = float(binding_energy) - (float(toehold_binding_energy) + float(trigger_energy))
    toehold_data.append(MFE_Difference)
    toehold_data.append(toehold_struct)
    toehold_data.append(toehold)

    # Add this row to the list of results
    results.append(toehold_data)
    
# Save as a dataframe #### Remove references no non-paired positions
results_df = pd.DataFrame(np.array(results), columns=['Non paired count', 'Structure', 'Sequence', 'Start position in trigger ', 'End position in trigger', 'Length unpaired trigger', 'Length paired trigger','Full trigger seq','full trigger struct','full_trigger_binding_energy', 'Index', 'Binding_energy_toehold_mRNA', 'Percentage_correct_matches', 'Binding_energy_toehold', 'GC content', 'deltaG_RBS_LINKER', 'MFE_Difference','toehold structure', 'toehold sequence'])

sorted_results = results_df.sort_values(['Binding_energy_toehold_mRNA', 'Percentage_correct_matches'], ascending = [False, True])
sorted_results.to_csv(path_or_buf=os.path.join(output_folder, 'all_toeholds_results.txt'), sep = '\t', index = False)

#### Add call to the script that aligns to the genomes ####
part1 = 'python alignments.py'
part2 = ' -w ' + output_folder
part3 = ' -r ' + reference_list
part4 = ' -p ' + str(pct_ident)
part5 = ' -e ' + str(evalue)
os.system(part1 + part2 + part3 + part4 + part5)
