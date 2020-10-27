# AUC-EGYPT
This Repo is maintained by Team AUC-EGYPT
# Toeholder 2.0

![](Figures/toeholder.png)

## Description

Toeholder is a tool that can efficiently design toehold riboswitches for the detection of a target gene in mammalian systems. A toehold riboswitch is an RNA molecule that contains the necessary elements for the expression of a reporter gene, i. e. kozak sequence (a ribosome binding sequence (RBS) and a start codon), and the reporter gene sequence. What makes toehold riboswitches special is that they fold into a secondary structure that blocks the access of the ribosome to the Kozak sequence, thus preventing gene expression. However, they are carefully designed so that they can bind to a trigger sequence within a target gene, which results in the unfolding of the secondary structure and allows the expression of the reporter gene (Green et al. 2014. Cell.). Applications for toehold riboswitches range from detection of sequences of pathogenic organisms (Pardee, et al. 2016. Cell; Ma, et al. 2018. Synthetic Biology) to creating logical systems (Green et al. 2014. Cell; Green, et al. 2017. Nature).

![](Figures/toehold_diagram.png)

Toeholder aims to facilitate the design of these molecules by offering the following capabilities:
- Testing all the different trigger sequences within a target gene.
- Simulating the secondary structure of the proposed riboswitch for each trigger.
- Simulating the binding of each toehold riboswitch to the target gene to test if it binds accurately to its corresponding trigger.
- Aligning to user-defined genomes in order to select generalist or specific riboswitches as needed

The above tests have allowed us to design toehold riboswitches for different organisms. Secondary structures obtained are very close to the ones observed for the riboswitches published by Green et al. (2014), and 70-79% of them are predicted to bind perfectly to their trigger sequences within the target gene. The remaining 21-30% are identified as binding partially or with a shift to the target or as having stop codons at undesirable positions, which allows discarding them. Finally, the alignment function has allowed us to design toeholds that are specific for the genome of a given strain, but also to select those that have matches in several strains of the same organism. All in all, these functions make toeholder a very versatile tool.

For further information, please refer to our wiki:
https://2020.igem.org/Team:AUC-EGYPT

## Dependencies and Installation

Installation in 6 simple steps!

To start working with the tool, users are opted to follow these 6 simple steps to install the compatible dependencies.

•	Download NUPACK 3.2.2 (http://nupack.org/downloads) and install as follows:

Go to the root directory type the following commands: 

mkdir build
cd build
cmake ../
make
make install
•	Install blast+: sudo apt-get install ncbi-blast+
•	Install Biopython 1.73: pip3 install biopython==1.73
•	Install Numpy 1.16.4: pip3 install numpy==1.16.4
•	Install pandas 0.24.2: pip3 install pandas==0.24.2
•	Download Toeholder (Git link) and go to the root directory and type 
python3 toeholder.py 

How to use the tool

To run the tool, users are required to fill in the inputs in the input_variables.py file. They are asked to provide the following:
•	Path to the input Fasta file
•	Path to the output folder
•	Length of the unpaired domain (the ‘a’ domain) 
•	Length of the paired domain (the ‘b’ domain) 
•	A path to reference genome for cross referencing if any and the percentage identity and the e-value.
•	The molecule type: “DNA” or “RNA” 
•	Reporter gene if any
•	Minimum number of unpaired bases in the trigger. 
After installing all the libraries and dependencies, users are opted to open the terminal in the root directory of Toeholder and type ‘python3 toeholder.py’. Results are incremented to the output folder path specified by the user.

## Scripts

All scripts are written in Python 3 and depend on the following libraries

- toeholder.py: Sweeps through the sequence of the target gene looking for suitable candidate recognition sequences. All candidate recognition sequences are evaluated based on the following parameters:
	- Secondary structure on the mRNA
	- ddG of the bound (toehold + target) and unbound state (toehold and target, separately)

- input_variables.py: Defines tunable parameters and input files for the toeholder script.

- toeholder_helper_functions.py: Contains several helper functions for the other scripts.

- alignments.py: Aligns the toeholds to a selected set of reference genomes to identify matches. Toeholds matching more than one sequence or matching sequences from other genomes would not be completely specific to the target.


## Output

The toeholder.py script generates an output folder with a subfolder for each of the candidate toeholds generated. When the candidate toehold contains a stop codon, its corresponding subfolder is empty. When it does not contain a stop codon, there are four files inside the subfolder:
- switch1_python.in: NUPACK-formatted input to test the toehold's secondary structure.
- switch1_python.mfe: NUPACK-formatted output with the most favorable secondary structure.
- toehold_mRNA.in: NUPACK-formatted input to test the toehold's ability to bind to the target mRNA.
- toehold_mRNA.mfe: NUPACK-formatted output with the most favorable structure of the toehold-mRNA complex.
- RBS_LINKER.in:NUPACK-formatted input to test the RBS-LINKER secondary structure.
- RBS_LINKER.mfe: NUPACK-formatted output with the most favorable secondary structure of the RBS-LINKER.

Outside those folders, the rest of the files are:
- input_variables.py: Copy of the input variables used for this run.
- mRNA.in: NUPACK-formatted input to test the mRNA's secondary structure. This file is not correctly formated in Toeholder.py as it contains information about the last trigger region only, not the full input mRNA
- mRNA.mfe: NUPACK-formatted output with the most favorable secondary structure for the mRNA. This file is not correctly formated in Toeholder.py as it contains information about the last trigger region only, not the full input mRNA
- toehold_candidates.txt: List of candidate triggers ranked by the number of non-paired positions in the secondary structure of the trigger of the target mRNA.
- toehold_seqs.fasta: FASTA-formatted file containing the recognition sequences for each of the toeholds.
- all_toeholds_results.txt: Results of the tests performed on the toeholds. It adds the following columns to the toehold_candidates.txt file:
	- Toehold index
	- Binding_energy_toehold_mRNA: the free energy of binding between the toehold and the 200-nt Trigger region 
	- Percentage of paired bases of the toehold-mRNA complex that correspond the intended base pairing
	- Binding_energy_toehold:Free energy of the toehold secondary structure
	- GC content: in Toehold 
	- deltaG_RBS_LINKER: free energy of the RBS-Linker
	- MFE_Difference: free energy difference between the bound and the unbound state
	- toehold structure: Minimal free energy structure in the .( format
	- toehold sequence
	
	
- \<tag\>_toeholds_alignment.aln: Output of the BLAST alignment of the library of toeholds to the genome referenced with the corresponding tag in the genome list.
- all_toeholds_results_genome_matches.txt: Adds the counts of matches for each toehold in each of the genomes referenced in the genome list.




