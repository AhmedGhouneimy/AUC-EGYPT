# AUC-EGYPT
This Repo is maintained by Team AUC-EGYPT
AUC-Egypt Software

Background on Toehold switches 
The design process described in the original publications of toehold switches is very lengthy(Green et al., 2014; Pardee et al., 2016). It involves the search for triggers that are structurally predicted to display a minimum number of unpaired bases. Based on these triggers, toehold switches are designed and checked for the presence of stop codons. Candidate toehold switches are assessed in silico based on free energy metrics, structural stability, and orthogonality (See Modeling). Accordingly, multiple iGEM teams sought to automate the design process of toehold switches (Ulaval, CUHK, EPFL). The original class of toehold switches was compatible with prokaryotic and cell free expression systems, but not with mammalian contexts. In our case, we are utilizing toehold switches to sense SARS-CoV-2 mRNA in human host cells which requires some modifications (See Engineering). 

Workflow
To our knowledge, there is no available tool customized to generate mammalian toehold switches. Accordingly, we sought to utilize an available software and customize it to fit our context. Ulaval’s Toeholder was the software of choice as it was developed recently in 2019, and the source code was publicly available. We troubleshooted some errors, enhanced the processing time and processing capacity as described in the Contributions. We also changed the toehold design scheme. Here, we detail the workflow of the improved algorithm. 
The tool identifies candidate triggers of user-specified length, typically around 30 nucleotides. A Candidate trigger is enlisted if it features 2 weak base pairs (A-T) at the stem base of the hairpin(Green et al., 2014). Using a sliding window of width 200 nucleotides, the tool utilizes NUPACK free energy functionalities to predict the secondary structure of the of the 30-nucleotide trigger region, flanked in a 200-nucleotide window. A trigger region is rejected if the number of unpaired bases is lower than that specified by the user. Subsequently, a toehold switch is designed for each candidate trigger as follows:
•	To a triad of GGG, the sensing domain (a & b), reverse complement of the trigger, is appended. The GGG triad improves the transcription efficiency (Green et al., 2014)
•	A rationally engineered loop structure featuring the Kozac consensus sequence is added. 
•	The hairpin is closed with the reverse complement of the ‘b’ domain of the sensing region.
•	A 21-nucleotide linker of low molecular weight amino acids is appended to ensure minimal cross talk between the toehold and downstream gene.
•	If specified by the user, a reporter is added downstream of the toehold switch. 
Candidate toeholds are checked for the presence of a stop codon after the Kozac sequence (the start codon is the last three nucleotides of the Kozac sequence). Using the NUPACK suite, minimal free energy (MFE) structure of the toehold switch is predicted, ΔGtoehold, ΔGRBS-Linker, binding energy, and MFE difference between the bound and the unbound states are calculated and outputted in a single CSV file. After that, toeholds are cross-referenced against user-specified sequences to ensure specificity and minimal crosstalk between the toehold and the host system.
Insert original Workflow design and toehold if possible
 
Improved processing time
We improved the processing time of Toeholder by redefining the trigger regions as described in the Contributions. Briefly, instead of passing the whole input sequence as a trigger region, the sequence is parsed, and only a region of 200 nucleotides spanning the 30-nucleotide trigger is passed to NUPACK minimal free energy (mfe) function. This step was of complexity O(N3), where N is the length of the sequence (Zadeh et al., 2011). Now that we are passing a trigger region of constant length, the complexity becomes O(1). In other words, irrespective of the size of the input sequence the processing time of this step constant. 
To test the processing time of the improved version, we performed a benchmarking experiment. We fed Fasta files of different sizes to both the improved and the original versions of Toeholder.  Our improved version proved it can process at less time (Fig). When the input file size was 1 kb, the new version processed it in 2.5 minutes while the original tool processed ~13 minutes. Moving beyond 1kb, the original tool consumed a full session of 12 hours on Google Colab and did not yield a full output. Comparatively, the same file was processes in ~5 minutes. Moreover, it was able to process a 30 kb file in less than 80 minutes. It is evident from the chart that the processing time of the new version grows linearly rather than exponentially. 
  
Figure 1 Benchmarking 
In the file input_variables.py, the minimum number of paired bases was set to 15, and the length of both the paired and unpaired was set to 15. We passed an empty file for the cross referencing test as it was not necessary for our benchmarking. It is clear that our improved version processes the file size in less amount of time.
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
والسلام عليكم ورحمة الله وبركاته

References 
Green, A. A., Silver, P. A., Collins, J. J., & Yin, P. (2014). Toehold Switches: De-Novo-Designed Regulators of Gene Expression. Cell, 159(4), 925–939. https://doi.org/10.1016/j.cell.2014.10.002
Pardee, K., Green, A. A., Takahashi, M. K., Braff, D., Lambert, G., Lee, J. W., Ferrante, T., Ma, D., Donghia, N., Fan, M., Daringer, N. M., Bosch, I., Dudley, D. M., O’Connor, D. H., Gehrke, L., & Collins, J. J. (2016). Rapid, Low-Cost Detection of Zika Virus Using Programmable Biomolecular Components. Cell, 165(5), 1255–1266. https://doi.org/10.1016/j.cell.2016.04.059
Zadeh, J. N., Wolfe, B. R., & Pierce, N. A. (2011). Nucleic acid sequence design via efficient ensemble defect optimization. Journal of Computational Chemistry, 32(3), 439–452. https://doi.org/10.1002/jcc.21633


