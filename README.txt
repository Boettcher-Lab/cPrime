Wir geben alle Rechte zur Nutzung frei




Count_edits_in_fastq: 

Contains bash and python script to count edits in fastq files as well as reference files to search edits. 

Count_edits.sh -> you have to change in line 119 the reference file depending on if you want to search for hits in:
 - gDNA (with synonymous marker)
 - cDNA (with synonymous marker)
 - gDNA (with wild type sequence (WT); no intended or synonymous marker edit)
 - gDNA (without synonymous marker; only intended edit)


Line 119:  awk '{ print $2 }' Searching_sequences_gDNA.txt > sequences.txt


-----------------------------------------------------------------------------------

Design_pegRNAs_additional_STOP_codon: 

Contains R script to design PrimeDesign input files for nonsense mutations and synonymous mutations on codon positions in exon 27 that are not described in cBioPortal. Also includes R scripts to generate the reference files with the search sequences, to determine the number of each intended edit in the fastq files.

- Design_input_sequences_for_PrimeDesign_STOP_and_synonymous_edit.R -> design input sequences for PrimeDesign
- Design_searching_sequences_for_STOP_and_synonymous_edit_cDNA.R -> design searching sequences for additional nonsense (STOP) and synonymous mutations for cDNA
- Design_searching_sequences_for_STOP_and_synonymous_edit_gDNA.R -> design searching sequences for additional nonsense (STOP) and synonymous mutations for gDNA
- Design_WT_sequences_for_STOP_and_synonymous_edit_cDNA.R -> design WT sequences for additional nonsense (STOP) and synonymous mutations for cDNA
- Design_WT_sequences_for_STOP_and_synonymous_edit_gDNA.R -> design WT sequences for additional nonsense (STOP) and synonymous mutations for gDNA

- Chr17_full_sequence.txt -> Sequence of chromosome 17 from hg19 genome
- NF1-Exons-position.txt -> Information about start and end position of each exon from NF1 in chromosome 17

Additional nonsense mutations or synonymous edits that are missing in the final Overview_names.xlsx list are removed by hand for example all InDels are removed.

-----------------------------------------------------------------------------------

Design_pegRNAs_cBioPortal_mutations:

Contains R script to design PrimeDesign input files for the mutations listed in cBioPortal for exon 27. Also includes R scripts to generate the reference files with the search sequences, to determine the number of each intended edit in the fastq files.

- cBioPortal_27Jan2022.tsv -> cBioPortal mutations for NF1 (Download: 27 January 2022)
- Chr17_full_sequence.txt -> Sequence of chromosome 17 from hg19 genome
- NF1-Exons-position.txt -> Information about start and end position of each exon from NF1 in chromosome 17
- Exon_splice_information.txt -> Information about number of bases of the first and last codon of each exon. 
For example: Start_codon_bp_count 3 means that the first codon of the exon starts with three bases; Start_codon_bp_count 2 means that the exon starts with two bases while the "third" base for this codon is located on the exon before (her the last base).
- Design_input_sequences_mutations_for_PrimeDesign.R -> design input sequences for PrimeDesign for mutations
- Design_input_sequences_synonymous_edit_for_PrimeDesign.R -> design input sequences for PrimeDesign for synonymous edit
- Design_searching_sequences_for mutations_cDNA.R -> design searching sequences for mutations for cDNA
- Design_searching_sequences_for mutations_gDNA.R -> design searching sequences for mutations for gDNA
- Design_searching_sequences_for_synonymous_edit_cDNA.R -> design searching sequences for synonymous edit for cDNA
- Design_searching_sequences_for_synonymous_edit_gDNA.R -> design searching sequences for synonymous edit for gDNA
- Design_WT_searching_sequences_for_mutations_cDNA.R -> design WT sequences for mutations for cDNA
- Design_WT_searching_sequences_for_mutations_gDNA.R -> design WT sequences for mutations for gDNA
- Design_WT_searching_sequences_for_synonymous_edit_cDNA.R -> design WT sequences for synonymous edit for cDNA
- Design_WT_searching_sequences_for_synonymous_edit_gDNA.R -> design WT sequences for synonymous edit for gDNA


Mutations or synonymous edits that are missing in the final Overview_names.xlsx list are removed by hand for example all InDels are removed.

-----------------------------------------------------------------------------------

Odds_Ratio:

Contains R script and reference files to calculate the Odds Ratio.
