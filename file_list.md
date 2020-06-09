# YX-AUG Motif Ancestry Data Files

## __Complete Genome Files__
(Uploaded April 24 to Slack by Helena)  

- __picornaviridae_genome_with_outgroup_aligned.fasta__: 228 aligned nucleotide sequences of the complete genome sequence of members of the Picornaviridae family. 
Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. 
- __potyviridae_genome_with_outgroup_aligned.fasta__: 220 aligned nucleotide sequences of the complete genome sequence of members of the Potyviridae family. 
Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __picornaviridae_potyviridae_genome_aligned.fasta__: 446 aligned nucleotide sequences of the complete genome of members of the Potyviridae family and the Picornaviridae family. Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups.

## __5' UTR Files__
(Uploaded April 24 to Slack by Helena)  
- __picornaviridae_5UTR_with_outgroup_aligned.fasta__: 131 aligned nucleotide sequences of the 5' UTR sequence of members of the Picornaviridae family. 
Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. 
- __potyviridae_5UTR_with_outgroup_aligned.fasta__: 154 aligned nucleotide sequences of the 5' UTR sequence of members of the Potyviridae family. 
Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __potyviridae_picornaviridae_5UTR.fasta__: 283 aligned nucleotide sequences of the 5' UTR sequence of members of the Potyviridae family and the Picornaviridae family. Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups.

## __Replicase Files__
(Uploaded May 18 to Slack by Helena)  

- __replicase_potyvirus_nucleotide_aligned.fasta__: 120 aligned nucleotide sequences of the replicase gene of members of the Potyviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_potyvirus_protein_aligned.fasta__: 120 aligned amino acid sequences of the replicase protein of members of the Potyviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_nucleotide_aligned.fasta__: 86 aligned nucleotide sequences of the replicase gene of members of the Picornaviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_protein_aligned.fasta__: 86 aligned amino acid sequences of the replicase protein of members of the Picornaviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_and_potyvirus_nucleotide_aligned.fasta__: 204 aligned nucleotide sequences of the replicase gene of members of the Potyviridae and Picornaviridae families. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_and_potyvirus_protein_aligned.fasta__: 204 aligned amino acid sequences of the replicase protein of members of the Potyviridae and Picornaviridae families. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 


# Script/results files
## Initial batch
Analysis of 4 files in `scripts/virus_project.Rmd`
- picornaviridae_aligned
- potyviridae_aligned
- picornaviridae_5UTR_aligned, and 
- potyviridae_5UTR_aligned

with output files in `results` folder:
- raxml-trees/RAxML_bestTree.T1
- raxml-trees/RAxML_bestTree.T2
- raxml-trees/RAxML_bestTree.T3
- raxml-trees/RAxML_bestTree.T4
- ancestral-sequence/tree3_N0-ancestors_GRASP.fasta

These analysis do not have an outgroup and we use midpoint rooting to root the estimated trees.
