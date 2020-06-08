# YX-AUG Motif Ancestry Data Files

## __Replicase Files__
(Uploaded May 18 by Helena)

- __replicase_potyvirus_nucleotide_aligned.fasta__: 120 aligned nucleotide sequences of the replicase gene of members of the Potyviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_potyvirus_protein_aligned.fasta__: 120 aligned amino acid sequences of the replicase protein of members of the Potyviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_nucleotide_aligned.fasta__: 86 aligned nucleotide sequences of the replicase gene of members of the Picornaviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_protein_aligned.fasta__: 86 aligned amino acid sequences of the replicase protein of members of the Picornaviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_and_potyvirus_nucleotide_aligned.fasta__: 204 aligned nucleotide sequences of the replicase gene of members of the Potyviridae and Picornaviridae families. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 
- __replicase_picornavirus_and_potyvirus_protein_aligned.fasta__: 204 aligned amino acid sequences of the replicase protein of members of the Potyviridae and Picornaviridae families. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. 


# Script/results files
1. Initial batch: analysis of 4 files in `scripts/virus_project.Rmd`
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
