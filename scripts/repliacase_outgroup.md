replicase virus with outgroup
================

1. Construct the tree using RAXML
=================================

Each tree file is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF.

Note: In order to make the output files more clear and display it correctly, I change R1, R2, R3, R4 to .tree in the final file names.

Nucleotide files
----------------

##### picornavirus and potyvirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_and\_potyvirus\_nucleotide\_aligned.fasta -n R1

This command creates the output file: "replicase\_picornavirus\_and\_potyvirus\_nucleotide\_aligned.tree"

##### picornavirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_nucleotide\_aligned.fasta -n R2

This command creates the output file: "replicase\_picornavirus\_nucleotide\_aligned.tree"

##### potyvirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_potyvirus\_nucleotide\_aligned.fasta -n R3

This command creates the output file: "replicase\_potyvirus\_nucleotide\_aligned.tree"

Protein files
-------------

##### picornavirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_protein\_aligned.fasta -n R4

##### potyvirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_potyvirus\_protein\_aligned.fasta -n R5

##### picornavirus and potyvirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_and\_potyvirus\_protein\_aligned.fasta -n R6

### 2. Root the trees at outgroup

##### potyvirus

replicase\_potyvirus\_nucleotide\_aligned.fasta has 120 aligned nucleotide sequences of the replicase gene of members of the Potyviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. replicase\_potyvirus\_protein\_aligned.fasta has 120 aligned amino acid sequences of the replicase protein of members of the Potyviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups.

``` r
library(ape) 
tree_replicase_potyvirus_nucleotide_aligned=read.tree("replicase_potyvirus_nucleotide_aligned.tree")
tree.1=root(tree_replicase_potyvirus_nucleotide_aligned, outgroup=c("NC001747.1:308-17741774-3495Potatoleafrollviruscompletegenome","NC004750.1:142-11551155-2744Barleyyellowdwarfvirus-PAVcompletegenome"))
#write.tree(tree.1,file="replicase_potyvirus_nucleotide_aligned_rooted.tree")

tree_replicase_potyvirus_protein_aligned=read.tree("replicase_potyvirus_protein_aligned.tree")
#tree.2=root(tree_replicase_potyvirus_protein_aligned, outgroup=c("NP_056748.3_RNA-dependent_RNA_polymerase_Potato_leafroll_virus","NP_840014.2_RNA-dependent_RNA_polymerase_P1-P2_fusion_protein_Barley_yellow_dwarf_virus_PAV"))
#write.tree(tree.2,file="replicase_protein_nucleotide_aligned_rooted.tree")
```

##### picornavirus

Next, root the tree of replicase picornavirus nucleotide and protein at outgroup.

replicase\_picornavirus\_nucleotide\_aligned.fasta has 86 aligned nucleotide sequences of the replicase gene of members of the Picornaviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. replicase\_picornavirus\_protein\_aligned.fasta has 86 aligned amino acid sequences of the replicase protein of members of the Picornaviridae family. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups.

``` r
library(ape) 
tree_replicase_picornvirus_nucleotide_aligned=read.tree("replicase_picornavirus_nucleotide_aligned.tree")
tree.3=root(tree_replicase_picornvirus_nucleotide_aligned, outgroup=c("NC004750.1:142-11551155-2744Barleyyellowdwarfvirus-PAVcompletegenome","NC001747.1:308-17741774-3495Potatoleafrollviruscompletegenome"))
#write.tree(tree.3,file="replicase_picornavirus_nucleotide_aligned_rooted.tree")

tree_replicase_picornvirus_protein_aligned=read.tree("replicase_picornavirus_protein_aligned.tree")
tree.4=root(tree_replicase_picornvirus_protein_aligned, outgroup=c("NP_056748.3_RNA-dependent_RNA_polymerase_Potato_leafroll_virus","NP_840014.2_RNA-dependent_RNA_polymerase_P1-P2_fusion_protein_Barley_yellow_dwarf_virus_PAV"))
#write.tree(tree.4,file="replicase_picornavirus_protein_aligned_rooted.tree")
```

##### picornavirus and potyvirus

Last, root the tree of replicase picornavirus and potyvirus nucleotide and protein at outgroup.

replicase\_picornavirus\_and\_potyvirus\_nucleotide\_aligned.fasta has 204 aligned nucleotide sequences of the replicase gene of members of the Potyviridae and Picornaviridae families. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. replicase\_picornavirus\_and\_potyvirus\_protein\_aligned.fasta has 204 aligned amino acid sequences of the replicase protein of members of the Potyviridae and Picornaviridae families. Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups.

``` r
library(ape) 
tree_replicase_picornavirus_and_potyvirus_nucleotide_aligned=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_aligned.tree")
#tree.5=root(tree_replicase_picornavirus_and_potyvirus_nucleotide_aligned, outgroup=c("NC001747.1:308-17741774-3495Potatoleafrollviruscompletegenome","NC004750.1:142-11551155-2744Barleyyellowdwarfvirus-PAVcompletegenome"))
#write.tree(tree.5,file="replicase_picornavirus_nucleotide_aligned_rooted.tree")

tree_replicase_picornavirus_and_potyvirus_protein_aligned=read.tree("replicase_picornavirus_and_potyvirus_protein_aligned.tree")
#tree.6=root(tree_replicase_picornavirus_and_potyvirus_protein_aligned, outgroup=c("NP_840014.2_RNA-dependent_RNA_polymerase_P1-P2_fusion_protein_Barley_yellow_dwarf_virus_PAV"  ,"NP_056748.3_RNA-dependent_RNA_polymerase_Potato_leafroll_virus"))
#write.tree(tree.6,file="replicase_picornavirus_nucleotide_aligned_rooted.tree")
```
