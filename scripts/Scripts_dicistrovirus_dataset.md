Scripts\_dicistrovirus\_dataset
================

New Dataset (Supplemental data dicistrovirus)
---------------------------------------------

First part: dicistrovirus I

<file:%22dicistrovirus_I_w_outgroup_aligned.fasta>": 7 aligned nucleotide sequences of the 5' UTR of members of Dicistroviridae family with a characterized IRES plus 1 outgroup (8 sequence total: 4 Cripavirus, 3 Aparavirus, 1 outgroup). Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s dicistrovirus\_I\_w\_outgroup\_aligned.fasta -n T1

This command creates the output file: "dicistrovirus\_I\_w\_outgroup\_aligned\_woroot.tree"

``` r
library(ape)
library(phylotools)
data=read.fasta("dicistrovirus_I_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("dicistrovirus_I_w_outgroup_aligned_woroot.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")
index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"dicistrovirus_I_wo_outgroup_aligned_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"dicistrovirus_I_wo_outgroup_aligned.fasta")
```

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t dicistrovirus\_I\_wo\_outgroup\_aligned\_w\_root.tree -s dicistrovirus\_I\_wo\_outgroup\_aligned.fasta -m GTRGAMMA -n S1

This command product output file: "dicistrovirus\_I\_wo\_outgroup\_AncestralSeq.fasta"

Joint ancestors file: "dicistrovirus\_I\_wo\_outgroup\_aligned\_joint\_ancestors.fasta" Tree with label <file:%22dicistrovirus_I_wo_outgroup_aligned_w_root_w_label.tree>"

Second part: dicistrovirus II <file:%22dicistrovirus_II_w_outgroup_aligned.fasta>"

7 aligned nucleotide sequences of the intergenic region of members of Dicistroviridae family with a characterized IRES plus 1 outgroup (8 sequence total: 4 Cripavirus, 3 Aparavirus, 1 outgroup). Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s dicistrovirus\_II\_w\_outgroup\_aligned.fasta -n T2

This command creates the output file: "dicistrovirus\_II\_w\_outgroup\_aligned\_woroot.tree"

``` r
rm(list=ls())
library(ape)
library(phylotools)
data=read.fasta("dicistrovirus_II_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("dicistrovirus_II_w_outgroup_aligned_woroot.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")


index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"dicistrovirus_II_wo_outgroup_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"dicistrovirus_II_wo_outgroup_aligned.fasta")
```

Note: the gaps here are removed

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t dicistrovirus\_II\_wo\_outgroup\_w\_root.tree -s dicistrovirus\_II\_wo\_outgroup\_aligned.fasta -m GTRGAMMA -n S2

This command product output file: "dicistrovirus\_II\_wo\_outgroup\_aligned\_AncestralSeq.fasta"

Joint ancestors file: "dicistrovirus\_II\_wo\_outgroup\_aligned\_joint\_ancestors.fasta" Tree with label <file:%22dicistrovirus_II_wo_outgroup_aligned_w_root_w_label.tree>"

##### Constructing the ancestral sequences and their confidence value at every internal nodes

Note: The method I use could only generatet the confidence value at one internal nodes. I will combine the confidence value into a single file. (e.g, cat *ConfValue* &gt;&gt; combined.txt) (N0,N1,...,)

1.dicistrovirus I

java -jar bnkit.jar -aln dicistrovirus\_I\_wo\_outgroup\_aligned.fasta -nwk dicistrovirus\_I\_wo\_outgroup\_w\_root.tree -out dicistrovirus\_I\_wo\_outgroup\_aligned\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "dicistrovirus\_I\_wo\_outgroup\_aligned\_joint\_ancestors\_ConfValue"

1.  dicistrovirus II

java -jar bnkit.jar -aln dicistrovirus\_II\_wo\_outgroup\_aligned.fasta -nwk dicistrovirus\_II\_wo\_outgroup\_w\_root.tree -out dicistrovirus\_II\_wo\_outgroup\_aligned\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "dicistrovirus\_II\_wo\_outgroup\_aligned\_joint\_ancestors\_ConfValue"
