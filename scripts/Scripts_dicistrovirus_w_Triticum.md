Scripts\_dicistrovirus\_w\_Triticum
================

New Dataset (Supplemental data dicistrovirus with Triticum and Picornaviridae)
------------------------------------------------------------------------------

*First part*: dicistrovirus I (with Triticum)

<file:%22Dicistrovirus_I_w_triticum_w_outgroup_aligned.fasta>": 7 aligned nucleotide sequences of the 5' UTR of members of Dicistroviridae family with a characterized IRES, plus Triticum mosaic virus (TriMV), plus 1 outgroup (9 sequence total: 4 Cripavirus, 3 Aparavirus, 1 TriMV, 1 outgroup). Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s Dicistrovirus\_I\_w\_triticum\_w\_outgroup\_aligned.fasta -n T1

This command creates the output file: "Dicistrovirus\_I\_w\_triticum\_w\_outgroup\_aligned\_wo\_root.tree"

``` r
library(ape)
library(phylotools)
data=read.fasta("Dicistrovirus_I_w_triticum_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("Dicistrovirus_I_w_triticum_w_outgroup_aligned_wo_root.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")
index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"Dicistrovirus_I_w_triticum_wo_outgroup_aligned_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"Dicistrovirus_I_w_triticum_wo_outgroup_aligned.fasta")
```

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_aligned\_w\_root.tree -s Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_aligned.fasta -m GTRGAMMA -n S1

This command produces output file: "Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_AncestralSeq.fasta"

\*\*Joint ancestors file: "Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_joint\_ancestors.fasta"

\*\*Tree with label file: "Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_w\_root\_w\_label.tree"

*Second part*: dicistrovirus II with Triticum <file:%22Dicistrovirus_II_w_triticum_w_outgroup_aligned.fasta>"

7 aligned nucleotide sequences of the intergenic region of members of Dicistroviridae family with a characterized IRES, plus Triticum mosaic virus (TriMV), plus 1 outgroup (9 sequence total: 4 Cripavirus, 3 Aparavirus, 1 TriMV, 1 outgroup). Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s Dicistrovirus\_II\_w\_triticum\_w\_outgroup\_aligned.fasta -n T2

This command creates the output file: "Dicistrovirus\_II\_w\_triticum\_w\_outgroup\_aligned\_wo\_root.tree"

``` r
rm(list=ls())
library(ape)
library(phylotools)
data=read.fasta("Dicistrovirus_II_w_triticum_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("Dicistrovirus_II_w_triticum_w_outgroup_aligned_wo_root.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")


index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"Dicistrovirus_II_w_triticum_wo_outgroup_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"Dicistrovirus_II_w_triticum_wo_outgroup_aligned.fasta")
```

Note: the gaps here are removed

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_aligned\_w\_root.tree -s Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_aligned.fasta -m GTRGAMMA -n S2

This command produces output file: "Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_AncestralSeq.fasta"

\*\*Joint ancestors file: "Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_joint\_ancestors.fasta"

\*\*Tree with label file: "Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_w\_root\_w\_label.tree"

*Third part*: Picornviridae with Triticum <file:%22Picornaviridae_w_triticum_w_outgroup_aligned.fasta>"

10 aligned nucleotide sequences of the 5' UTR of members of Picornaviridae family with a characterized IRES, plus Triticum mosaic virus (TriMV), plus 1 outgroup (12 sequence total: 5 Enterovirus, 2 Aphotvirus, 1 Hepatovirus, 2 Cardiovirus, 1 TriMV, 1 outgroup). Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s Picornaviridae\_w\_triticum\_w\_outgroup\_aligned.fasta -n T3

This command creates the output file: "Picornaviridae\_w\_triticum\_w\_outgroup\_aligned\_wo\_root.tree"

``` r
rm(list=ls())
library(ape)
library(phylotools)
data=read.fasta("Picornaviridae_w_triticum_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("Picornaviridae_w_triticum_w_outgroup_aligned_wo_root.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")


index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"Picornaviridae_w_triticum_wo_outgroup_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"Picornaviridae_w_triticum_wo_outgroup_aligned.fasta")
```

Note: the gaps here are removed

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t Picornaviridae\_w\_triticum\_wo\_outgroup\_w\_root.tree -s Picornaviridae\_w\_triticum\_wo\_outgroup\_aligned.fasta -m GTRGAMMA -n S3

This command produces output file: "Picornaviridae\_w\_triticum\_wo\_outgroup\_AncestralSeq.fasta"

\*\*Joint ancestors file: "Picornaviridae\_w\_triticum\_wo\_outgroup\_joint\_ancestors.fasta"

\*\*Tree with label file: "Picornaviridae\_w\_triticum\_wo\_outgroup\_w\_root\_w\_label.tree"

##### Constructing the ancestral sequences and their confidence value at every internal nodes (sequences with Triticum

Note: The method I use could only generatet the confidence value at one internal nodes. I will combine the confidence value into a single file. (e.g, cat *ConfValue* &gt;&gt; combined.txt) (N0,N1,...,)

1.  dicistrovirus I (with Triticum)

java -jar bnkit.jar -aln Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_aligned.fasta -nwk Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_aligned\_w\_root.tree -out Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "Dicistrovirus\_I\_w\_triticum\_wo\_outgroup\_joint\_ancestors\_ConfValue"

1.  dicistrovirus II with Triticum

java -jar bnkit.jar -aln Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_aligned.fasta -nwk dicistrovirus\_II\_wo\_outgroup\_w\_root.tree -out Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "Dicistrovirus\_II\_w\_triticum\_wo\_outgroup\_joint\_ancestors\_ConfValue"

1.  Picornviridae with Triticum

java -jar bnkit.jar -aln Picornaviridae\_w\_triticum\_wo\_outgroup\_aligned.fasta -nwk Picornaviridae\_w\_triticum\_wo\_outgroup\_w\_root.tree -out Picornaviridae\_w\_triticum\_wo\_outgroup\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "Picornaviridae\_w\_triticum\_wo\_outgroup\_joint\_ancestors\_ConfValue"
