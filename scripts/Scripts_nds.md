Scripts\_New2021
================

New Dataset (2021)
------------------

First part: Flaviridae

<file:%22flaviviridae_w_outgroup_wo_flavivirus_aligned.fasta>"

Ophiostoma mitovirus 4 (Narnaviridae family) is the outgroup

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s flaviviridae\_w\_outgroup\_wo\_flavivirus\_aligned.fasta -n T1

This command creates the output file: "flaviviridae\_w\_outgroup\_wo\_flavivirus\_woroot.tree"

``` r
library(ape)
library(phylotools)
data=read.fasta("flaviviridae_w_outgroup_wo_flavivirus_aligned.fasta")

# Root at outgroup

tree.1=read.tree("flaviviridae_w_outgroup_wo_flavivirus_woroot.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")
index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"flaviviridae_wo_outgroup_wo_flavivirus_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"flaviviridae_wo_outgroup_wo_flavivirus_aligned.fasta")
```

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t flaviviridae\_wo\_outgroup\_wo\_flavivirus\_w\_root.tree -s flaviviridae\_wo\_outgroup\_wo\_flavivirus\_aligned.fasta -m GTRGAMMA -n S1

This command product output file: "flaviviridae\_wo\_outgroup\_wo\_flavivirus\_AncestralSeq.fasta"

Second part: Picornviridae

10 aligned sequences Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s picornaviridae\_w\_outgroup\_aligned.fasta -n T2

This command creates the output file: "picornaviridae\_w\_outgroup\_woroot.tree"

``` r
rm(list=ls())
library(ape)
library(phylotools)
data=read.fasta("picornaviridae_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("picornaviridae_w_outgroup_woroot.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")


index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"picornaviridae_wo_outgroup_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"picornaviridae_wo_outgroup_aligned.fasta")
```

Note: the gaps here are removed

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t picornaviridae\_wo\_outgroup\_w\_root.tree -s picornaviridae\_wo\_outgroup\_aligned.fasta -m GTRGAMMA -n S2

This command product output file: "picornaviridae\_wo\_outgroup\_aligned\_AncestralSeq.fasta"

Third part: Potyviridae and Poacevirus

Ophiostoma mitovirus 4 (Narnaviridae family) is added as outgroup.

Construct the phylogenetic tree:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 500 -s potyviridae\_poacevirus\_w\_outgroup\_aligned.fasta -n T3

This command creates the output file: "potyviridae\_poacevirus\_w\_outgroup\_woroot.tree"

``` r
rm(list=ls())
library(ape)
library(phylotools)
data=read.fasta("potyviridae_poacevirus_w_outgroup_aligned.fasta")

# Root at outgroup

tree.1=read.tree("potyviridae_poacevirus_w_outgroup_woroot.tree")
tree.1_unroot=unroot(tree.1)
tree.2=root(tree.1_unroot,"NC004052.1:1-207Ophiostomamitovirus4completegenome")

index=which(tree.2[["tip.label"]]=="NC004052.1:1-207Ophiostomamitovirus4completegenome") 
name_without_outgroup=tree.2[["tip.label"]][-(index)]
tree.3=keep.tip(tree.2,name_without_outgroup)
#write.tree(tree.3,"potyviridae_poacevirus_wo_outgroup_w_root.tree")

# Remove the name in the ancestral sequences
index2=which(data$seq.name=="NC_004052.1:1-207_Ophiostoma_mitovirus_4_complete_genome")
seq_deletion=data[-(index2),]
#dat2fasta(seq_deletion,"potyviridae_poacevirus_wo_outgroup_aligned.fasta")
```

##### Ancestral Reconstruction (the one at the rooting point)

raxmlHPC -f A -p 23456 -t potyviridae\_poacevirus\_w\_outgroup\_w\_root.tree -s potyviridae\_poacevirus\_w\_outgroup\_aligned.fasta -m GTRGAMMA -n S3

This command product output file: "potyviridae\_poacevirus\_w\_outgroup\_AncestralSeq.fasta"

##### Constructing the ancestral sequences at every internal nodes

Note: The method I use could only generatet the confidence value at one internal nodes. I will combine the confidence value into a single file. (e.g, cat *ConfValue* &gt;&gt; combined.txt) (N0,N1,...,)

1.Flaviridae

java -jar bnkit.jar -aln flaviviridae\_wo\_outgroup\_wo\_flavivirus\_aligned.fasta -nwk flaviviridae\_wo\_outgroup\_wo\_flavivirus\_w\_root.tree -out flaviviridae\_wo\_outgroup\_wo\_flavivirus\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "flaviviridae\_wo\_outgroup\_wo\_flavivirus\_joint\_ancestors\_ConfValue"

1.  Picornviridae

java -jar bnkit.jar -aln picornaviridae\_wo\_outgroup\_aligned.fasta -nwk picornaviridae\_wo\_outgroup\_w\_root.tree -out picornaviridae\_wo\_outgroup\_aligned\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "picornaviridae\_wo\_outgroup\_aligned\_joint\_ancestors\_ConfValue"

1.  Potyviridae and Poacevirus

java -jar bnkit.jar -aln potyviridae\_poacevirus\_wo\_outgroup\_aligned.fasta -nwk potyviridae\_poacevirus\_wo\_outgroup\_w\_root.tree -out potyviridae\_poacevirus\_wo\_outgroup\_joint\_ancestors\_ConfValue -marg N0 -format DISTRIB -model Yang

This step produces files: "potyviridae\_poacevirus\_wo\_outgroup\_joint\_ancestors\_ConfValue.fasta"
