Analysis\_Without\_Outgroup
================

#### Section 1

##### 5UTR potyviridae (replicase 2 tree rooted at outgroup)

Part1:

original file: "potyviridae\_5UTR\_mitovirus\_outgroup\_aligned.fasta" original tree: "replicase\_potyvirus\_nucleotide\_mito\_outgroup.fasta"

``` r
library(phylotools)
```

    ## Loading required package: ape

``` r
library(ape)
poty_5UTR_mito=read.fasta("potyviridae_5UTR_mitovirus_outgroup_ReplicaseID_nospacegap.fas")
tree_rep_poty_mito_or=read.tree("potyvirus_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")

# Remove the outgroup in the tree

index1=which(tree_rep_poty_mito_or[["tip.label"]]=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index2=which(tree_rep_poty_mito_or[["tip.label"]]=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
name_without_outgroup=tree_rep_poty_mito_or[["tip.label"]][-c(index1,index2)]
tree_rep_poty_mito_or_deletion=keep.tip(tree_rep_poty_mito_or,name_without_outgroup)
#write.tree(tree_rep_poty_mito_or_deletion,"potyvirus_nucleotide_mitovirus_replicasetree_outgrouprooting_nooutgroup.tree")

# Remove the outgroup in the name

index3=which(poty_5UTR_mito$seq.name=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index4=which(poty_5UTR_mito$seq.name=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
poty_5UTR_mito_deletion=poty_5UTR_mito[-c(index3,index4),]
#dat2fasta(poty_5UTR_mito_deletion,"potyviridae_5UTR_mitovirus_outgroup_ReplicaseID_nooutgroup.fas")
```

Part 2: Ancestral Reconstruction

raxmlHPC -f A -p 12345 -t potyvirus\_nucleotide\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -s potyviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nooutgroup.fas -m GTRGAMMA -n S1

This command product output file: "potyviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup\_AncestralSeq.fasta"

#### Section 2

##### 5UTR picornaviridae + potyviridae (replicase 2 tree rooted at outgroup)

original file: "potyviridae\_picornaviridae\_5UTR\_mitovirus\_outgroup\_aligned.fasta.fasta" original tree: "replicase\_picornavirus\_and\_potyvirus\_nucleotide\_mito\_outgroup.tree"

``` r
picorna_poty_5UTR_mito=read.fasta("potyviridae_picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID_nospacegap.fas")
tree_rep_picorna_poty_mito_or=read.tree("picornavirus_and_potyvirus_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")

# Remove the outgroup in the tree

index1=which(tree_rep_picorna_poty_mito_or[["tip.label"]]=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index2=which(tree_rep_picorna_poty_mito_or[["tip.label"]]=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
name_without_outgroup=tree_rep_picorna_poty_mito_or[["tip.label"]][-c(index1,index2)]
tree_rep_picorna_poty_mito_or_deletion=keep.tip(tree_rep_picorna_poty_mito_or,name_without_outgroup)
#write.tree(tree_rep_picorna_poty_mito_or_deletion,"potyviridae_picornaviridae_5UTR_mitovirus_replicasetree_outgrouprooting_nooutgroup.tree")

# Remove the outgroup in the name

index3=which(picorna_poty_5UTR_mito$seq.name=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index4=which(picorna_poty_5UTR_mito$seq.name=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
picorna_poty_5UTR_mito_deletion=picorna_poty_5UTR_mito[-c(index3,index4),]
#dat2fasta(picorna_poty_5UTR_mito_deletion,"potyviridae_picornaviridae_5UTR_mitovirus_ReplicaseID_nooutgroup.fas")
```

Part 2: Ancestral Reconstruction

raxmlHPC -f A -p 12345 -t potyviridae\_picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -s potyviridae\_picornaviridae\_5UTR\_mitovirus\_ReplicaseID\_nooutgroup.fas -m GTRGAMMA -n S2

This command product output file: "potyviridae\_picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup\_AncestralSeq.fasta"

#### Section 3

##### 5UTR picornaviridae (replicase 2 tree rooted at outgroup)

original file: "picornaviridae\_5UTR\_mitovirus\_outgroup\_aligned.fasta" original tree: "replicase\_picornavirus\_nucleotide\_mito\_outgrouprooting.tree"

``` r
picorna_5UTR_mito=read.fasta("picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID_nospacegap.fas")
tree_rep_picorna_mito_or=read.tree("picornavirus_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")

# Remove the outgroup in the tree

index1=which(tree_rep_picorna_mito_or[["tip.label"]]=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index2=which(tree_rep_picorna_mito_or[["tip.label"]]=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
name_without_outgroup=tree_rep_picorna_mito_or[["tip.label"]][-c(index1,index2)]
tree_rep_picorna_mito_or_deletion=keep.tip(tree_rep_picorna_mito_or,name_without_outgroup)
#write.tree(tree_rep_picorna_mito_or_deletion,"picornaviridae_5UTR_mitovirus_replicasetree_outgrouprooting_nooutgroup.tree")

# Remove the outgroup in the name

index3=which(picorna_5UTR_mito$seq.name=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index4=which(picorna_5UTR_mito$seq.name=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
picorna_5UTR_mito_deletion=picorna_5UTR_mito[-c(index3,index4),]
#dat2fasta(picorna_5UTR_mito_deletion,"picornaviridae_5UTR_mitovirus_ReplicaseID_nooutgroup.fas")
```

Part 2: Ancestral Reconstruction

raxmlHPC -f A -p 12345 -t picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -s picornaviridae\_5UTR\_mitovirus\_ReplicaseID\_nooutgroup.fas -m GTRGAMMA -n S3

This command product output file: "picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup\_AncestralSeq.fasta"

Section 4
=========

Now, after constructing the ancestral sequences, we want to construct the confidence value of of the inferred ancestral sequence

We should note that the inferred ancestral sequence derives from the most probable character in the maarginal reconstruction, and now the confidence value (support value) of four nucleotide bases are displayed in a table

The following code is operated in the Linux (Terminal)

##### 5UTR potyviridae (replicase 2 tree rooted at outgroup)

java -jar bnkit.jar -aln potyviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nooutgroup.fas -nwk potyvirus\_nucleotide\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -out potyviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup\_Ancestral\_Confidence -marg N0 -format DISTRIB -model Yang

##### 5UTR picornaviridae + potyviridae (replicase 2 tree rooted at outgroup)

java -jar bnkit.jar -aln potyviridae\_picornaviridae\_5UTR\_mitovirus\_ReplicaseID\_nooutgroup.fas -nwk potyviridae\_picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -out potyviridae\_picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup\_Ancestral\_Confidence -marg N0 -format DISTRIB -model Yang

##### 5UTR picornaviridae (replicase 2 tree rooted at outgroup)

java -jar bnkit.jar -aln picornaviridae\_5UTR\_mitovirus\_ReplicaseID\_nooutgroup.fas -nwk picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -out picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup\_Ancestral\_Confidence -marg N0 -format DISTRIB -model Yang
