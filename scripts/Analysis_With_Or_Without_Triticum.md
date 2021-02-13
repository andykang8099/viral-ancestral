Analysis\_Without\_Outgroup\_Triticum
================

Section 1
=========

Now, construct the ancestral reconstruction with the picornviridae, but the outgroup is now converted to the Triticum virus

##### Picornaviridae Nucleotide and Protein files

###### Part 1. Construct the Ancestral Sequences without the outgroup (mitovirus), keeping the Triticum

###### Nucleotide files

``` r
rm(list=ls())
library(phylotools)
```

    ## Loading required package: ape

``` r
library(ape)
picorna_tri_5UTR_mito=read.fasta("picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID_nospacegap.fas")
tree_rep_picorna_tri_nucle_mito_or=read.tree("picornavirus_triticum_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")

# Remove the outgroup mitovirus in the tree

index1=which(tree_rep_picorna_tri_nucle_mito_or[["tip.label"]]=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index2=which(tree_rep_picorna_tri_nucle_mito_or[["tip.label"]]=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
name_without_outgroup=tree_rep_picorna_tri_nucle_mito_or[["tip.label"]][-c(index1,index2)]
tree_rep_picorna_tri_nucle_mito_or_deletion=keep.tip(tree_rep_picorna_tri_nucle_mito_or,name_without_outgroup)
#write.tree(tree_rep_picorna_tri_nucle_mito_or_deletion,"picornavirus_triticum_nucleotide_mitovirus_replicasetree_outgrouprooting_nooutgroup.tree")

# Remove the outgroup in the name

index3=which(picorna_tri_5UTR_mito$seq.name=="NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome")
index4=which(picorna_tri_5UTR_mito$seq.name=="NC004052.1-205-2556Ophiostomamitovirus4completegenome")
picorna_tri_5UTR_mito_deletion=picorna_tri_5UTR_mito[-c(index3,index4),]
#dat2fasta(picorna_tri_5UTR_mito_deletion,"picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID_nooutgroup.fas")
```

Ancestral reconstruction on nucleotide files:

raxmlHPC -f A -p 12345 -t picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -s picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nooutgroup.fas -m GTRGAMMA -n S1

This command product output file: "picornaviridae\_triticum\_5UTR\_mitovirus\_nucleotide\_outgrouprooting\_with\_Triticum\_nooutgroup\_AncestralSeq.fasta"

###### Protein files

``` r
rm(list=ls());library(phylotools);library(ape)
picorna_tri_protein_mito=read.fasta("replicase_picornavirus_triticum_protein_mitovirus_outgroup_nospaceX.fas")
tree_rep_picorna_tri_pro_mito_or=read.tree("replicase_picornavirus_triticum_protein_mitovirus_outgrouprooting.tree")

# Remove the outgroup in the tree

index1=which(tree_rep_picorna_tri_pro_mito_or[["tip.label"]]=="NP660174.1putativeRNA-dependentRNApolymerase_Cryphonectriaparasiticamitovirus1-NB631")
index2=which(tree_rep_picorna_tri_pro_mito_or[["tip.label"]]=="NP660179.1RNA-dependentRNApolymeraseputative_Ophiostomamitovirus4")
name_without_outgroup=tree_rep_picorna_tri_pro_mito_or[["tip.label"]][-c(index1,index2)]
tree_rep_picorna_tri_pro_mito_or_deletion=keep.tip(tree_rep_picorna_tri_pro_mito_or,name_without_outgroup)
#write.tree(tree_rep_picorna_tri_pro_mito_or_deletion,"replicase_picornavirus_triticum_protein_mitovirus_outgrouprooting_nooutgroup.tree")

# Remove the outgroup in the name

index3=which(picorna_tri_protein_mito$seq.name=="NP660174.1putativeRNA-dependentRNApolymerase_Cryphonectriaparasiticamitovirus1-NB631")
index4=which(picorna_tri_protein_mito$seq.name=="NP660179.1RNA-dependentRNApolymeraseputative_Ophiostomamitovirus4")
picorna_tri_protein_mito_deletion=picorna_tri_protein_mito[-c(index3,index4),]
#dat2fasta(picorna_tri_protein_mito_deletion,"replicase_picornavirus_triticum_protein_mitovirus_outgroup_nooutgroupfile.fas")
```

Note: the gap-only sites are also deleted

Ancestral reconstruction on nucleotide files:

raxmlHPC -f A -p 23456 -t replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting\_nooutgroup.tree -s replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgroup\_nooutgroupfile.fas -m PROTGAMMAJTT -n S2

This command product output file: "replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting\_with\_Triticum\_nooutgroup\_AncestralSeq.fasta"

###### Part 2. Construct the Ancestral Sequences without the Triticum and Outgroup (delete mitovirus and Triticum)

###### Nucleotide files

``` r
rm(list=ls());library(phylotools);library(ape)
seq1=read.fasta("picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID_nooutgroup.fas")
tree1=read.tree("picornavirus_triticum_nucleotide_mitovirus_replicasetree_outgrouprooting_nooutgroup.tree")

# Remove the Triticum in the tree

index1=which(tree1[["tip.label"]]=="NC012799TriticummosaicviruscompletegenomeNIbreplicase-77259195")
name_without_tri=tree1[["tip.label"]][-(index1)]
tree1_deletion=keep.tip(tree1,name_without_tri)
#write.tree(tree1_deletion,"picornavirus_triticum_nucleotide_mitovirus_replicasetree_outgrouprooting_noTriticum.tree")

# Remove the Triticum in the sequences name

index2=which(seq1$seq.name=="NC012799TriticummosaicviruscompletegenomeNIbreplicase-77259195")
seq1_deletion=seq1[-(index2),]
#dat2fasta(seq1_deletion,"picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID_noTriticum.fas")
```

Note: the gap-only sites are deleted

Ancestral reconstruction on nucleotide files:

raxmlHPC -f A -p 12345 -t picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_outgrouprooting\_noTriticum.tree -s picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_noTriticum.fas -m GTRGAMMA -n S3

This command product output file: "picornaviridae\_triticum\_5UTR\_mitovirus\_nucleotide\_outgrouprooting\_noTriticum\_noOutgroup\_AncestralSeq.fasta"

###### Protein files

``` r
rm(list=ls());library(phylotools);library(ape)
seq2=read.fasta("replicase_picornavirus_triticum_protein_mitovirus_outgroup_nooutgroupfile.fas")
tree2=read.tree("replicase_picornavirus_triticum_protein_mitovirus_outgrouprooting_nooutgroup.tree")

# Remove the Triticum in the tree

index1=which(tree2[["tip.label"]]=="rf1NC012799TriticummosaicviruscompletegenomeNIbreplicase-77259195")
name_without_tri=tree2[["tip.label"]][-(index1)]
tree2_deletion=keep.tip(tree2,name_without_tri)
#write.tree(tree2_deletion,"replicase_picornavirus_triticum_protein_mitovirus_outgrouprooting_noTriticum.tree")

# Remove the Triticum in the name

index2=which(seq2$seq.name=="rf1NC012799TriticummosaicviruscompletegenomeNIbreplicase-77259195")
seq2_deletion=seq2[-(index2),]
#dat2fasta(seq2_deletion,"replicase_picornavirus_triticum_protein_mitovirus_outgroup_noTriticumfile.fas")
```

Note: the gap-only sites are also deleted

Ancestral reconstruction on the protein files:

raxmlHPC -f A -p 23456 -t replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting\_noTriticum.tree -s replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgroup\_noTriticumfile.fas -m PROTGAMMAJTT -n S4

This command produces output file: "replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting\_noTriticum\_noOutgroup\_AncestralSeq.fasta"

Part 3: Ancestral Confidence value Construction
-----------------------------------------------

### Picornaviridae

1.  Without outgroup, keep the Triticum

java -jar bnkit.jar -aln picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nooutgroup.fas -nwk picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_outgrouprooting\_nooutgroup.tree -out picornaviridae\_triticum\_5UTR\_mitovirus\_nucleotide\_outgrouprooting\_With\_Triticum\_nooutgroup\_Ancestral\_Confidence -marg N0 -format DISTRIB -model Yang

This step produces files: "picornaviridae\_triticum\_5UTR\_mitovirus\_nucleotide\_outgrouprooting\_With\_Triticum\_nooutgroup\_Ancestral\_Confidence"

1.  Without outgroup and Triticum

java -jar bnkit.jar -aln picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_noTriticum\_nogap.fas -nwk picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_outgrouprooting\_noTriticum.tree -out picornaviridae\_triticum\_5UTR\_mitovirus\_nucleotide\_outgrouprooting\_noTriticum\_noOutgroup\_Ancestral\_Confidence -marg N0 -format DISTRIB -model Yang

This step produces files:"picornaviridae\_triticum\_5UTR\_mitovirus\_nucleotide\_outgrouprooting\_noTriticum\_noOutgroup\_Ancestral\_Confidence"
