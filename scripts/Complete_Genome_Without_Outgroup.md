Complete\_Genome\_Mitovirus\_Without\_Outgroup
================

In this part, we want to construct the ancestral sequences and states with the complete genome without the outgroup. The outgroup is Mitovirus.

#### Picornavirus

``` r
rm(list=ls())
library(phylotools)
```

    ## Loading required package: ape

``` r
library(ape)
seq_picorn=read.fasta("picornaviridae_genome_mitovirus_outgroup_aligned_nospacegap_deleted.fas")
tree_picorn_outgroup=read.tree("replicase_picornavirus_nucleotide_mito_outgrouprooting_ID.tree")

# Remove the outgroup in the tree

index1=which(tree_picorn_outgroup[["tip.label"]]=="NC004052")
index2=which(tree_picorn_outgroup[["tip.label"]]=="NC004046")
name_without_outgroup=tree_picorn_outgroup[["tip.label"]][-c(index1,index2)]
tree_picorn_outgroup_deletion=keep.tip(tree_picorn_outgroup,name_without_outgroup)
#write.tree(tree_picorn_outgroup_deletion,"replicase_picornavirus_nucleotide_mito_outgrouprooting_ID_nooutgroup.tree")

# Remove the outgroup in the sequences file

index3=which(seq_picorn$seq.name=="NC004046")
index4=which(seq_picorn$seq.name=="NC004052")
seq_picorn_deletion=seq_picorn[-c(index3,index4),]
#dat2fasta(seq_picorn_deletion,"picornaviridae_genome_mitovirus_aligned_nospacegap_deleted_nooutgroup.fas")
```

Before the ancestral reconstruction, the gap-only columns are removed using MEGAX Special characters are treated as the missing (-) at each site.

##### Ancestral Sequences Reconstruction using the nucleotide files (with RAxML):

raxmlHPC -f A -p 12345 -t replicase\_picornavirus\_nucleotide\_mito\_outgrouprooting\_ID\_nooutgroup.tree -s picornaviridae\_genome\_mitovirus\_aligned\_nospacegap\_deleted\_nooutgroupgap.fas -m GTRGAMMA -n S1

This command produces output file: "picornaviridae\_genome\_mitovirus\_aligned\_outgrouprooting\_nooutgroup\_AncestralSeq.fasta"

##### The confidence values of most probable sequences reconstruction (with JAVA):

java -jar bnkit.jar -aln picornaviridae\_genome\_mitovirus\_aligned\_nospacegap\_deleted\_nooutgroupgap.fas -nwk replicase\_picornavirus\_nucleotide\_mito\_outgrouprooting\_ID\_nooutgroup.tree -out picornaviridae\_genome\_mitovirus\_aligned\_outgrouprooting\_nooutgroup\_Ancestral\_ConfValue -marg N0 -format DISTRIB -model Yang

This command produces the following file: "picornaviridae\_genome\_mitovirus\_aligned\_outgrouprooting\_nooutgroup\_Ancestral\_ConfValue"

#### Potyvirus

The outgroup is Mitovirus

``` r
rm(list=ls())
library(phylotools)
library(ape)
seq_poty=read.fasta("potyviridae_genome_mitovirus_outgroup_aligned_nospace.fasta")
tree_poty_outgroup=read.tree("potyviridae_genome_mitovirus_outgroup_outgrouprooting.tree")

# Remove the outgroup in the tree

index1=which(tree_poty_outgroup[["tip.label"]]=="NC004052.1_Ophiostoma_mitovirus_4_complete_genome")
index2=which(tree_poty_outgroup[["tip.label"]]=="NC004046.1_Cryphonectria_parasitica_mitovirus_1-NB631_complete_genome")
name_without_outgroup=tree_poty_outgroup[["tip.label"]][-c(index1,index2)]
tree_poty_outgroup_deletion=keep.tip(tree_poty_outgroup,name_without_outgroup)
#write.tree(tree_poty_outgroup_deletion,"potyviridae_genome_mitovirus_outgrouprooting_nooutgroup.tree")

# Remove the outgroup in the sequences file

index3=which(seq_poty$seq.name=="NC004052.1_Ophiostoma_mitovirus_4_complete_genome")
index4=which(seq_poty$seq.name=="NC004046.1_Cryphonectria_parasitica_mitovirus_1-NB631_complete_genome")
seq_poty_deletion=seq_poty[-c(index3,index4),]
#dat2fasta(seq_poty_deletion,"potyviridae_genome_mitovirus_aligned_nospace_nooutgroup.fas")
```

Before the ancestral reconstruction, the gap-only columns are removed using MEGAX Special characters are treated as the missing (-) at each site.

##### Ancestral Sequences Reconstruction using the nucleotide files (with RAxML):

raxmlHPC -f A -p 12345 -t potyviridae\_genome\_mitovirus\_outgrouprooting\_nooutgroup.tree -s potyviridae\_genome\_mitovirus\_aligned\_nospacegap\_nooutgroup.fas -m GTRGAMMA -n S2

This command produces output file: "potyviridae\_genome\_mitovirus\_aligned\_outgrouprooting\_nooutgroup\_AncestralSeq.fasta"

##### The confidence values of most probable sequences reconstruction (with JAVA):

java -jar bnkit.jar -aln potyviridae\_genome\_mitovirus\_aligned\_nospacegap\_nooutgroup.fas -nwk potyviridae\_genome\_mitovirus\_outgrouprooting\_nooutgroup.tree -out potyviridae\_genome\_mitovirus\_aligned\_outgrouprooting\_nooutgroup\_Ancestral\_ConfValue -marg N0 -format DISTRIB -model Yang

This command produces the following file: "potyviridae\_genome\_mitovirus\_aligned\_outgrouprooting\_nooutgroup\_Ancestral\_ConfValue"

###### Note:

This part with the confidence values of Potyviridae doesn't work because of my previous method has limitation on the memory. I am trying to get this result with an alternative approach.

#### Picornavirus + Potyvirus

``` r
rm(list=ls())
library(phylotools)
library(ape)
seq_picorn_poty=read.fasta("picornaviridae_potyviridae_genome_mitovirus_outgroup_aligned_nospacegap_deleted.fas")
tree_picorn_poty_outgroup=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_mito_outgroup_ID.tree")

# Remove the outgroup in the tree

index1=which(tree_picorn_poty_outgroup[["tip.label"]]=="NC004052")
index2=which(tree_picorn_poty_outgroup[["tip.label"]]=="NC004046")
name_without_outgroup=tree_picorn_poty_outgroup[["tip.label"]][-c(index1,index2)]
tree_picorn_poty_outgroup_deletion=keep.tip(tree_picorn_poty_outgroup,name_without_outgroup)
#write.tree(tree_picorn_poty_outgroup_deletion,"replicase_picornavirus_potyvirus_nucleotide_mitovirus_outgrouprooting_ID_nooutgroup.tree")

# Remove the outgroup in the sequences file

index3=which(seq_picorn_poty$seq.name=="NC004046")
index4=which(seq_picorn_poty$seq.name=="NC004052")
seq_picorn_poty_deletion=seq_picorn_poty[-c(index3,index4),]
#dat2fasta(seq_picorn_poty_deletion,"picornaviridae_potyviridae_genome_mitovirus_nospacegap_deleted_nooutgroup.fas")
```

Before the ancestral reconstruction, the gap-only columns are removed using MEGAX Special characters are treated as the missing (-) at each site.

##### Ancestral Sequences Reconstruction using the nucleotide files (with RAxML):

raxmlHPC -f A -p 12345 -t replicase\_picornavirus\_potyvirus\_nucleotide\_mitovirus\_outgrouprooting\_ID\_nooutgroup.tree -s picornaviridae\_potyviridae\_genome\_mitovirus\_nospacegap\_deleted\_nooutgroup.fas -m GTRGAMMA -n S3

This command produces output file: "picornaviridae\_potyviridae\_genome\_mitovirus\_outgrouprooting\_nooutgroup\_AncestralSeq.fasta"

##### The confidence values of most probable sequences reconstruction (with JAVA):

java -jar bnkit.jar -aln picornaviridae\_potyviridae\_genome\_mitovirus\_nospacegap\_deleted\_nooutgroup.fas -nwk replicase\_picornavirus\_potyvirus\_nucleotide\_mitovirus\_outgrouprooting\_ID\_nooutgroup.tree -out picornaviridae\_potyviridae\_genome\_mitovirus\_outgrouprooting\_nooutgroup\_Ancestral\_ConfValue -marg N0 -format DISTRIB -model Yang

This command produces the following file: "picornaviridae\_potyviridae\_genome\_mitovirus\_outgrouprooting\_nooutgroup\_Ancestral\_ConfValue"

###### Note:

This part with the confidence values of Potyviridae doesn't work because of my previous method has limitation on the memory. I am trying to get this result with an alternative approach.
