replicase\_triticum\_outgroup
================

This part we also focus on the the replicase file of picornaviridae, but the outgroup is now triticum
-----------------------------------------------------------------------------------------------------

Section 1
---------

This parts show how to remove the space and , in the ID of each virus because RAXML cannot discern space and comma when RAXML try to build the tree.

##### picornviridae

``` r
library(phylotools)
```

    ## Loading required package: ape

``` r
library(stringr)
seq=read.fasta("replicase_picornavirus_triticum_nucleotide_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_picornavirus_triticum_nucleotide_mitovirus_outgroup_nospace.fasta")

# Convert the name of ID of the protein file
seq1=read.fasta("replicase_picornavirus_triticum_protein_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}

#dat2fasta(seq1, outfile = "replicase_picornavirus_triticum_protein_mitovirus_outgroup_nospace.fasta")
```

Two files are generated: "replicase\_picornavirus\_triticum\_nucleotide\_mitovirus\_outgroup\_nospace.fasta" "replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgroup\_nospace.fasta"

Section 2
---------

This part matches the ID between the replicase files and the 5UTR files exactly. The reason is to make sure that the ancestral sequences could be constructed successfully.

The difference between the ID of the repicase file and the 5UTR file is that the replicase file contains the location, and the ID is more specified.

##### picornavirus

``` r
rm(list=ls())
library(phylotools)
library(stringr)

# 132 virus
seq_picorna_trit_5UTR_mito=read.fasta("picornaviridae_triticum_5UTR_mitovirus_outgroup_aligned.fasta")

# 87 virus
seq_replicase_picorna_trit_nucleotide_mito=read.fasta("replicase_picornavirus_triticum_nucleotide_mitovirus_outgroup_aligned.fasta")
replicase_picorna_trit_nucleotide_mito_name=seq_replicase_picorna_trit_nucleotide_mito$seq.name

seq=seq_picorna_trit_5UTR_mito
name=replicase_picorna_trit_nucleotide_mito_name

m=0
for (i in 1:nrow(seq)) {
  name_seq=strsplit(seq[i,1], split=" ")[[1]][1]
  name_seqid=str_extract(name_seq, pattern = "\\d+")
   for ( q in 1:length(name)) {
    if (grepl(name_seqid, name[q],fixed=TRUE)==TRUE) {
      seq[i,1]=name[q]
      m=1
    }
   }
  if (m==0) {
    seq[i,1]="not found"
  }
  m=0
}

# Delete extra virus that only in 5UTR
index=which(seq[,1]=="not found")
seq=seq[-(index),]
# Create a new row name
rownames(seq)=1:nrow(seq)

# now, the sequence has 73 virus, less than 87

# Detect and display the extra virus in replicase but not in 5UTR

seq1=seq_replicase_picorna_trit_nucleotide_mito
for ( i in 1: nrow(seq1)) {
   if (sum(seq1[i,1]==seq[,1])==0) {
     seq1[i,1]="not match 5UTR"
   }
}
index=which(seq1[,1]=="not match 5UTR")
seq1=seq1[-(index),]

# Create a new row name
rownames(seq1)=1:nrow(seq1)

for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}

#dat2fasta(seq1, outfile = "replicase_picornavirus_triricum_nucleotide_5UTR_mito_combination.fasta")
#dat2fasta(seq, outfile = "picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
```

Now, there are 73 virus sequences.

This command generates: "picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID.fasta"

Section 3
---------

Now, construct the phylogenetic tree based on the aligned sequence alignment files derived in the first section

Each tree file is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF. According to some research, JTTF model is more reliable in the construction of phylogenetic tree of protein sequence.

Note: In order to make the output files more clear and display it correctly, I change T1 and T2 to .tree in the final file names.

In all the RAXML commands, I use T1 and T2 to differentiate two tree files in order to use a uniform name to show the tree successfully in the phylogenetic tree software.The first part of the name of the files remains same (output name in the RAXML).

### Nucleotide files

##### picornavirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_triticum\_nucleotide\_mitovirus\_outgroup\_nospace.fasta -n T1

This command creates the output file: "replicase\_picornavirus\_triticum\_nucleotide\_mitovirus\_outgroup\_noroot.tree"

### Protein files

##### picornavirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_picorna\_protein\_mito\_nospace.fasta -n T1

This command creates the output file: "replicase\_picornavirus\_protein\_mito.tree"

Section 4
---------

Then root two phylogenetic trees derived in section 3 using midpoint rooting and outgroup rooting. The outgroup is mitovirus.

### 1. Root at Outgroup

#### picornavirus

``` r
rm(list=ls())
library(ape) 
tree_replicase_picorna_triticum_nucleotide_mito=read.tree("replicase_picornavirus_triticum_nucleotide_mitovirus_outgroup_noroot.tree")
tree_replicase_picorna_triticum_nucleotide_mito=unroot(tree_replicase_picorna_triticum_nucleotide_mito)
tree.1=root(tree_replicase_picorna_triticum_nucleotide_mito, outgroup=c("NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome"  ,"NC004052.1-205-2556Ophiostomamitovirus4completegenome"))
#write.tree(tree.1,file="replicase_picornavirus_triticum_nucleotide_mitovirus_outgrouprooting.tree")

tree_replicase_potyvirus_triticum_protein_mito=read.tree("replicase_picornavirus_triticum_protein_mitovirus_outgroup_noroot.tree")
tree_replicase_potyvirus_triticum_protein_mito=unroot(tree_replicase_potyvirus_triticum_protein_mito)
tree.2=root(tree_replicase_potyvirus_triticum_protein_mito, outgroup=c("NP660174.1putativeRNA-dependentRNApolymerase_Cryphonectriaparasiticamitovirus1-NB631","NP660179.1RNA-dependentRNApolymeraseputative_Ophiostomamitovirus4"))
#write.tree(tree.2,file="replicase_picornavirus_triticum_protein_mitovirus_outgrouprooting.tree")
```

This command produces: "replicase\_picornavirus\_triticum\_nucleotide\_mitovirus\_outgrouprooting.tree" and "replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting.tree"

### 2. Root at Midpoint

#### potyvirus

``` r
library(phytools)
```

    ## Loading required package: maps

``` r
tree_replicase_picorna_triticum_nucleotide_mito_mp=midpoint.root(tree_replicase_picorna_triticum_nucleotide_mito)
#write.tree(tree_replicase_picorna_triticum_nucleotide_mito_mp,file="replicase_picornavirus_triticum_nucleotide_mitovirus_midpoint.tree")

tree_replicase_potyvirus_triticum_protein_mito_mp=midpoint.root(tree_replicase_potyvirus_triticum_protein_mito)
#write.tree(tree_replicase_potyvirus_triticum_protein_mito_mp,file="replicase_picornavirus_triticum_protein_mitovirus_midpoint.tree")
```

This command produces: "replicase\_picornavirus\_triticum\_nucleotide\_mitovirus\_midpoint.tree" and "replicase\_picornavirus\_triticum\_protein\_mitovirus\_midpoint.tree".

Section 5
---------

Before the ancestral sequence construction, delete the taxas which are not contained in the 5UTR files in the tree. We now focus on the nucleotide files.

##### picornavirus

``` r
# Midpoint
rm(list=ls())
library(ape)
library(phylotools)
picorna_triticum_replicase_match=read.fasta("replicase_picornavirus_triricum_nucleotide_5UTR_mito_combination.fasta")
seqname_picorna_triticum_replicase_match=picorna_triticum_replicase_match$seq.name
picorna_triticum_replicasetree_midpoint=read.tree("replicase_picornavirus_triticum_nucleotide_mitovirus_midpoint.tree")

picorna_triticum_replicasetree_midpoint_deletion=keep.tip(picorna_triticum_replicasetree_midpoint,seqname_picorna_triticum_replicase_match)

#write.tree(picorna_triticum_replicasetree_midpoint_deletion,"picornavirus_triticum_nucleotide_mitovirus_replicasetree_midpoint_deletion.tree")

# Outgroup
picorna_triticum_replicasetree_outgroup=read.tree("replicase_picornavirus_triticum_nucleotide_mitovirus_outgrouprooting.tree")
picorna_triticum_replicasetree_outgroup_deletion=keep.tip(picorna_triticum_replicasetree_outgroup,seqname_picorna_triticum_replicase_match)

#write.tree(picorna_triticum_replicasetree_outgroup_deletion,"picornavirus_triticum_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")
```

Two trees are produced:

Midpont: "picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree" Outgroup: "picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree"

Section 6
---------

With the purpose of executing the ancestral sequence, it is necessary to assure the name of the taxa ID of the replicase tree matches the name of the 5UTR files. Therefore, this step removes the space in the ID of virus and gap in the sequence alignment.

the columns of missing gap only are deleted

The software to delete the columns in the aligment is MEGAX.

##### picornavirus:

``` r
rm(list=ls())
library(phylotools)
library(stringr)
seq=read.fasta("picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "picornaviridae_triticum_5UTR_mitovirus_outgroup_ReplicaseID_nospace.fasta")
```

The files of the potyvirus 5UTR without space and gap:

"picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas"

Section 7
---------

#### Ancestral sequence reconstruction

Next, use the RAXML to construct the ancestral sequence.

The sequence alignment file is the 5 UTR files, and the phylogenetic tree is constructed from the replicase file. In order to infer the ancestral sequences successfully, the 5UTR files and replicase files have the exact same taxas. Their ID(the name of the taxon) also matches.

Note: In order to make the output files more clear and display it correctly, I change S1 and S2 to .fasta in the final file names to differential sequence alignment files in order to show a uniform name to show the ancestor sequence. The first part of the name of the files remains same.

##### picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree -s picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S1

This command product output file: "picornaviridae\_triticum\_5UTR\_mitovirus\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t picornavirus\_triticum\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree -s picornaviridae\_triticum\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S2

This command product output file: "picornaviridae\_triticum\_5UTR\_mitovirus\_replicasetree\_outgroup\_AncestralSeq.fasta"

Section 8
---------

This part focuses on the ancestral construction of the protein file. Since the protein file does not contain 5UTR files, we use the replicase tree and replicase sequence alignment for ancestral sequence construction.

Since the protein files contain X, which cannot be discerned by RAXML, so I use - instead.

The new files are:

"replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgroup\_nospaceX.fas"

Note: In order to make the output files more clear and display it correctly, I change P1 and P2 to .fasta in the final file names to differential sequence alignment files in order to show a uniform name to show the ancestor sequence. The first part of the name of the files remains same.

The model is JTT model.

Three unique ID proteins rf1NC040673Pelodiscussinensispicornavirus1strainCNSR2011completegenomeputative3DRNA-dependentRNApolymerase-64047856 NP660179.1RNA-dependentRNApolymeraseputativeOphiostomamitovirus4 NP660174.1putativeRNA-dependentRNApolymeraseCryphonectriaparasiticamitovirus1-NB631

##### picornavirus

Midpoint:

raxmlHPC -f A -p 23456 -t replicase\_picornavirus\_triticum\_protein\_mitovirus\_midpoint.tree replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgroup\_nospaceX.fas -m PROTGAMMAJTT -n P1

This command product output file: "replicase\_picornavirus\_triticum\_protein\_mitovirus\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 23456 -t replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting.tree -s replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgroup\_nospaceX.fas -m PROTGAMMAJTT -n P2

This command product output file: "replicase\_picornavirus\_triticum\_protein\_mitovirus\_outgrouprooting\_AncestralSeq.fasta"
