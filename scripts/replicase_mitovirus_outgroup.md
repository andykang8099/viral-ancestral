replicase\_mitovirus\_outgroup
================

This part we focus on the second part of the replicase file, the difference between the second and first sets of the replicase file is the outgroup
---------------------------------------------------------------------------------------------------------------------------------------------------

Section 1
---------

This parts show how to remove the space and , in the ID of each virus because RAXML cannot discern space and comma when RAXML try to build the tree.

##### potyvirus

``` r
library(phylotools)
```

    ## Loading required package: ape

``` r
library(stringr)
seq=read.fasta("replicase_potyvirus_nucleotide_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_poty_nucleotide_mito_nospace.fasta")

# Convert the name of the protein file
seq1=read.fasta("replicase_potyvirus_protein_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}
seq1[119,1]="NC004046.1putativeRNA-dependentRNApolymerase-Cryphonectriaparasiticamitovirus1-NB631" 
seq1[120,1]="NC004052.1RNA-dependentRNApolymeraseputative-Ophiostomamitovirus4"

#dat2fasta(seq1, outfile = "replicase_poty_protein_mito_nospace.fasta")
```

Two files are generated: "replicase\_poty\_nucleotide\_mito\_nospace.fasta" "replicase\_poty\_protein\_mito\_nospace.fasta"

##### picornavirus

``` r
seq=read.fasta("replicase_picornavirus_nucleotide_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_picorna_nucleotide_mito_nospace.fasta")

# Convert the name of the protein file
seq1=read.fasta("replicase_picornavirus_protein_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}
seq1[85,1]="NC004046.1putativeRNA-dependentRNApolymerase-Cryphonectriaparasiticamitovirus1-NB631" 
seq1[86,1]="NC004052.1RNA-dependentRNApolymeraseputative-Ophiostomamitovirus4"

#dat2fasta(seq1, outfile = "replicase_picorna_protein_mito_nospace.fasta")
```

Two files are generated: "replicase\_picrona\_nucleotide\_mito\_nospace.fasta" "replicase\_picorna\_protein\_mito\_nospace.fasta"

##### potyvirus + picornavirus

``` r
seq=read.fasta("replicase_picornavirus_and_potyvirus_nucleotide_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_picornapoty_nucleotide_mito_nospace.fasta")

# Convert the name of the protein file
seq1=read.fasta("replicase_picornavirus_and_potyvirus_protein_mitovirus_outgroup_aligned.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}
seq1[203,1]="NC004046.1putativeRNA-dependentRNApolymerase-Cryphonectriaparasiticamitovirus1-NB631" 
seq1[204,1]="NC004052.1RNA-dependentRNApolymeraseputative-Ophiostomamitovirus4"

#dat2fasta(seq1, outfile = "replicase_picornapoty_protein_mito_nospace.fasta")
```

Two files are generated: "replicase\_picornapoty\_nucleotide\_mito\_nospace.fasta" "replicase\_picornapoty\_protein\_mito\_nospace.fasta"

Section 2
---------

This part matches the ID between the replicase files and the 5UTR files exactly. The reason is to make sure that the ancestral sequences could be constructed successfully.

The difference between the ID of the repicase file and the 5UTR file is that the replicase file contains the location, and the ID is more specified.

##### potyvirus

``` r
rm(list=ls())
library(phylotools)
library(stringr)
# 154 virus
seq_poty_5UTR_mito=read.fasta("potyviridae_5UTR_mitovirus_outgroup_aligned.fasta")
# 120 virus
seq_replicase_poty_nucleotide_mito=read.fasta("replicase_potyvirus_nucleotide_mitovirus_outgroup_aligned.fasta")
replicase_poty_nucleotide_mito_name=seq_replicase_poty_nucleotide_mito$seq.name

seq=seq_poty_5UTR_mito
name=replicase_poty_nucleotide_mito_name

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

# now, the sequene has 114 virus, less than 120

# Detect and display the extra virus in replicase but not in 5UTR

seq1=seq_replicase_poty_nucleotide_mito
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

#dat2fasta(seq1, outfile = "replicase_potyvirus_5UTR_mito_comb.fasta")
#dat2fasta(seq, outfile = "potyviridae_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
```

This command generates: "potyviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID.fasta"

##### picornavirus

``` r
rm(list=ls())
library(phylotools)
library(stringr)
# 131 virus
seq_picorna_5UTR_mito=read.fasta("picornaviridae_5UTR_mitovirus_outgroup_aligned.fasta")
# 86 virus
seq_replicase_picorna_nucleotide_mito=read.fasta("replicase_picornavirus_nucleotide_mitovirus_outgroup_aligned.fasta")
replicase_picorna_nucleotide_mito_name=seq_replicase_picorna_nucleotide_mito$seq.name

seq=seq_picorna_5UTR_mito
name=replicase_picorna_nucleotide_mito_name

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

# now, the sequene has 72 virus, less than 86

# Detect and display the extra virus in replicase but not in 5UTR

seq1=seq_replicase_picorna_nucleotide_mito
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


#dat2fasta(seq1, outfile = "replicase_picornavirus_5UTR_mito_comb.fasta")
#dat2fasta(seq, outfile = "picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
```

This command generates: "picornaviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID.fasta"

##### picornavirus and potyvirus

``` r
rm(list=ls())
library(phylotools)
library(stringr)
# 283 virus
seq_picornapoty_5UTR_mito=read.fasta("potyviridae_picornaviridae_5UTR_mitovirus_outgroup_aligned.fasta")
# 204 virus
seq_replicase_picornapoty_nucleotide_mito=read.fasta("replicase_picornavirus_and_potyvirus_nucleotide_mitovirus_outgroup_aligned.fasta")
replicase_picornapoty_nucleotide_mito_name=seq_replicase_picornapoty_nucleotide_mito$seq.name

seq=seq_picornapoty_5UTR_mito
name=replicase_picornapoty_nucleotide_mito_name

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

# now, the sequene has 184 virus, less than 204

# Detect and display the extra virus in replicase but not in 5UTR

seq1=seq_replicase_picornapoty_nucleotide_mito
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


#dat2fasta(seq1, outfile = "replicase_picornapotyvirus_5UTR_mito_comb.fasta")
#dat2fasta(seq, outfile = "potyviridae_picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
```

This command generates: "potyviridae\_picornaviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID.fasta"

Section 3
---------

Now, construct the phylogenetic tree based on the aligned sequence alignment files derived in the first section

Each tree file is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF. According to some research, JTTF model is more reliable in the construction of phylogenetic tree of protein sequence.

Note: In order to make the output files more clear and display it correctly, I change M1, M2, M3, M4 to .tree in the final file names.

In all the RAXML commands, I use M1, M2, ..., M6 to differential six tree files in order to use a uniform name to show the tree successfully in the phylogenetic tree software.The first part of the name of the files remains same (output name in the RAXML).

### Nucleotide files

##### picornavirus and potyvirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picornapoty\_nucleotide\_mito\_nospace.fasta -n M1

This command creates the output file: "replicase\_picornavirus\_and\_potyvirus\_nucleotide\_mito.tree"

##### picornavirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picorna\_nucleotide\_mito\_nospace.fasta -n M2

This command creates the output file: "replicase\_picornavirus\_nucleotide\_mito.tree"

##### potyvirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_poty\_nucleotide\_mito\_nospace.fasta -n M3

This command creates the output file: "replicase\_potyvirus\_nucleotide\_mito.tree"

### Protein files

##### picornavirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_picorna\_protein\_mito\_nospace.fasta -n M4

This command creates the output file: "replicase\_picornavirus\_protein\_mito.tree"

##### potyvirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_poty\_protein\_mito\_nospace.fasta -n M5

This command creates the output file: "replicase\_potyvirus\_protein\_mito.tree"

##### picornavirus and potyvirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_picornapoty\_protein\_mito\_nospace.fasta -n M6

This command creates the output file: "replicase\_picornavirus\_and\_potyvirus\_protein\_mito.tree"

Section 4
---------

Then root six phylogenetic tree derived in section 3 using midpoint rooting and outgroup rooting. The outgroup is mitovirus.

### 1. Root at Outgroup

#### potyvirus

``` r
rm(list=ls())
library(ape) 
tree_replicase_potyvirus_nucleotide_mito=read.tree("replicase_potyvirus_nucleotide_mito.tree")
tree_replicase_potyvirus_nucleotide_mito=unroot(tree_replicase_potyvirus_nucleotide_mito)
tree.1=root(tree_replicase_potyvirus_nucleotide_mito, outgroup=c("NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome","NC004052.1-205-2556Ophiostomamitovirus4completegenome"))
#write.tree(tree.1,file="replicase_potyvirus_nucleotide_mito_outgroup.tree")

tree_replicase_potyvirus_protein_mito=read.tree("replicase_potyvirus_protein_mito.tree")
tree_replicase_potyvirus_protein_mito=unroot(tree_replicase_potyvirus_protein_mito)
tree.2=root(tree_replicase_potyvirus_protein_mito, outgroup=c("NC004046.1putativeRNA-dependentRNApolymerase-Cryphonectriaparasiticamitovirus1-NB631","NC004052.1RNA-dependentRNApolymeraseputative-Ophiostomamitovirus4"))
#write.tree(tree.2,file="replicase_potyvirus_protein_mito_outgroup.tree")
```

This command produces: "replicase\_potyvirus\_nucleotide\_mito\_outgroup.tree" and "replicase\_potyvirus\_protein\_mito\_outgroup.tree"

#### picornavirus

``` r
library(ape)
tree_replicase_picornavirus_nucleotide_mito=read.tree("replicase_picornavirus_nucleotide_mito1.tree")
tree_replicase_picornavirus_nucleotide_mito=unroot(tree_replicase_picornavirus_nucleotide_mito)
tree.3=root(tree_replicase_picornavirus_nucleotide_mito, outgroup=c("NC004052.1-205-2556Ophiostomamitovirus4completegenome", "NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome"))
#write.tree(tree.3,file="replicase_picornavirus_nucleotide_mito_outgrouprooting.tree")

tree_replicase_picornavirus_protein_mito=read.tree("replicase_picornavirus_protein_mito.tree")
tree_replicase_picornavirus_protein_mito=unroot(tree_replicase_picornavirus_protein_mito)
tree.4=root(tree_replicase_picornavirus_protein_mito, outgroup=c("NC004046.1putativeRNA-dependentRNApolymerase-Cryphonectriaparasiticamitovirus1-NB631","NC004052.1RNA-dependentRNApolymeraseputative-Ophiostomamitovirus4"))
#write.tree(tree.4,file="replicase_picornavirus_protein_mito_outgroup.tree")
```

This command produces: "replicase\_picornavirus\_protein\_mito\_outgroup.tree" and "replicase\_picornavirus\_nucleotide\_mito\_outgrouprooting.tree"

#### potyvirus + picornavirus

``` r
library(ape) 
tree_replicase_picornapotyvirus_nucleotide_mito=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_mito.tree")
tree_replicase_picornapotyvirus_nucleotide_mito=unroot(tree_replicase_picornapotyvirus_nucleotide_mito)
tree.5=root(tree_replicase_picornapotyvirus_nucleotide_mito, outgroup=c("NC004046.1-87-2516Cryphonectriaparasiticamitovirus1-NB631completegenome","NC004052.1-205-2556Ophiostomamitovirus4completegenome"))

#write.tree(tree.5,file="replicase_picornavirus_and_potyvirus_nucleotide_mito_outgroup.tree")

tree_replicase_picornapotyvirus_protein_mito=read.tree("replicase_picornavirus_and_potyvirus_protein_mito.tree")
tree_replicase_picornapotyvirus_protein_mito=unroot(tree_replicase_picornapotyvirus_protein_mito)
tree.6=root(tree_replicase_picornapotyvirus_protein_mito, outgroup=c("NC004046.1putativeRNA-dependentRNApolymerase-Cryphonectriaparasiticamitovirus1-NB631","NC004052.1RNA-dependentRNApolymeraseputative-Ophiostomamitovirus4"))

#write.tree(tree.6,file="replicase_picornavirus_and_potyvirus_protein_mito_outgroup.tree")
```

This command produces: "replicase\_potyvirus\_nucleotide\_mito\_outgroup.tree" and "replicase\_potyvirus\_protein\_mito\_outgroup.tree"

### 2. Root at Midpoint

#### potyvirus

``` r
library(phytools)
```

    ## Loading required package: maps

``` r
tree_replicase_potyvirus_nucleotide_mito_mp=midpoint.root(tree_replicase_potyvirus_nucleotide_mito)
#write.tree(tree_replicase_potyvirus_nucleotide_mito_mp,file="replicase_potyvirus_nucleotide_mito_midpoint.tree")

tree_replicase_potyvirus_protein_mito_mp=midpoint.root(tree_replicase_potyvirus_protein_mito)
#write.tree(tree_replicase_potyvirus_protein_mito_mp,file="replicase_potyvirus_protein_mito_midpoint.tree")
```

This command produces: "replicase\_potyvirus\_nucleotide\_mito\_midpoint.tree" and "replicase\_potyvirus\_protein\_mito\_midpoint.tree"

#### picornavirus

``` r
library(phytools)
tree_replicase_picornavirus_nucleotide_mito_mp=midpoint.root(tree_replicase_picornavirus_nucleotide_mito)
#write.tree(tree_replicase_picornavirus_nucleotide_mito_mp,file="replicase_picornavirus_nucleotide_mito_midpoint.tree")

tree_replicase_picornavirus_protein_mito_mp=midpoint.root(tree_replicase_picornavirus_protein_mito)
#write.tree(tree_replicase_picornavirus_protein_mito_mp,file="replicase_picornavirus_protein_mito_midpoint.tree")
```

This command produces: "replicase\_picornavirus\_nucleotide\_mito\_midpoint.tree" and "replicase\_picornavirus\_protein\_mito\_midpoint.tree"

#### picornavirus + potyvirus

``` r
#non coding and coding region
library(phytools)
tree_replicase_picorna_poty_nucleotide_mito_mp=midpoint.root(tree_replicase_picornapotyvirus_nucleotide_mito)
#write.tree(tree_replicase_picorna_poty_nucleotide_mito_mp,file="replicase_picornavirus_and_potyvirus_nucleotide_mito_midpoint.tree")

tree_replicase_picorna_poty_protein_mito_mp=midpoint.root(tree_replicase_picornapotyvirus_protein_mito)
#write.tree(tree_replicase_picorna_poty_protein_mito_mp,file="replicase_picornavirus_and_potyvirus_protein_mito_midpoint.tree")
```

This command produces: "replicase\_picornavirus\_and\_potyvirus\_nucleotide\_mito\_midpoint.tree" and "replicase\_picornavirus\_and\_potyvirus\_protein\_mito\_midpoint.tree"

Section 5
---------

Before the ancestral sequence construction, delete the taxas which are not contained in the 5UTR files in the tree. We now focus on the nucleotide files.

##### potyvirus

``` r
# Midpoint
rm(list=ls())
library(ape)
library(phylotools)
poty_replicase_match=read.fasta("replicase_potyvirus_5UTR_mito_comb.fasta")
seqname_poty_replicase_match=poty_replicase_match$seq.name
poty_replicasetree_midpoint=read.tree("replicase_potyvirus_nucleotide_mito_midpoint.tree")

poty_replicasetree_midpoint_deletion=keep.tip(poty_replicasetree_midpoint,seqname_poty_replicase_match)

#write.tree(poty_replicasetree_midpoint_deletion,"potyvirus_nucleotide_mitovirus_replicasetree_midpoint_deletion.tree")

# Outgroup
poty_replicasetree_outgroup=read.tree("replicase_potyvirus_nucleotide_mito_outgroup.tree")
poty_replicasetree_outgroup_deletion=keep.tip(poty_replicasetree_outgroup,seqname_poty_replicase_match)

#write.tree(poty_replicasetree_outgroup_deletion,"potyvirus_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")
```

Two trees are produces:

Midpont: "potyvirus\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree" Outgroup: "potyvirus\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree"

##### picornavirus

``` r
# Midpoint
rm(list=ls())
library(ape)
library(phylotools)
picorna_replicase_match=read.fasta("replicase_picornavirus_5UTR_mito_comb.fasta")
seqname_picorna_replicase_match=picorna_replicase_match$seq.name
picorna_replicasetree_midpoint=read.tree("replicase_picornavirus_nucleotide_mito_midpoint.tree")

picorna_replicasetree_midpoint_deletion=keep.tip(picorna_replicasetree_midpoint,seqname_picorna_replicase_match)

#write.tree(picorna_replicasetree_midpoint_deletion,"picornavirus_nucleotide_mitovirus_replicasetree_midpoint_deletion.tree")


# Outgroup

picorna_replicasetree_outgroup=read.tree("replicase_picornavirus_nucleotide_mito_outgrouprooting.tree")
picorna_replicasetree_outgroup_deletion=keep.tip(picorna_replicasetree_outgroup,seqname_picorna_replicase_match)

#write.tree(picorna_replicasetree_outgroup_deletion,"picornavirus_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")
```

Two trees are produced:

Midpont: "picornavirus\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree" Outgroup: "picornavirus\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree"

##### picornavirus + potyvirus

``` r
# Midpoint
rm(list=ls())
library(ape)
library(phylotools)
picorna_poty_replicase_match=read.fasta("replicase_picornapotyvirus_5UTR_mito_comb.fasta")
seqname_picorna_poty_replicase_match=picorna_poty_replicase_match$seq.name
picorna_poty_replicasetree_midpoint=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_mito_midpoint.tree")

picorna_poty_replicasetree_midpoint_deletion=keep.tip(picorna_poty_replicasetree_midpoint,seqname_picorna_poty_replicase_match)

#write.tree(picorna_poty_replicasetree_midpoint_deletion,"picornavirus_and_potyvirus_nucleotide_mitovirus_replicasetree_midpoint_deletion.tree")

# Outgroup
picorna_poty_replicasetree_outgroup=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_mito_outgroup.tree")
picorna_poty_replicasetree_outgroup_deletion=keep.tip(picorna_poty_replicasetree_outgroup,seqname_picorna_poty_replicase_match)

#write.tree(picorna_poty_replicasetree_outgroup_deletion,"picornavirus_and_potyvirus_nucleotide_mitovirus_replicasetree_outgroup_deletion.tree")
```

Two trees are produced:

Midpont: "picornavirus\_and\_potyvirus\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree" Outgroup: "picornavirus\_and\_potyvirus\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree"

Section 6
---------

With the purpose of executing the ancestral sequence, it is necessary to assure the name of the taxa ID of the replicase tree matches the name of the 5UTR files. Therefore, this step removes the space in the ID of virus and gap in the sequence alignment.

the columns of missing gap only are deleted

The software to delete the columns in the aligment is MEGAX.

##### potyvirus:

``` r
rm(list=ls())
library(phylotools)
library(stringr)
seq=read.fasta("potyviridae_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "potyviridae_5UTR_mitovirus_outgroup_ReplicaseID_nospace.fasta")
```

The files of the potyvirus 5UTR without space and gap:

"potyviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas"

##### picornavirus:

``` r
rm(list=ls())
library(phylotools)
library(stringr)
seq=read.fasta("picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}

dat2fasta(seq, outfile = "picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID_nospace.fasta")
```

    ## picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID_nospace.fasta has been saved to  /Users/yuzhuokang/Desktop/replicase2

The files of the potyvirus 5UTR without space and gap:

"picornaviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas"

##### potyvirus + picornavirus

``` r
rm(list=ls())
library(phylotools)
library(stringr)
seq=read.fasta("potyviridae_picornaviridae_5UTR_mitovirus_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "potyviridae_picornavirida_5UTR_mitovirus_outgroup_ReplicaseID_nospace.fasta")
```

The files of the potyvirus 5UTR without space and gap:

"potyviridae\_picornavirida\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas"

Section 7
---------

#### Ancestral sequence reconstruction

Next, use the RAXML to construct the ancestral sequence.

The sequence alignment file is the 5 UTR files, and the phylogenetic tree is constructed from the replicase file. In order to infer the ancestral sequences successfully, the 5UTR files and replicase files have the exact same taxas. Their ID(the name of the taxon) also matches.

Note: In order to make the output files more clear and display it correctly, I change S1, S2, S3, S4, and S5 to .fasta in the final file names to differential sequence alignment files in order to show a uniform name to show the ancestor sequence. The first part of the name of the files remains same.

##### potyvirus

Midpoint:

raxmlHPC -f A -p 12345 -t potyvirus\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree -s potyviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S1

This command product output file: "potyviridae\_5UTR\_mitovirus\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t potyvirus\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree -s potyviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S2

This command product output file: "potyviridae\_5UTR\_mitovirus\_replicasetree\_outgroup\_AncestralSeq.fasta"

##### picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t picornavirus\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree -s picornairidae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S3

This command product output file: "picornaviridae\_5UTR\_mitovirus\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t picornavirus\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree -s picornairidae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S4

This command product output file: "picornaviridae\_5UTR\_mitovirus\_replicasetree\_outgrouprooting\_AncestralSeq.fasta"

##### potyvirus + picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t picornavirus\_and\_potyvirus\_nucleotide\_mitovirus\_replicasetree\_midpoint\_deletion.tree -s potyviridae\_picornaviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -n S5

This command product output file: "picornaviridae\_and\_potyviridae\_5UTR\_mitovirus\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t picornavirus\_and\_potyvirus\_nucleotide\_mitovirus\_replicasetree\_outgroup\_deletion.tree -s potyviridae\_picornaviridae\_5UTR\_mitovirus\_outgroup\_ReplicaseID\_nospacegap.fas -m GTRGAMMA -n S6

This command product output file: "picornaviridae\_and\_potyviridae\_5UTR\_mitovirus\_replicasetree\_outgroup\_AncestralSeq.fasta"

Section 8
---------

This part focuses on the ancestral construction of the protein file. Since the protein file does not contain 5UTR files, we use the replicase tree and replicase sequence alignment for ancestral sequence construction.

Since the protein files contain X, which cannot be discerned by RAXML, so I use - instead.

The new files are:

"replicase\_poty\_protein\_mito\_nospaceX.fas" "replicase\_picorna\_protein\_mito\_nospaceX.fas" "replicase\_picornapoty\_protein\_mito\_nospaceX.fas"

Note: In order to make the output files more clear and display it correctly, I change P1, P2, P3, P4, P5, and P6 to .fasta in the final file names to differential sequence alignment files in order to show a uniform name to show the ancestor sequence. The first part of the name of the files remains same.

The model is JTT model.

##### potyvirus

Midpoint:

raxmlHPC -f A -p 23456 -t replicase\_potyvirus\_protein\_mito\_midpoint.tree -s replicase\_poty\_protein\_mito\_nospaceX.fas -m PROTGAMMAJTT -n P1

This command product output file: "potyviridae\_mitovirus\_protein\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 23456 -t replicase\_potyvirus\_protein\_mito\_outgroup.tree -s replicase\_poty\_protein\_mito\_nospaceX.fas -m PROTGAMMAJTT -n P2

This command product output file: "potyviridae\_mitovirus\_protein\_outgrouproot\_AncestralSeq.fasta"

##### picornavirus

Midpoint:

raxmlHPC -f A -p 23456 -t replicase\_picornavirus\_protein\_mito\_midpoint.tree -s replicase\_picorna\_protein\_mito\_nospaceX.fass -m PROTGAMMAJTT -n P3

This command product output file: "picornaviridae\_mitovirus\_protein\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 23456 -t replicase\_picornavirus\_protein\_mito\_outgroup.tree -s replicase\_picorna\_protein\_mito\_nospaceX.fas -m PROTGAMMAJTT -n P4

This command product output file: "picornaviridae\_mitovirus\_protein\_midpoint\_AncestralSeq.fasta"
