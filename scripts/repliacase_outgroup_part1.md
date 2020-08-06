replicase\_picorn\_poty\_outgroup\_part1
================

1. Construct the tree using RAXML
=================================

Each tree file is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF. According to some research, JTTF model is more reliable in the construction of phylogenetic tree of protein sequence.

Note: In order to make the output files more clear and display it correctly, I change R1, R2, R3, R4 to .tree in the final file names.

In all the RAXML commands, I use R1, R2, ..., R6 to differential six tree files. In order to use a uniform name to show the tree successfully in the phylogenetic tree software, I change the suffix to tree for all the six tree files. The first part of the name of the files remains same.

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

This command creates the output file: "replicase\_picornavirus\_protein\_aligned.tree"

##### potyvirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_potyvirus\_protein\_aligned.fasta -n R5

This command creates the output file: "replicase\_potyvirus\_protein\_aligned.tree"

##### picornavirus and potyvirus

raxmlHPC -f a -m PROTCATJTTF -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_and\_potyvirus\_protein\_aligned.fasta -n R6

This command creates the output file: "replicase\_picornavirus\_and\_potyvirus\_protein\_aligned.tree"

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

3. Root the tree at the midpoint
================================

##### potyvirus

``` r
library(phytools)
```

    ## Loading required package: maps

``` r
library(ape) 
tree_replicase_potyvirus_nucleotide_aligned=read.tree("replicase_potyvirus_nucleotide_aligned.tree")
tree_replicase_potyvirus_nucleotide_aligned_midpoint=midpoint.root(tree_replicase_potyvirus_nucleotide_aligned)
#write.tree(tree_replicase_potyvirus_nucleotide_aligned_midpoint,file="replicase_potyvirus_nucleotide_aligned_midpoint.tree")

#Output file: "replicase_potyvirus_nucleotide_aligned_midpoint.tree"

tree_replicase_potyvirus_protein_aligned=read.tree("replicase_potyvirus_protein_aligned.tree")
tree_replicase_potyvirus_protein_aligned_midpoint=midpoint.root(tree_replicase_potyvirus_protein_aligned)

#write.tree(tree_replicase_potyvirus_protein_aligned_midpoint,file="replicase_potyvirus_protein_aligned_midpoint.tree")

#Output file: "replicase_potyvirus_protein_aligned_midpoint.tree"
```

##### picornavirus

``` r
tree_replicase_picornvirus_nucleotide_aligned=read.tree("replicase_picornavirus_nucleotide_aligned.tree")
tree_replicase_picornvirus_nucleotide_aligned_midpoint=midpoint.root(tree_replicase_picornvirus_nucleotide_aligned)
#write.tree(tree_replicase_picornvirus_nucleotide_aligned_midpoint,file="replicase_picornavirus_nucleotide_aligned_midpoint.tree")
#
#Output file: "replicase_picornavirus_nucleotide_aligned_midpoint.tree"

tree_replicase_picornvirus_protein_aligned=read.tree("replicase_picornavirus_protein_aligned.tree")
tree_replicase_picornvirus_protein_aligned_midpoint=midpoint.root(tree_replicase_picornvirus_protein_aligned)
#write.tree(tree_replicase_picornvirus_protein_aligned_midpoint,file="replicase_picornvirus_protein_aligned_midpoint.tree")
#
#Output file: "replicase_picornavirus_protein_aligned_midpoint.tree"
```

##### picornavirus and potyvirus

``` r
library(ape) 
tree_replicase_picornavirus_and_potyvirus_nucleotide_aligned=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_aligned.tree")
tree_replicase_picornavirus_and_potyvirus_nucleotide_aligned_midpoint=midpoint.root(tree_replicase_picornavirus_and_potyvirus_nucleotide_aligned) 
#write.tree(tree_replicase_picornavirus_and_potyvirus_nucleotide_aligned_midpoint,file="replicase_picornavirus_and_potyvirus_nucleotide_aligned_midpoint.tree")

tree_replicase_picornavirus_and_potyvirus_protein_aligned=read.tree("replicase_picornavirus_and_potyvirus_protein_aligned.tree")
tree_replicase_picornavirus_and_potyvirus_protein_aligned_midpoint=midpoint.root(tree_replicase_picornavirus_and_potyvirus_protein_aligned) 
#write.tree(tree_replicase_picornavirus_and_potyvirus_protein_aligned_midpoint,file="replicase_picornavirus_and_potyvirus_protein_aligned_midpoint.tree")
```

### Before the ancestral sequence, match the ID of each virus in the replicase file with the ID of each virus in 5UTR name. The difference between the ID of the repicase file and the 5UTR file is that the replicase file contains the location, and the ID is more specified.

#### potyvirus

``` r
ReadFasta<-function(file) {
  # Read the file line by line
  fasta<-readLines(file)
  # Identify header lines
  ind<-grep(">", fasta)
  # Identify the sequence lines
  s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
  # Process sequence lines
  seqs<-rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
  }
  # Create a data frame 
  DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
  # Return the data frame as a result object from the function
  return(DF)
}

library(phylotools)
# 154 virus
seq_potyviridae_5UTR_with_outgroup_aligned=read.fasta("potyviridae_5UTR_with_outgroup_aligned.fasta")
# 120 virus
seq_replicase_potyvirus_nucleotide_aligned=read.fasta("replicase_potyvirus_nucleotide_aligned.fasta")
replicase_potyvirus_nucleotide_name=seq_replicase_potyvirus_nucleotide_aligned$seq.name
# Convert the name and ID of the replicase file to csv 
#replicase_potyvirus_nucleotide_namedf=data.frame(replicase_potyvirus_nucleotide_name)
#write.csv(replicase_potyvirus_nucleotide_namedf,"replicase_potyvirus_nucleotide_name.csv")

seq=seq_potyviridae_5UTR_with_outgroup_aligned
name=replicase_potyvirus_nucleotide_name

library(stringr)
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

# Delete extra virus
index=which(seq[,1]=="not found")
seq=seq[-(index),]
# Create a new row name
rownames(seq)=1:nrow(seq)


# Delete the extra virus in replicase but not in 5UTR

seq1=seq_replicase_potyvirus_nucleotide_aligned
for ( i in 1: nrow(seq1)) {
   if (sum(seq1[i,1]==seq[,1])==0) {
     seq1[i,1]="not match 5UTR"
   }
}
index=which(seq1[,1]=="not match 5UTR")
seq1=seq1[-(index),]
# Create a new row name
rownames(seq1)=1:nrow(seq1)


#dat2fasta(seq1, outfile = "replicase_potyvirus_5UTR_combination.fasta")
#dat2fasta(seq, outfile = "potyviridae_5UTR_with_outgroup_ReplicaseID.fasta")
```

##### picornavirus

``` r
seq_picornaviridae_5UTR_with_outgroup_aligned=read.fasta("picornaviridae_5UTR_with_outgroup_aligned.fasta")
# 120 virus
seq_replicase_picornavirus_nucleotide_aligned=read.fasta("replicase_picornavirus_nucleotide_aligned.fasta")
replicase_picornavirus_nucleotide_name=seq_replicase_picornavirus_nucleotide_aligned$seq.name
# Convert the name and ID of the replicase file to csv 
#replicase_picornavirus_nucleotide_namedf=data.frame(replicase_picornavirus_nucleotide_name)
#write.csv(replicase_picornavirus_nucleotide_namedf,"replicase_picornavirus_nucleotide_name.csv")

seq=seq_picornaviridae_5UTR_with_outgroup_aligned
name=replicase_picornavirus_nucleotide_name

library(stringr)
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

# Delete extra virus in 5UTR but not in virus
index=which(seq[,1]=="not found")
seq=seq[-(index),]
# Create a new row name
rownames(seq)=1:nrow(seq)

# Delete the extra virus in replicase but not in 5UTR

seq1=seq_replicase_picornavirus_nucleotide_aligned
for ( i in 1: nrow(seq1)) {
   if (sum(seq1[i,1]==seq[,1])==0) {
     seq1[i,1]="not match 5UTR"
   }
}
index=which(seq1[,1]=="not match 5UTR")
seq1=seq1[-(index),]
# Create a new row name
rownames(seq1)=1:nrow(seq1)

#dat2fasta(seq1, outfile = "replicase_picornavirus_5UTR_combination.fasta")
#dat2fasta(seq, outfile = "picornaviridae_5UTR_with_outgroup_ReplicaseID.fasta")
```

##### picornavirus and potyvirus

``` r
seq_picorn_poty_5UTR_with_outgroup_aligned=read.fasta("picornaviridae_potyviridae_5UTR_aligned.fasta")
seq_replicase_picorn_poty_nucleotide_aligned=read.fasta("replicase_picornavirus_and_potyvirus_nucleotide_aligned.fasta")
replicase_picorn_poty_nucleotide_name=seq_replicase_picorn_poty_nucleotide_aligned$seq.name
# Convert the name and ID of the replicase file to csv 
#replicase_picorn_poty_nucleotide_namedf=data.frame(replicase_picorn_poty_nucleotide_name)
#write.csv(replicase_picorn_poty_nucleotide_namedf,"replicase_picorn_poty_nucleotide_name.csv")

seq=seq_picorn_poty_5UTR_with_outgroup_aligned
name=replicase_picorn_poty_nucleotide_name

library(stringr)
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

# Delete extra virus in 5UTR but not in virus
index=which(seq[,1]=="not found")
seq=seq[-(index),]
# Create a new row name
rownames(seq)=1:nrow(seq)

# Delete the extra virus in replicase but not in 5UTR

seq1=seq_replicase_picorn_poty_nucleotide_aligned
for ( i in 1: nrow(seq1)) {
   if (sum(seq1[i,1]==seq[,1])==0) {
     seq1[i,1]="not match 5UTR"
   }
}
index=which(seq1[,1]=="not match 5UTR")
seq1=seq1[-(index),]
# Create a new row name
rownames(seq1)=1:nrow(seq1)

#dat2fasta(seq1, outfile = "replicase_picornavirus_potyvirus_5UTR_combination.fasta")
#dat2fasta(seq, outfile = "picornaviridae_potyviridae_5UTR_with_outgroup_ReplicaseID.fasta")
```

In summary,

### The name of nucleotide replicase files with mathched taxa in 5UTR files:

"replicase\_potyvirus\_5UTR\_combination.fasta" "replicase\_picornavirus\_5UTR\_combination.fasta" "replicase\_picornavirus\_potyvirus\_5UTR\_combination.fasta"
