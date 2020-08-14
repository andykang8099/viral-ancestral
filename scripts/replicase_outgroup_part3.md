replicase\_outgroup\_part3
================

Note: This part focuses on the nucleotide
=========================================

#### Changing the name first

Before the ancestral reconstruction, this parts show how to remove the space and , in the ID of each virus because RAXML cannot discern space and comma

##### potyvirus

``` r
library(phylotools)
```

    ## Loading required package: ape

``` r
library(stringr)
seq=read.fasta("replicase_potyvirus_5UTR_combination.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_poty5UTR_comb_wospace.fasta")

# Now, convert the original 5UTR file
seq1=read.fasta("potyviridae_5UTR_with_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq1, outfile = "poty5UTR_replicaseID_wospace.fasta")
```

Two files are generated: "replicase\_poty5UTR\_comb\_wospace.fasta" "poty5UTR\_replicaseID\_wospace.fasta"

##### picornavirus

``` r
rm(list=ls())
seq=read.fasta("replicase_picornavirus_5UTR_combination.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_picorna5UTR_comb_wospace.fasta")

# Now, convert the original 5UTR file
seq1=read.fasta("picornaviridae_5UTR_with_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq1, outfile = "picorna5UTR_replicaseID_wospace.fasta")
```

Two files are generated: "replicase\_picorna5UTR\_comb\_wospace.fasta" "picorna5UTR\_replicaseID\_wospace.fasta"

##### potyvirus + picornavirus

``` r
rm(list=ls())
seq=read.fasta("replicase_picornavirus_potyvirus_5UTR_combination.fasta")
for ( i in 1:nrow(seq)) {
  name=strsplit(seq[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq, outfile = "replicase_picorna_poty_5UTR_comb_wospace.fasta")

# Now, convert the original 5UTR file
seq1=read.fasta("picornaviridae_potyviridae_5UTR_with_outgroup_ReplicaseID.fasta")
for ( i in 1:nrow(seq1)) {
  name=strsplit(seq1[i,1], split=" ")[[1]]
  name1=str_replace(name, "([(])", "-")
  name2=str_replace(name1, "([)])", "")
  name3=str_replace(name2, "([,])", "")
  name4=str_replace(name3, "([_])", "")
  name5=str_replace(name4, "([:])", "-")
  seq1[i,1]=paste(name5, collapse = "")
}
#dat2fasta(seq1, outfile = "picorna_poty_5UTR_replicaseID_wospace.fasta")
```

Two files are generated: "replicase\_picorna5UTR\_comb\_wospace.fasta" "picorna5UTR\_replicaseID\_wospace.fasta"

### 5. Before the ancestral sequence construction, delete the taxas which are not displayed in the 5UTR files in the tree.

We now focus on the nucleotide files. Before the deletion, the trees derived from the replicase files are already root with both midpoint rooting and outgroup rooting.

##### potyvirus

``` r
# Midpoint

library(ape)
library(phylotools)
poty_nucleotide_replicase_match=read.fasta("replicase_poty5UTR_comb_wospace.fasta")
seqname_poty_nucleotide_replicase_match=poty_nucleotide_replicase_match$seq.name
poty_nucleotide_midpoint_replicasetree=read.tree("replicase_potyvirus_nucleotide_aligned_midpoint.tree")

poty_nucleotide_midpoint_replicasetree_deletion=keep.tip(poty_nucleotide_midpoint_replicasetree,seqname_poty_nucleotide_replicase_match)

#write.tree(poty_nucleotide_midpoint_replicasetree_deletion,"potyvirus_nucleotide_midpoint_replicasetree_deletion.tree")

# Outgroup
poty_nucleotide_outgroup_replicasetree=read.tree("replicase_potyvirus_nucleotide_aligned_rooted.tree")

poty_nucleotide_outgroup_replicasetree_deletion=keep.tip(poty_nucleotide_outgroup_replicasetree,seqname_poty_nucleotide_replicase_match)

#write.tree(poty_nucleotide_outgroup_replicasetree_deletion,"potyvirus_nucleotide_outgroup_replicasetree_deletion.tree")
```

Two trees are produces:

Midpont: "potyvirus\_nucleotide\_midpoint\_replicasetree\_deletion.tree" Outgroup: "potyvirus\_nucleotide\_outgroup\_replicasetree\_deletion.tree"

##### picornavirus

``` r
picorna_nucleotide_replicase_match=read.fasta("replicase_picorna5UTR_comb_wospace.fasta")
seqname_picorna_nucleotide_replicase_match=picorna_nucleotide_replicase_match$seq.name
picorna_nucleotide_midpoint_replicasetree=read.tree("replicase_picornavirus_nucleotide_aligned_midpoint.tree")

picorna_nucleotide_midpoint_replicasetree_deletion=keep.tip(picorna_nucleotide_midpoint_replicasetree,seqname_picorna_nucleotide_replicase_match)

#write.tree(picorna_nucleotide_midpoint_replicasetree_deletion,"picornavirus_nucleotide_midpoint_replicasetree_deletion.tree")

# Outgroup
picorna_nucleotide_outgroup_replicasetree=read.tree("replicase_picornavirus_nucleotide_aligned_rooted.tree")

picorna_nucleotide_outgroup_replicasetree_deletion=keep.tip(picorna_nucleotide_outgroup_replicasetree,seqname_picorna_nucleotide_replicase_match)

#write.tree(picorna_nucleotide_outgroup_replicasetree_deletion,"picornavirus_nucleotide_outgroup_replicasetree_deletion.tree")
```

Two trees are produced:

Midpont: "picornavirus\_nucleotide\_midpoint\_replicasetree\_deletion.tree" Outgroup: "picornavirus\_nucleotide\_outgroup\_replicasetree\_deletion.tree"

##### potyvirus + picornavirus

``` r
#Midpoint

picornapoty_nucleotide_replicase_match=read.fasta("replicase_picorna_poty_5UTR_comb_wospace.fasta")
seqname_picornapoty_nucleotide_replicase_match=picornapoty_nucleotide_replicase_match$seq.name
picornapoty_nucleotide_midpoint_replicasetree=read.tree("replicase_picornavirus_and_potyvirus_nucleotide_aligned_midpoint.tree")

picornapoty_nucleotide_midpoint_replicasetree_deletion=keep.tip(picornapoty_nucleotide_midpoint_replicasetree,seqname_picornapoty_nucleotide_replicase_match)

write.tree(picornapoty_nucleotide_midpoint_replicasetree_deletion,"picornavirus_potyvirus_nucleotide_midpoint_replicasetree_deletion.tree")
```

One tree is produced:

Midpont: "picornavirus\_potyvirus\_nucleotide\_midpoint\_replicasetree\_deletion.tree"

6. Delete the columns of the sequence alignment with missing gap only
---------------------------------------------------------------------

In order to construct the ancestral sequences successfully in RAXML, the column of missing gap should be deleted

The software to delete the columns in the aligment is MEGAX.

Now, the files of replicase nucleotide file with the removal of missing gap columns are below. For all these files, the columns with gaps only are deleted.

##### potyvirus:

"poty5UTR\_replicaseID\_wospacegap.fas"

##### picornavirus:

"picorna5UTR\_replicaseID\_wospacegap.fas"

##### potyvirus + picornavirus:

"picorna\_poty\_5UTR\_replicaseID\_wospacegap.fas"

7. Ancestral sequence reconstruction
====================================

### Next, use the RAXML to construct the ancestral sequence.

The sequence alignment file is the 5 UTR files, and the phylogenetic tree is constructed from the replicase file. In order to infer the ancestral sequences successfully, the 5UTR files and replicase files have the exact same taxas. Their ID(the name of the taxon) also matches.

Note: In order to make the output files more clear and display it correctly, I change S1, S2, S3, S4, and S5 to .fasta in the final file names to differential sequence alignment files in order to show a uniform name to show the ancestor sequence. The first part of the name of the files remains same.

##### potyvirus

Midpoint:

raxmlHPC -f A -p 12345 -t potyvirus\_nucleotide\_midpoint\_replicasetree\_deletion.tree -s poty5UTR\_replicaseID\_wospacegap.fas -m GTRGAMMA -n S1

This command product output file: "potyviridae\_5UTR\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t potyvirus\_nucleotide\_outgroup\_replicasetree\_deletion.tree -s potyviridae\_5UTR\_replicaseID\_nospacegap.fas -m GTRGAMMA -n S2

This command product output file: "potyviridae\_5UTR\_replicasetree\_outgroup\_AncestralSeq.fasta"

##### picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t replicase\_picornavirus\_nucleotide\_aligned\_midpoint.tree -s picorna5UTR\_replicaseID\_wospacegap.fas -m GTRGAMMA -n S3

This command product output file: "picornaviridae\_5UTR\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t replicase\_picornavirus\_nucleotide\_aligned\_rooted.tree -s picorna5UTR\_replicaseID\_wospacegap.fas -m GTRGAMMA -n S4

This command product output file: "picornaviridae\_5UTR\_replicasetree\_outgroup\_AncestralSeq.fasta"

##### potyvirus + picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t picornavirus\_potyvirus\_nucleotide\_midpoint\_replicasetree\_deletion.tree -s picorna\_poty\_5UTR\_replicaseID\_wospacegap.fas -m GTRGAMMA -n S5

This command product output file: "picornaviridae\_potyviridae\_5UTR\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup: Not available, because the selected outgroup is not monophyletic.

Note: this part focuses on the ancestral construction of the protein
====================================================================
