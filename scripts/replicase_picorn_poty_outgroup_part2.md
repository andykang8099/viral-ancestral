replicase\_picorn\_poty\_outgroup\_part2
================

#### 5. Construct the phylogenetic tree based on the replicase file with matched taxa in the 5UTR file

The name of nucleotide replicase files with mathched taxa in 5UTR files:

"replicase\_potyvirus\_5UTR\_combination.fasta" "replicase\_picornavirus\_5UTR\_combination.fasta" "replicase\_picornavirus\_potyvirus\_5UTR\_combination.fasta"

Each tree is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF. According to some research, JTTF model is more reliable in the construction of phylogenetic tree of protein sequence.

Note: In order to make the output files more clear and display it correctly, I change Q1, Q2, Q3 to .tree in the final file names. In all the RAXML commands, I use Q1, Q2, Q3 to differential tree files in order to use a uniform name to show the tree successfully in the phylogenetic tree. The first part of the name of the files remains same. Unroot means random rooting or no rooting derived from RAXML.

### Nucleotide files

##### potyvirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_potyvirus\_5UTR\_combination.fasta -n Q1

("replicase\_poty5UTR\_comb\_nospace.fasta") is coverted one.

This command creates the output file: "replicase\_potyvirus\_5UTR\_comb\_noroot.tree"

##### picornavirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_5UTR\_combination.fasta -n Q2

("replicase\_picorna5UTR\_comb\_nospace.fasta") is coverted one.

This command creates the output file: "replicase\_picornavirus\_5UTR\_comb\_noroot.tree"

##### potyvirus + picornavirus

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 50 -s replicase\_picornavirus\_potyvirus\_5UTR\_combination.fasta -n Q3

("replicase\_picornapoty5UTR\_comb\_nospace.fasta") is coverted one.

This command creates the output file: "replicase\_picornavirus\_potyvirus\_5UTR\_comb\_ noroot.tree"

### 6. Root the tree by midpoint and outgroup

The aim is to compare the ancestral sequence by different way of rooting

Note: noroot means some random rooting from the RAXML.

##### potyvirus

``` r
library(phytools)
```

    ## Loading required package: ape

    ## Loading required package: maps

``` r
library(ape)
rep_poty_5UTR=read.tree("./file_for_raxml/tree/replicase_potyvirus_5UTR_comb_noroot.tree")
rep_poty_5UTR_unroot=unroot(rep_poty_5UTR)
#write.tree(rep_poty_5UTR_unroot,"replicase_potyvirus_5UTR_comb_unrooted.tree")

# Midpoint rooting

rep_poty_5UTR_midpoint=midpoint.root(rep_poty_5UTR_unroot)
#write.tree(rep_poty_5UTR_midpoint,"replicase_potyvirus_5UTR_comb_midpoint.tree")

# This part generate the file: "replicase_potyvirus_5UTR_comb_midpoint.tree"

# Outgroup Rooting
rep_poty_5UTR_outgroup=root(rep_poty_5UTR_unroot,outgroup=c("NC_004750.1.142-11551155-2744Barleyyellowdwarfvirus-PAVcompletegenome", "NC_001747.1.308-17741774-3495Potatoleafrollviruscompletegenome"))

#write.tree(rep_poty_5UTR_outgroup,"replicase_potyvirus_5UTR_comb_outgroup.tree")

# This part generate the file: "replicase_potyvirus_5UTR_comb_outgroup.tree"
```

This way generate the file:

Midpoint Rooting: "replicase\_potyvirus\_5UTR\_comb\_midpoint.tree" Root at outgroup: "replicase\_potyvirus\_5UTR\_comb\_outgroup.tree"

##### picornavirus

``` r
rep_picorna_5UTR=read.tree("./file_for_raxml/tree/replicase_picornavirus_5UTR_comb_noroot.tree")
rep_picorna_5UTR_unroot=unroot(rep_picorna_5UTR)
#write.tree(rep_picorna_5UTR_unroot,"replicase_picornavirus_5UTR_comb_unrooted.tree")

# Midpoint rooting

rep_picorna_5UTR_midpoint=midpoint.root(rep_picorna_5UTR_unroot)
#write.tree(rep_picorna_5UTR_midpoint,"replicase_picornavirus_5UTR_comb_midpoint.tree")

# This part generate the file: "replicase_picornavirus_5UTR_comb_midpoint.tree"
```

This way generate the file:

Midpoint Rooting: "replicase\_picornavirus\_5UTR\_comb\_midpoint.tree"

##### potyvirus + picornavirus

``` r
rep_picorna_poty_5UTR=read.tree("./file_for_raxml/tree/replicase_picornavirus_potyvirus_5UTR_comb_ noroot.tree")
rep_picorna_poty_5UTR_unroot=unroot(rep_picorna_poty_5UTR)
#write.tree(rep_picorna_poty_5UTR_unroot,"replicase_picornavirus_potyvirus_5UTR_comb_unrooted.tree")

# Midpoint rooting

rep_picorna_poty_5UTR_midpoint=midpoint.root(rep_picorna_poty_5UTR_unroot)
#write.tree(rep_picorna_poty_5UTR_midpoint,"replicase_picornavirus_potyvirus_5UTR_comb_midpoint.tree")

# This part generate the file: "replicase_picornavirus_potyvirus_5UTR_comb_midpoint.tree"


# Outgroup Rooting
#rep_picorna_poty_5UTR_outgroup=root(rep_picorna_poty_5UTR_unroot,outgroup=c("NC_004750.1.142-11551155-2744Barleyyellowdwarfvirus-PAVcompletegenome", "NC_001747.1.308-17741774-3495Potatoleafrollviruscompletegenome"))
#write.tree(rep_poty_5UTR_outgroup,"replicase_potyvirus_5UTR_comb_outgroup.tree")
# This part generate the file: "replicase_potyvirus_5UTR_comb_outgroup.tree"
```

This way generate the file:

Midpoint Rooting: "replicase\_picornavirus\_potyvirus\_5UTR\_comb\_midpoint.tree"

7. Delete the columns of the sequence alignment with missing gap only
---------------------------------------------------------------------

In order to construct the ancestral sequences successfully in RAXML, the column of missing gap should be deleted

The software to delete the columns in the aligment is MEGAX.

Now, the files of replicase nucleotide file with the removal of missing gap columns are below. For all these files, the columns with gaps only are deleted.

##### potyvirus:

"potyviridae\_5UTR\_replicaseID\_nospacegap.fas"

##### picornavirus:

"picornaviridae\_5UTR\_replicaseID\_nospacegap.fas"

##### potyvirus + picornavirus:

"picornaviridae\_potyviridae\_5UTR\_replicaseID\_nospacegap.fas"

**These files are prepared for later ancestral reconstruction.**

8. Ancestral sequence reconstruction
====================================

### Next, use the RAXML to construct the ancestral sequence.

The sequence alignment file is the 5 UTR files, and the phylogenetic tree is constructed from the replicase file. In order to infer the ancestral sequences successfully, the 5UTR files and replicase files have the exact same taxas. Their ID(the name of the taxon) also matches.

Note: In order to make the output files more clear and display it correctly, I change S1, S2, S3, S4 to .fasta in the final file names to differential sequence alignment files in order to show a uniform name to show the ancestor sequence. The first part of the name of the files remains same.

##### potyvirus

Midpoint:

raxmlHPC -f A -p 12345 -t replicase\_potyvirus\_5UTR\_comb\_midpoint.tree -s potyviridae\_5UTR\_replicaseID\_nospacegap.fas -m GTRGAMMA -n S1

This command product output file: "potyviridae\_5UTR\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

raxmlHPC -f A -p 12345 -t replicase\_potyvirus\_5UTR\_comb\_outgroup.tree -s potyviridae\_5UTR\_replicaseID\_nospacegap.fas -m GTRGAMMA -n S2

This command product output file: "potyviridae\_5UTR\_replicasetree\_outgroup\_AncestralSeq.fasta"

##### picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t replicase\_picornavirus\_5UTR\_comb\_midpoint.tree -s picornaviridae\_5UTR\_replicaseID\_nospacegap.fas -m GTRGAMMA -n S3

This command product output file: "picornaviridae\_5UTR\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:

##### potyvirus + picornavirus

Midpoint:

raxmlHPC -f A -p 12345 -t replicase\_picornavirus\_potyvirus\_5UTR\_comb\_midpoint.tree -s picornaviridae\_potyviridae\_5UTR\_replicaseID\_nospacegap.fas -m GTRGAMMA -n S4

This command product output file: "picornaviridae\_potyviridae\_5UTR\_replicasetree\_midpoint\_AncestralSeq.fasta"

Outgroup:
