viral ancestral reconstruction project
complete genomes and 5UTR files
================

# Analysis on the complete genome files and 5 UTR files. 

They are 

- potyviridae alone
- picornaviridae alone
- both potyviridae+picornaviridae

All files have outgroup.

The construction of the best tree on the virus genome is based on total 20 maximized likelihood search and 1 bootstrap

The following are the code and result in the Raxml.

# 1. Estimate RAxML trees

## 5UTR files

##### potyviridae+picornaviridae with 5UTR

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_potyviridae\_5UTR\_aligned.fasta -n T7

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) = 0.254346, 0.251361, 0.231342, 0.262952, and the transition rate are A &lt;-&gt; C: 1.142242, rate A &lt;-&gt; G: 2.531875, rate A &lt;-&gt; T: 1.412656, rate C &lt;-&gt; G: 0.963031, rate C &lt;-&gt; T: 3.510598, rate G &lt;-&gt; T: 1.000000

This command creates the output file: `RAxML_bestTree_picornaviridae_potyviridae_5UTR_aligned_unrooted.tree.T7`

##### potyviridae with 5UTR

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s potyviridae\_5UTR\_with\_outgroup\_aligned.fasta -n T8

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) = 0.344122, 0.237556, 0.154793, 0.263528, and the transition rate are rate A &lt;-&gt; C: 1.302173, rate A &lt;-&gt; G: 3.270232, rate A &lt;-&gt; T: 1.338095, rate C &lt;-&gt; G: 0.456309, rate C &lt;-&gt; T: 3.213243, rate G &lt;-&gt; T: 1.000000.

##### picornaviridae with 5UTR

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_5UTR\_with\_outgroup\_aligned.fasta -n T9

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.229238, 0.256425, 0.243877,and 0.270460, and the transition rate are A &lt;-&gt; C: 1.146259, rate A &lt;-&gt; G: 3.233354, rate A &lt;-&gt; T: 1.514142, rate C &lt;-&gt; G: 0.974224, rate C &lt;-&gt; T: 3.528592, rate G &lt;-&gt; T: 1.000000.

## Complete genomes files

##### potyviridae with complete genome and outgroup

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s potyviridae\_genome\_with\_outgroup\_aligned.fasta -n T10

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.313372, 0.196132, 0.233898, and 0.256599, and the transition rate are A &lt;-&gt; C: 1.302173, rate A &lt;-&gt; G: 3.270232, rate A &lt;-&gt; T: 1.338095, rate C &lt;-&gt; G: 0.456309, rate C &lt;-&gt; T: 3.213243, rate G &lt;-&gt; T: 1.000000.

##### picornaviridae with compelete genome and outgroup

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_genome\_with\_outgroup\_aligned.fasta -n T11

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.270188, 0.222885, 0.225654, and 0.281272, and the transition rate are A &lt;-&gt; C: 1.479092, rate A &lt;-&gt; G: 3.273508, rate A &lt;-&gt; T: 1.354061, rate C &lt;-&gt; G: 1.291275, rate C &lt;-&gt; T: 3.307613, rate G &lt;-&gt; T: 1.000000.

##### picornaviridae + potyviridae with compelete genome and outgroup

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_potyviridae\_genome\_aligned.fasta -n T12

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.282127, 0.217977, 0.232850, and 0.267047, and the transition rate are rate A &lt;-&gt; C: 1.689051, rate A &lt;-&gt; G: 3.261407, rate A &lt;-&gt; T: 1.388764, rate C &lt;-&gt; G: 1.257696, rate C &lt;-&gt; T: 3.595042, rate G &lt;-&gt; T: 1.000000.

# 2. Root the trees with the outgroups

Next, root the unrooted trees for later ancestral construction with midpoint rooting first. It locates the midpoint of the longest path between any two tips and putting the root in that location. Then write the rooted trees out. This method requires you to make the assumption that all of your sequences are evolving at the same rate - you should do so cautiously because this assumption does not hold for many biological datasets. In this case, the root is positioned at the midpoint between the two longest branches. If you have taxa that were not sampled at the same time point then a slight modification of this method would be required to take into account the time elapsed between samples.

## 5UTR files:
``` r
library(phytools)
tree_string_picornaviridae_potyviridae_5UTR_aligned=read.tree("RAxML_bestTree_picornaviridae_potyviridae_5UTR_aligned_unrooted.T7")

tree_string_picornaviridae_potyviridae_5UTR_aligned_rooted=midpoint.root(tree_string_picornaviridae_potyviridae_5UTR_aligned)

#write.tree(tree_string_picornaviridae_potyviridae_5UTR_aligned_rooted,file="RAxML_bestTree_picornaviridae_potyviridae_5UTR_aligned_midpoint.T7")

tree_string_potyviridae_5UTR_with_outgroup_aligned=read.tree("RAxML_bestTree_potyviridae_5UTR_with_outgroup_aligned_unrooted.T8.T2")

tree_string_potyviridae_5UTR_with_outgroup_aligned_rooted=midpoint.root(tree_string_potyviridae_5UTR_with_outgroup_aligned)

#write.tree(tree_string_potyviridae_5UTR_with_outgroup_aligned_rooted,file="RAxML_bestTree_potyviridae_5UTR_with_outgroup_aligned_midpoint.T8")

tree_string_picornaviridae_5UTR_with_outgroup_aligned=read.tree("RAxML_bestTree_picornaviridae_5UTR_with_outgroup_aligned_unrooted.T9.T3")

tree_string_picornaviridae_5UTR_with_outgroup_aligned_rooted=midpoint.root(tree_string_picornaviridae_5UTR_with_outgroup_aligned)

#write.tree(tree_string_picornaviridae_5UTR_with_outgroup_aligned_rooted,file="RAxML_bestTree_picornaviridae_5UTR_with_outgroup_aligned_midpoint.T9")
```

## Complete genomes files:

``` r
tree_string_potyviridae_genome_with_outgroup_aligned=read.tree("RAxML_bestTree_potyviridae_genome_with_outgroup_aligned_unrooted.T10.T4")

tree_string_potyviridae_genome_with_outgroup_aligned_rooted=midpoint.root(tree_string_potyviridae_genome_with_outgroup_aligned)

#write.tree(tree_string_potyviridae_genome_with_outgroup_aligned_rooted,file="RAxML_bestTree_potyviridae_genome_with_outgroup_aligned_midpoint.T10")

tree_string_picornaviridae_genome_with_outgroup_aligned=read.tree("RAxML_bestTree_picornaviridae_genome_with_outgroup_aligned_unrooted.T11.T5")

tree_string_picornaviridae_genome_with_outgroup_aligned_rooted=midpoint.root(tree_string_picornaviridae_genome_with_outgroup_aligned)

#write.tree(tree_string_picornaviridae_genome_with_outgroup_aligned_rooted,file="RAxML_bestTree_picornaviridae_genome_with_outgroup_aligned_midpoint.T11")


tree_string_picornaviridae_potyviridae_genome_aligned=read.tree("RAxML_bestTree_picornaviridae_potyviridae_genome_aligned_unrooted.T12")

tree_string_picornaviridae_potyviridae_genome_aligned_rooted=midpoint.root(tree_string_picornaviridae_potyviridae_genome_aligned)

#write.tree(tree_string_picornaviridae_potyviridae_genome_aligned_rooted,file="RAxML_bestTree_picornaviridae_potyviridae_genome_aligned_midpoint.T12")
```

Then, we want to compare the result of mid-point rooting with rooting with the outgroup to investigate the differences between them. The best possible outgroups are those available which are most closely related to our sequences of interest. If outgroups are too distantly related then they can be unreliable as they may be difficult to align reliably or have become saturated with substitutions.

For the `potyviridae\_genome\_with\_outgroup\_aligned.fasta`, Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. The IDs of these two virus are "NC\_001747.1","NC\_004750.1". For the `potyviridae\_5UTR\_with\_outgroup\_aligned.fasta`, Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. The IDs of these two virus are "NC\_001747.1.1-177","NC\_004750.1.1-144".

For the `picornaviridae\_genome\_with\_outgroup\_aligned.fasta`, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. Their IDs are NC\_002551.1 and NC\_008311.1. For the `picornaviridae\_5UTR\_with\_outgroup\_aligned.fasta`, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. The IDs are NC\_002551.1.1-22 and NC\_008311.1.1-8.

For the `picornaviridae\_potyviridae\_genome\_aligned.fasta`, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. The IDs are NC\_002551.1 and NC\_008311.1. For the `potyviridae\_picornaviridae\_5UTR.fasta`, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. The IDs are NC\_002551.1.1-22 and NC\_008311.1.1-8

**Problem:** The result shows that the specified outgroup is not monophyletic, so we need to select one outgroup.

``` r
tree_potyviridae_genome_with_outgroup_aligned=read.tree("RAxML_bestTree_potyviridae_genome_with_outgroup_aligned_unrooted.T10.T4")
#tree1=root(tree_potyviridae_genome_with_outgroup_aligned, outgroup=c("NC_001747.1","NC_004750.1"))
tree_potyviridae_5UTR_with_outgroup_aligned=read.tree("RAxML_bestTree_potyviridae_5UTR_with_outgroup_aligned_unrooted.T8.T2")
#tree2=root(tree_potyviridae_5UTR_with_outgroup_aligned, outgroup=c("NC_001747.1.1-177","NC_004750.1.1-144"))
tree_picornaviridae_genome_with_outgroup_aligned=read.tree("RAxML_bestTree_picornaviridae_genome_with_outgroup_aligned_unrooted.T11.T5")
#tree3=root(tree_picornaviridae_genome_with_outgroup_aligned, outgroup=c("NC_002551.1","NC_008311.1"))
tree_picornaviridae_5UTR_with_outgroup_aligned=read.tree("RAxML_bestTree_picornaviridae_5UTR_with_outgroup_aligned_unrooted.T9.T3")
#tree4=root(tree_picornaviridae_5UTR_with_outgroup_aligned, outgroup=c("NC_002551.1.1-22","NC_008311.1.1-8"))
tree_picornaviridae_potyviridae_genome_aligned=read.tree("RAxML_bestTree_picornaviridae_potyviridae_genome_aligned_unrooted.T12")
#tree5=root(tree_picornaviridae_potyviridae_genome_aligned, outgroup=c("NC_002551.1","NC_008311.1"))
tree_picornaviridae_potyviridae_5UTR_aligned=read.tree("RAxML_bestTree_picornaviridae_potyviridae_5UTR_aligned_unrooted.T7")
#tree6=root(tree_picornaviridae_potyviridae_5UTR_aligned, outgroup=c("NC_002551.1.1-22","NC_008311.1.1-8"))
```

# 3. Ancestral sequence reconstruction

Next, use the RAXML to construct the ancestral sequence.

First, construct the ancestral sequence using mid-point root. With the GRASP function, it infers ancestral proteins from homologous protein sequences--a process known as ancestral sequence reconstruction. GRASP efficiently determines ancestral character states, and the most supported insertions and deletions. GRASP presents all this as partial-order graphs using a visual interface, connected to a phylogenetic tree, with which the user can interact. Detailed description see <http://grasp.scmb.uq.edu.au/guide>.

GRASP can perform reconstruction based on a variety of evolutionary models including JTT, Dayoff, LG, and WAG. JTT (Jones-Taylor-Thornton) is commonly used for the analysis of amino acid sequences.

Following code shows how to use the GRASP function with RAXML:

./raxmlHPC-PTHREADS -m PROTGAMMAJTT -p 23456 -n GRASPTutorial.nwk -s GRASPTutorial.aln.

PROTGAMMAJTT specifies to optimise substitution rates. it use a GAMMA model of rate heterogeneity and use the JTT amino acid substitution matrix. -p 23456 is the random number seed for inference. -n GRASPTutorial.nwk represents the name of the output file, and -s GRASPTutorial.aln is the name of the input alignment file.
