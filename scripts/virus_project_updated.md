virus\_project
================

### First, read four virus files, picornaviridae\_aligned, potyviridae\_aligned, picornaviridae\_5UTR\_aligned, and potyviridae\_5UTR\_aligned into dataframes, which are prepared for later ancestral construction.

``` r
library(ape)
library(phyloch) 
```

    ## Loading required package: colorspace

    ## Loading required package: XML

``` r
library(phylotools)

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

# The sequences of picornaviridae, 226 virus
seq_picornaviridae=ReadFasta("picornaviridae_aligned.fasta")
nrow(seq_picornaviridae)
```

    ## [1] 226

``` r
# The sequences of potyviridae, 218 virus
seq_potyviridae=ReadFasta("potyviridae_aligned.fasta")
nrow(seq_potyviridae)
```

    ## [1] 218

``` r
# The sequences of picornaviridae_5UTR, 129 virus
seq_picornaviridae_5UTR=ReadFasta("picornaviridae_5UTR_aligned.fas")
nrow(seq_picornaviridae_5UTR)
```

    ## [1] 129

``` r
# The sequences of potyviridae_5UTR, 152 virus
seq_potyviridae_5UTR=ReadFasta("potyviridae_5UTR_aligned.fas")
nrow(seq_potyviridae_5UTR)
```

    ## [1] 152

### Second, use Raxml ( 20 ML search + bootstrapping ) to get the tree with best likelihood. This part of commands is from bash (Terminal).

Note1: the file size of picornaviridae\_aligned and potyviridae\_aligned is too large. As a result, I limit the bootstrap search time to 1. The search times for Maximal Likelihood (ML) are always 20 times. The estimated best tree could be further improved if the bootstrap times increase. The bootstrap search times for picornaviridae\_5UTR\_aligned and potyviridae\_5UTR\_aligned are 6 and 10 times.

Note2: For picornaviridae\_5UTR\_aligned and potyviridae\_5UTR\_aligned, the names of some virus sequences contain colon (:), so I tried to delete the colon when the virus name contain it in order to make the fas file executable in Raxml. The files with deleted colon are picornaviridae\_5UTR\_aligned1 and potyviridae\_5UTR\_aligned1.

In the picornaviridae\_5UTR\_aligned, names of five virus sequences contain colon (:). They are NC\_004421.1\_Aichivirus\_B\_genomic\_RNA\_complete\_genome\_strain\_U-1,NC\_033695.1\_Enterovirus\_AN12\_genomic\_RNA\_complete\_genome\_strain\_AN12,NC\_027919.1\_Kobuvirus\_cattle/Kagoshima-1-22-KoV/2014/JPN\_genomic\_RNA\_nearly\_complete\_genome\_strain\_Kagoshima-1-22-KoV/2014/JPN,NC\_027918.1\_Kobuvirus\_cattle/Kagoshima-2-24-KoV/2015/JPN\_genomic\_RNA\_complete\_genome\_strain\_Kagoshima-2-24-KoV/2015/JPN,NC\_025474.1\_Crohivirus\_gene\_for\_polyprotein\_complete\_cds\_strain\_ZM54. For potyviridae\_5UTR\_aligned, the names of NC\_014327.1\_Pepper\_yellow\_mosaic\_virus\_RNA\_complete\_genome\_strain\_DF,NC\_030391.1\_Wild\_onion\_symptomless\_virus\_gene\_for\_polyprotein\_complete\_cds\_isolate\_TUR256-1,NC\_021786.1\_Habenaria\_mosaic\_virus\_genomic\_RNA\_complete\_genome\_isolate\_Ha-1 contain colon (:).

### Bash code for rapid bootstrapping in Raxml:

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_aligned.fasta -n T1

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s potyviridae\_aligned.fasta -n T2

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\#6 -s picornaviridae\_5UTR\_aligned1.fas -n T3

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 10 -s potyviridae\_5UTR\_aligned1.fas -n T4

### Relevant parameters for each ML search + bootstrapping best tree estimation

picornaviridae\_aligned.fasta: alpha: 4.934338, Tree-Length: 126.910679, rate A &lt;-&gt; C: 1.448674, rate A &lt;-&gt; G: 3.216025, rate A &lt;-&gt; T: 1.314619, rate C &lt;-&gt; G: 1.269701, rate C &lt;-&gt; T: 3.239043, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.270869, 0.221802, 0.224860, 0.282469).

potyviridae\_aligned.fasta: alpha: 0.932693, Tree-Length: 132.828870, rate A &lt;-&gt; C: 2.486990, rate A &lt;-&gt; G: 3.739170, rate A &lt;-&gt; T: 1.369314, rate C &lt;-&gt; G: 1.513431, rate C &lt;-&gt; T: 5.417214, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.313614,0.195751,0.233819,0.256815)

picornaviridae\_5UTR\_aligned.fas: alpha: 5.094113, Tree-Length: 110.802014, rate A &lt;-&gt; C: 1.172756, rate A &lt;-&gt; G: 2.585556, rate A &lt;-&gt; T: 1.493185, rate C &lt;-&gt; G: 0.866442, rate C &lt;-&gt; T: 2.625791, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.225632,0.257601,0.242809,0.273959)

potyviridae\_5UTR\_aligned.fas: alpha: 1.039252, Tree-Length: 115.318660, rate A &lt;-&gt; C: 1.037454, rate A &lt;-&gt; G: 1.791470, rate A &lt;-&gt; T: 1.154272, rate C &lt;-&gt; G: 0.616361, rate C &lt;-&gt; T: 2.463065, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.354004,0.236219,0.145031,0.264746)

### Read the best tree into R for each virus file. ( Use plot( ) to check. Display tree better with Dendroscope because there are lots of taxas)

``` r
tree_string_picornaviridae = read.tree("RAxML_bestTree.T1")
#plot(tree_string_picornaviridae)

tree_string_potyviridae = read.tree("RAxML_bestTree.T2")
#plot(tree_string_potyviridae)

tree_string_picornaviridae_5UTR=read.tree("RAxML_bestTree.T3")
#plot(tree_string_picornaviridae_5UTR)

tree_string_potyviridae_5UTR=read.tree("RAxML_bestTree.T4")
#plot(tree_string_potyviridae_5UTR)

#write.tree(tree_string_picornaviridae_5UTR,file="tree_3.newick")
```

### Next, root the unrooted trees prepared for later ancestral construction. Midpoint rooting locates the midpoint of the longest path between any two tips and putting the root in that location. After using midpoint root on the tree string, write.tree ( ) function in R could help check how new rooted trees look like.

``` r
library(phytools)
```

    ## Loading required package: maps

``` r
tree_string_picornaviridae_rooted=midpoint.root(tree_string_picornaviridae)

tree_string_potyviridae_rooted=midpoint.root(tree_string_potyviridae)

tree_string_picornaviridae_5UTR_rooted=midpoint.root(tree_string_picornaviridae_5UTR)

tree_string_potyviridae_5UTR_rooted=midpoint.root(tree_string_potyviridae_5UTR)

#write.tree(tree_string_picornaviridae_5UTR_rooted,file="tree5.newick")
```

### Ancestral state reconstruction

``` r
library(phytools)
library(magrittr)
library(kableExtra)
tree=read.tree("RAxML_bestTree.T3")
tree=midpoint.root(tree)
try1=fastBM(tree_string_picornaviridae_5UTR_rooted)
# The fitted paramter of reconstruction
fit<-fastAnc(tree_string_picornaviridae_5UTR_rooted,try1,vars=TRUE,CI=TRUE)

x<-fastBM(tree,internal=TRUE)
## ancestral states
a<-x[as.character(1:tree$Nnode+Ntip(tree))]
## tip data
x<-x[tree$tip.label]
fit<-fastAnc(tree,x,CI=TRUE)

plot(a,fit$ace,xlab="true states",ylab="estimated states")
```

![](virus_project_updated_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
phenogram(tree,x,ftype="off",spread.labels=FALSE)
```

![](virus_project_updated_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
#phenogram(tree_string_picornaviridae_5UTR_rooted,try1,fsize=0.6,spread.costs=c(1,0))
```

Show the ancestral construction using grasp
===========================================

``` r
#ancestral3=read.dna("tree3_N0-ancestors_GRASP.fasta",as.character=FALSE,format="fasta")
#ancestral3
#ancestral3_character=read.dna("tree3_N0-ancestors_GRASP.fasta",as.character=TRUE,format="fasta")

#ancestral4=read.dna("tree4_N0_ancestors_GRASP.fasta",as.character=FALSE,format="fasta")
#ancestral4
#ancestral4_character=read.dna("tree4_N0_ancestors_GRASP.fasta",as.character=TRUE,format="fasta")
```

### Analysis on the remaining complete genome files and 5 UTR files. They are \#\#\#potyviridae alone, picornaviridae alone, both potyviridae+picornaviridae

### All files have outgroup.

The construction of the best tree on the virus genome is based on total 20 maximized likelihood search and 1 bootstrap

The following are the code and result in the Raxml

##### potyviridae+picornaviridae with UTR

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_potyviridae\_5UTR\_aligned.fasta -n T7

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) = 0.254346, 0.251361, 0.231342, 0.262952, and the transition rate are A &lt;-&gt; C: 1.142242, rate A &lt;-&gt; G: 2.531875, rate A &lt;-&gt; T: 1.412656, rate C &lt;-&gt; G: 0.963031, rate C &lt;-&gt; T: 3.510598, rate G &lt;-&gt; T: 1.000000

##### potyviridae with UTR

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s potyviridae\_5UTR\_with\_outgroup\_aligned.fasta -n T8

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) = 0.344122, 0.237556, 0.154793, 0.263528, and the transition rate are rate A &lt;-&gt; C: 1.302173, rate A &lt;-&gt; G: 3.270232, rate A &lt;-&gt; T: 1.338095, rate C &lt;-&gt; G: 0.456309, rate C &lt;-&gt; T: 3.213243, rate G &lt;-&gt; T: 1.000000.

##### picornaviridae with UTR

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_5UTR\_with\_outgroup\_aligned.fasta -n T9

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.229238, 0.256425, 0.243877,and 0.270460, and the transition rate are A &lt;-&gt; C: 1.146259, rate A &lt;-&gt; G: 3.233354, rate A &lt;-&gt; T: 1.514142, rate C &lt;-&gt; G: 0.974224, rate C &lt;-&gt; T: 3.528592, rate G &lt;-&gt; T: 1.000000.

##### potyviridae with compelete genome and outgroup

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s potyviridae\_genome\_with\_outgroup\_aligned.fasta -n T10

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.313372, 0.196132, 0.233898, and 0.256599, and the transition rate are A &lt;-&gt; C: 1.302173, rate A &lt;-&gt; G: 3.270232, rate A &lt;-&gt; T: 1.338095, rate C &lt;-&gt; G: 0.456309, rate C &lt;-&gt; T: 3.213243, rate G &lt;-&gt; T: 1.000000.

##### picornaviridae with compelete genome and outgroup

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_genome\_with\_outgroup\_aligned.fasta -n T11

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.270188, 0.222885, 0.225654, and 0.281272, and the transition rate are A &lt;-&gt; C: 1.479092, rate A &lt;-&gt; G: 3.273508, rate A &lt;-&gt; T: 1.354061, rate C &lt;-&gt; G: 1.291275, rate C &lt;-&gt; T: 3.307613, rate G &lt;-&gt; T: 1.000000.

##### picornaviridae + potyviridae with compelete genome and outgroup

raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -\# 1 -s picornaviridae\_potyviridae\_genome\_aligned.fasta -n T12

The stationary distribution *π*<sub>(</sub>*A*, *C*, *G*, *T*) are 0.282127, 0.217977, 0.232850, and 0.267047, and the transition rate are rate A &lt;-&gt; C: 1.689051, rate A &lt;-&gt; G: 3.261407, rate A &lt;-&gt; T: 1.388764, rate C &lt;-&gt; G: 1.257696, rate C &lt;-&gt; T: 3.595042, rate G &lt;-&gt; T: 1.000000.

### Next, root the unrooted trees for later ancestral construction with midpoint rooting first. It locates the midpoint of the longest path between any two tips and putting the root in that location. Then write the rooted trees out. This method requires you to make the assumption that all of your sequences are evolving at the same rate - you should do so cautiously because this assumption does not hold for many biological datasets. In this case, the root is positioned at the midpoint between the two longest branches. If you have taxa that were not sampled at the same time point then a slight modification of this method would be required to take into account the time elapsed between samples.

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

# Without UTR region

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

### Then, we want to compare the result of mid-point rooting with rooting with the outgrouo to investigate the differences between them. The best possible outgroups are those available which are most closely related to our sequences of interest. If outgroups are too distantly related then they can be unreliable as they may be difficult to align reliably or have become saturated with substitutions.

##### For the potyviridae\_genome\_with\_outgroup\_aligned.fasta, Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. The IDs of these two virus are "NC\_001747.1","NC\_004750.1". For the potyviridae\_5UTR\_with\_outgroup\_aligned.fasta, Potato leafroll virus and Barley yellow dwarf virus (Luteoviridae family) are added as outgroups. The IDs of these two virus are "NC\_001747.1.1-177","NC\_004750.1.1-144".

##### For the picornaviridae\_genome\_with\_outgroup\_aligned.fasta, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. Their IDs are NC\_002551.1 and NC\_008311.1. For the picornaviridae\_5UTR\_with\_outgroup\_aligned.fasta, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. The IDs are NC\_002551.1.1-22 and NC\_008311.1.1-8.

##### For the picornaviridae\_potyviridae\_genome\_aligned.fasta, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. The IDs are NC\_002551.1 and NC\_008311.1. For the potyviridae\_picornaviridae\_5UTR.fasta, Vesicular exanthema of swine virus and Norovirus (Caliciviridae family) are added as outgroups. The IDs are NC\_002551.1.1-22 and NC\_008311.1.1-8

### The result shows that the specified outgroup is not monophyletic, so we need to select one outgroup.

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

### Next, use the RAXML to construct the ancestral sequence.

First, construct the ancestral sequence using mid-point root. With the GRASP function, it infers ancestral proteins from homologous protein sequences--a process known as ancestral sequence reconstruction. GRASP efficiently determines ancestral character states, and the most supported insertions and deletions. GRASP presents all this as partial-order graphs using a visual interface, connected to a phylogenetic tree, with which the user can interact. Detailed description see <http://grasp.scmb.uq.edu.au/guide>.

GRASP can perform reconstruction based on a variety of evolutionary models including JTT, Dayoff, LG, and WAG. JTT (Jones-Taylor-Thornton) is commonly used for the analysis of amino acid sequences.

Following code shows how to use the GRASP function with RAXML:

./raxmlHPC-PTHREADS -m PROTGAMMAJTT -p 23456 -n GRASPTutorial.nwk -s GRASPTutorial.aln.

PROTGAMMAJTT specifies to optimise substitution rates. it use a GAMMA model of rate heterogeneity and use the JTT amino acid substitution matrix. -p 23456 is the random number seed for inference. -n GRASPTutorial.nwk represents the name of the output file, and -s GRASPTutorial.aln is the name of the input alignment file.
