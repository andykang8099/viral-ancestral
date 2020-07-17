virus\_project\_initial\_files
================

### Reading files

First, read four virus files: - picornaviridae\_aligned - potyviridae\_aligned - picornaviridae\_5UTR\_aligned, and - potyviridae\_5UTR\_aligned

into dataframes, which are prepared for later ancestral construction.

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

### RAxML runs

Second, use Raxml ( 20 ML search + bootstrapping ) to get the tree with best likelihood. This part of commands is from bash (Terminal).

Note1: the file size of picornaviridae\_aligned and potyviridae\_aligned is too large. As a result, I limit the bootstrap search time to 1. The search times for Maximal Likelihood (ML) are always 20 times. The estimated best tree could be further improved if the bootstrap times increase. The bootstrap search times for picornaviridae\_5UTR\_aligned and potyviridae\_5UTR\_aligned are 6 and 10 times.

Note2: For picornaviridae\_5UTR\_aligned and potyviridae\_5UTR\_aligned, the names of some virus sequences contain colon (:), so I tried to delete the colon when the virus name contain it in order to make the fas file executable in Raxml. The files with deleted colon are picornaviridae\_5UTR\_aligned1 and potyviridae\_5UTR\_aligned1.

In the picornaviridae\_5UTR\_aligned, names of five virus sequences contain colon (:). They are - NC\_004421.1\_Aichivirus\_B\_genomic\_RNA\_complete\_genome\_strain\_U-1,NC\_033695.1\_Enterovirus\_AN12\_genomic\_RNA\_complete\_genome\_strain\_AN12, - NC\_027919.1\_Kobuvirus\_cattle/Kagoshima-1-22-KoV/2014/JPN\_genomic\_RNA\_nearly\_complete\_genome\_strain\_Kagoshima-1-22-KoV/2014/JPN, - NC\_027918.1\_Kobuvirus\_cattle/Kagoshima-2-24-KoV/2015/JPN\_genomic\_RNA\_complete\_genome\_strain\_Kagoshima-2-24-KoV/2015/JPN, - NC\_025474.1\_Crohivirus\_gene\_for\_polyprotein\_complete\_cds\_strain\_ZM54.

For potyviridae\_5UTR\_aligned, the names of virus sequences contain colon (:) are
- NC\_014327.1\_Pepper\_yellow\_mosaic\_virus\_RNA\_complete\_genome\_strain\_DF,NC\_030391.1\_Wild\_onion\_symptomless\_virus\_gene\_for\_polyprotein\_complete\_cds\_isolate\_TUR256-1, - NC\_021786.1\_Habenaria\_mosaic\_virus\_genomic\_RNA\_complete\_genome\_isolate\_Ha-1

### Bash code for rapid bootstrapping in Raxml:

Note: In the RAXML, I use T1, T2, T3, and T4 to differentiate four trees from four different files. In the folder, I change the suffix of all the tree files to .tree to make it more clear.

All the trees are constructed using the GTR and GAMMA model. 12345 refers to the randome seed via -p and -x. -n restricts the name of the output files. -f a represents the fast bootstraping in the RAXML. The advantage is that it allows to do a complete analysis (ML search + Bootstrapping) in one single step.

    raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 1 -s picornaviridae_aligned.fasta -n T1  

    Output file: RAxML_bestTree_picornaviridae_aligned.tree

    raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 1 -s potyviridae_aligned.fasta -n T2 

    Output file: RAxML_bestTree_potyviridae_aligned.tree

    raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 6 -s picornaviridae_5UTR_aligned.fas -n T3 

    Output file: RAxML_bestTree_picornaviridae_5UTR_aligned.tree

    raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 10 -s potyviridae_5UTR_aligned.fas -n T4 

    Output file: RAxML_bestTree_potyviridae_5UTR_aligned.tree

#### Relevant parameters for each ML search + bootstrapping best tree estimation

picornaviridae\_aligned.fasta: alpha: 4.934338, Tree-Length: 126.910679, rate A &lt;-&gt; C: 1.448674, rate A &lt;-&gt; G: 3.216025, rate A &lt;-&gt; T: 1.314619, rate C &lt;-&gt; G: 1.269701, rate C &lt;-&gt; T: 3.239043, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.270869, 0.221802, 0.224860, 0.282469).

potyviridae\_aligned.fasta: alpha: 0.932693, Tree-Length: 132.828870, rate A &lt;-&gt; C: 2.486990, rate A &lt;-&gt; G: 3.739170, rate A &lt;-&gt; T: 1.369314, rate C &lt;-&gt; G: 1.513431, rate C &lt;-&gt; T: 5.417214, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.313614,0.195751,0.233819,0.256815)

picornaviridae\_5UTR\_aligned.fas: alpha: 5.094113, Tree-Length: 110.802014, rate A &lt;-&gt; C: 1.172756, rate A &lt;-&gt; G: 2.585556, rate A &lt;-&gt; T: 1.493185, rate C &lt;-&gt; G: 0.866442, rate C &lt;-&gt; T: 2.625791, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.225632,0.257601,0.242809,0.273959)

potyviridae\_5UTR\_aligned.fas: alpha: 1.039252, Tree-Length: 115.318660, rate A &lt;-&gt; C: 1.037454, rate A &lt;-&gt; G: 1.791470, rate A &lt;-&gt; T: 1.154272, rate C &lt;-&gt; G: 0.616361, rate C &lt;-&gt; T: 2.463065, rate G &lt;-&gt; T: 1.000000, pi(A,C,G,T) = c(0.354004,0.236219,0.145031,0.264746)

### Estimate ancestral sequence

#### Read trees in R

Read the best tree into R for each virus file. ( Use `plot( )` to check. Display tree better with Dendroscope because there are lot of taxa)

``` r
tree_string_picornaviridae = read.tree("RAxML_bestTree_picornaviridae_aligned.tree")
#plot(tree_string_picornaviridae)

tree_string_potyviridae = read.tree("RAxML_bestTree_potyviridae_aligned.tree")
#plot(tree_string_potyviridae)

tree_string_picornaviridae_5UTR=read.tree("RAxML_bestTree_picornaviridae_5UTR_aligned.tree")
#plot(tree_string_picornaviridae_5UTR)

tree_string_potyviridae_5UTR=read.tree("RAxML_bestTree_potyviridae_5UTR_aligned.tree")
#plot(tree_string_potyviridae_5UTR)
```

#### Root trees in R (midpoint)

Next, root the unrooted trees prepared for later ancestral construction. Midpoint rooting locates the midpoint of the longest path between any two tips and putting the root in that location. After using midpoint root on the tree string, `write.tree( )` function in R could help check how new rooted trees look like.

``` r
library(phytools)
```

    ## Loading required package: maps

``` r
library(ape)
tree_string_picornaviridae_rooted=midpoint.root(tree_string_picornaviridae)
write.tree(tree_string_picornaviridae_rooted,file="RAxML_bestTree_picornaviridae_aligned_midpoint.tree")

#output file: RAxML_bestTree_picornaviridae_aligned_midpoint.tree

tree_string_potyviridae_rooted=midpoint.root(tree_string_potyviridae)
write.tree(tree_string_potyviridae_rooted,file="RAxML_bestTree_potyviridae_aligned_midpoint.tree")

#output file: RAxML_bestTree_potyviridae_aligned_midpoint.tree

tree_string_picornaviridae_5UTR_rooted=midpoint.root(tree_string_picornaviridae_5UTR)

write.tree(tree_string_picornaviridae_5UTR_rooted,file="RAxML_bestTree_picornaviridae_5UTR_aligned_midpoint.tree")

#output file: RAxML_bestTree_picornaviridae_5UTR_aligned_midpoint.tree

tree_string_potyviridae_5UTR_rooted=midpoint.root(tree_string_potyviridae_5UTR)
write.tree(tree_string_picornaviridae_5UTR_rooted,file="RAxML_bestTree_potyviridae_5UTR_aligned_midpoint.tree")

#output file: RAxML_bestTree_potyviridae_5UTR_aligned_midpoint.tree
```

#### Ancestral state reconstruction, use the aligned files of picornaviridae virus as an example to analyze the state of each nucleotide of the ancestral virus

``` r
library(phytools)
library(magrittr)
library(kableExtra)
result1=fastBM(tree_string_picornaviridae_rooted)
# The fitted paramter of reconstruction
fit1<-fastAnc(tree_string_picornaviridae_rooted,result1,vars=TRUE,CI=TRUE)
x<-fastBM(tree_string_picornaviridae_rooted,internal=TRUE)

## ancestral states
as<-x[as.character(1:tree_string_picornaviridae_rooted$Nnode+Ntip(tree_string_picornaviridae_rooted))]

## tip data
x<-x[tree_string_picornaviridae_rooted$tip.label]
fit<-fastAnc(tree_string_picornaviridae_rooted,x,CI=TRUE)

plot(as,fit$ace,xlab="true states",ylab="estimated states")
```

![](virus_project_initial_files_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
# Draw the graphs to show the state using phenogram function in R
phenogram(tree_string_picornaviridae_rooted,x,ftype="off",spread.labels=FALSE)
```

![](virus_project_initial_files_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
#phenogram(tree_string_picornaviridae_rooted,result1,fsize=0.6,spread.costs=c(1,0))
```

### Next step is to use RAXML to construct the ancestral sequence:

RAXML allows to choose differents models for ancestral sequences reconstruction, which include GTR, GAMMA, JTT, Dayoff, LG, and WAG model. Using the joint reconstrunction, it infers the most probable combination of all ancestral states. The assignment of ancestral states (across the whole phylogenetic tree) optimises at the joint state with the greatest likelihood. It creates the most probable ancestors, given the extant sequences and their phylogenetic relationships. The inference is valid for the whole tree, and it is practical to compare ancestors with one another as they represent snapshots from the same evolutionary history.

Note: In order to make the output files more clear, the suffix names ending with A will convert to fasta instead of the phylip files derived from the RAXML. By this way, all the ancestral sequences files will have the unified name to compare

#### The bash code to produce the ancestral sequence using RAXML:

raxmlHPC -f A -p 12345 -t RAxML\_bestTree\_picornaviridae\_5UTR\_aligned\_midpoint.tree -s picornaviridae\_5UTR\_aligned.fas -m GTRGAMMA -n A3

Output files: picornaviridae\_5UTR\_aligned\_midpoint\_AncestralSeq.fasta

raxmlHPC -f A -p 12345 -t RAxML\_bestTree\_potyviridae\_5UTR\_aligned\_midpoint.tree -s potyviridae\_5UTR\_aligned.fas -m GTRGAMMA -n A4

Output files: potyviridae\_5UTR\_aligned\_midpoint\_AncestralSeq.fasta

Note: For the first two files, maybe because the original files are large, it may take long time to get the result, I will try to run it again to get result.

Output files
============

In `initial-files` folder:

#### Tree files:

-   RAxML\_bestTree\_picornaviridae\_5UTR\_aligned\_midpoint.tree
-   RAxML\_bestTree\_picornaviridae\_5UTR\_aligned.tree
-   RAxML\_bestTree\_picornaviridae\_aligned\_midpoint.tree
-   RAxML\_bestTree\_picornaviridae\_aligned.tree
-   RAxML\_bestTree\_potyviridae\_5UTR\_aligned\_midpoint.tree
-   RAxML\_bestTree\_potyviridae\_5UTR\_aligned.tree
-   RAxML\_bestTree\_potyviridae\_aligned\_midpoint.tree
-   RAxML\_bestTree\_potyviridae\_aligned.tree

#### Ancestral Sequence files:

-picornaviridae\_5UTR\_aligned\_midpoint\_AncestralSeq.fasta -potyviridae\_5UTR\_aligned\_midpoint\_AncestralSeq.fasta

### Reference:

<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004763>

<http://vu-wp0.s3.amazonaws.com/wp-content/uploads/sites/191/pdfs/2011_Rokas_CPMB.pdf>
