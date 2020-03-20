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
```
