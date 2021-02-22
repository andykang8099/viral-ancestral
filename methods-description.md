# 5UTR ancestral sequence

Replicase 2: with mitovirus

- tree reconstructed from replicase 2 (mitovirus)
    - potyvirus
    - picornavirus
    - both
- the analyses is descrbed in replicase_mitovirus_outgroup.md
- we used raxml: "Each tree file is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF. According to some research, JTTF model is more reliable in the construction of phylogenetic tree of protein sequence."
- same settings for the 3 cases
- we do two trees: nucleotide and protein for each case
- we root the trees on the outgroup (mitovirus). we also rooted at midpoint, but not used


- for the ancestral sequence reconstruction, we use the 5UTR sequence
    - potyvirus
    - picornavirus
    - both
- we removed the outgroup prior to ancestral reconstruction
- we use raxml again under the same GTRGAMMA model
- we also reconstructed the ancestral sequence for the replicase protein with raxml and JTT model
- we used both trees (nucleotides/protein), but we are only focusing on the nucleotides tree for now
- we had to remove sites with only gaps


---------------
Replicase and triticum (only for picornaviridae)

- tree reconsructed from replicase `replicase_triticum_mitovirus_outgroup_code.md`
- we used raxml: "Each tree file is constructed using 20 ML search and total 50 bootstrap simulation. The model of nucleotide files are GTR and GAMMA models. The model of the protein files are JTTF. According to some research, JTTF model is more reliable in the construction of phylogenetic tree of protein sequence."

- for the ancestral sequence reconstruction for 5UTR, the analyses are in `Analysis_with_or_without_triticum.md`
- we removed the outgroup, and we did two reconstructions: 1) with triticum and 2) without triticum
- we use raxml with the same GTRGAMMA
- we had to remove sites with only gaps
- we want to compare this ancestal sequence to the one we got from the replicase 2 (mitovirus) tree


# Complete genome

- we use the tree reconstructed from replicase 2 (mitovirus)
- we reconstruct the ancestral sequence for the whole genome (after removing the outgroup) for the three cases
    - potyvirus
    - picornavirus
    - both
- the analyses are in `Complete_Genome_Without_Outgroup.md`
- we used raxml again with GTRGAMMA model
- we only have nucleotides


# Next steps
- We will receive a smaller data to reconstruct tree and reconstruct ancestral sequence
- We want to get the ancestral sequence at every internal node