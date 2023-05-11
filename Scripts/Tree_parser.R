#load packages
suppressMessages(suppressWarnings(library("ape")))
suppressMessages(suppressWarnings(library("dplyr")))
suppressMessages(suppressWarnings(library("phytools")))

#load the tree
mytree <- read.tree("Final_ALL_verification_alignment.aln.treefile")


#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")


#import the name of known target family genes
known_target_family_genes <- scan("Known_family_genes_id.txt", what="character", sep="\n")
outgroups_seqs <- scan("Known_outgroup_genes_id.txt", what="character", sep="\n")


#Root the tree at the common ancestor of outgroup sequences
MRCA_outgroup <- findMRCA(mytree_label, tips=outgroups_seqs, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_outgroup, resolve.root= TRUE) #root


#Check the MRCA of target family genes
MRCA_Target <- findMRCA(mytree_rooted, tips=known_target_family_genes, type="node")


#grep all tips from these MRCS
target_family_genes <- as.vector(extract.clade(mytree_rooted, MRCA_Target)$tip)


#Remove already known outgroups genes
Current_species_target_genes <- setdiff(target_family_genes, known_target_family_genes)


write(x=Current_species_target_genes, file="Current_species_target.txt")
