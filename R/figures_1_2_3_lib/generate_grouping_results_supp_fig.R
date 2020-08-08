# Goal: create New Supplementary Figure S3 which will highlight how the grouping
#       method can create genotypes with higher significance than the individual
#       constituent loci. 

# What we need: Single locus genotype transition plots for the 6 SNPs. 
# Pre-AR: need simluation genotypes #, #, #
# Post-AR: need simulation genotypes that map to GENE 1. 

# Import functions ----
source("plot_lib.R")

# Read in gene key
gene_key <- read.table("../../data/NEW_group_key_for_phyc_BM_tree_1_pheno_1.tsv")

# Read in genotype matrix
geno_mat <- read.table("../../data/simulated_genotype_for_discrete_pheno_BM_tree_1_pheno_1.tsv", 
                       sep = "\t")

# Read in tree
tree <- read.tree("../../data/simulated_discrete_tree_1.tree")

# Get SNPs of interest for pre- and post-ar approaches
post_ar_snp_ids <- as.character(gene_key[gene_key$group == "GENE1", 1])
post_ar_geno_mat <- geno_mat[, colnames(geno_mat) %in% post_ar_snp_ids, drop = FALSE]
# Let's plot a tree where the tips show the presence/absence of the tree snps

post_ar_snp1_tips <- post_ar_geno_mat[, 1, drop = FALSE]
post_ar_snp1_tips <- post_ar_snp1_tips[post_ar_snp1_tips == 1, , drop = FALSE]
post_ar_snp1_tips <- row.names(post_ar_snp1_tips)

post_ar_snp2_tips <- post_ar_geno_mat[, 2, drop = FALSE]
post_ar_snp2_tips <- post_ar_snp2_tips[post_ar_snp2_tips == 1, , drop = FALSE]
post_ar_snp2_tips <- row.names(post_ar_snp2_tips)

post_ar_snp3_tips <- post_ar_geno_mat[, 3, drop = FALSE]
post_ar_snp3_tips <- post_ar_snp3_tips[post_ar_snp3_tips == 1, , drop = FALSE]
post_ar_snp3_tips <- row.names(post_ar_snp3_tips)

pdf("../../img/supp_post_ar_GENE1_snps.pdf", 
    width = 8, 
    height = 8)
plot(tree, show.tip.label = FALSE, use.edge.length = FALSE)
tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% post_ar_snp1_tips],
          pch = 19, 
          cex = .5, 
          col = "blue", 
          offset = 1)
tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% post_ar_snp2_tips],
          pch = 19, 
          cex = .5, 
          col = "red", 
          offset = 2)
tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% post_ar_snp3_tips],
          pch = 19, 
          cex = .5, 
          col = "purple", 
          offset = 3)
dev.off()

# Repeat for pre-ar
# Get SNPs of interest for pre- and post-ar approaches
pre_ar_snp_ids <- as.character(gene_key[gene_key$group == "GENE1", 1])
pre_ar_geno_mat <- geno_mat[, colnames(geno_mat) %in% pre_ar_snp_ids, drop = FALSE]
# Let's plot a tree where the tips show the presence/absence of the tree snps

pre_ar_snp1_tips <- pre_ar_geno_mat[, 1, drop = FALSE]
pre_ar_snp1_tips <- pre_ar_snp1_tips[pre_ar_snp1_tips == 1, , drop = FALSE]
pre_ar_snp1_tips <- row.names(pre_ar_snp1_tips)

pre_ar_snp2_tips <- pre_ar_geno_mat[, 2, drop = FALSE]
pre_ar_snp2_tips <- pre_ar_snp2_tips[pre_ar_snp2_tips == 1, , drop = FALSE]
pre_ar_snp2_tips <- row.names(pre_ar_snp2_tips)

pre_ar_snp3_tips <- pre_ar_geno_mat[, 3, drop = FALSE]
pre_ar_snp3_tips <- pre_ar_snp3_tips[pre_ar_snp3_tips == 1, , drop = FALSE]
pre_ar_snp3_tips <- row.names(pre_ar_snp3_tips)

pdf("../../img/supp_pre_ar_GENE1_snps.pdf", 
    width = 8, 
    height = 8)
plot(tree, show.tip.label = FALSE, use.edge.length = FALSE)
tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% pre_ar_snp1_tips],
          pch = 19, 
          cex = .5, 
          col = "green3", 
          offset = 1)
tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% pre_ar_snp2_tips],
          pch = 19, 
          cex = .5, 
          col = "magenta3", 
          offset = 2)
tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% pre_ar_snp3_tips],
          pch = 19, 
          cex = .5, 
          col = "orange", 
          offset = 3)
dev.off()


