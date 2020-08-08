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
