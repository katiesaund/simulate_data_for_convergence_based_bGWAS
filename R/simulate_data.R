suppressWarnings(library(ape))
suppressWarnings(library(caper))
suppressWarnings(library(phytools))
suppressWarnings(library(optparse))
suppressWarnings(library(phylolm))
suppressWarnings(library(dplyr))

source("generate_disc_traits.R")


# pseudocode / outline

# Initialize variables / read in user input
num_trees <- 1 # change to user defined input

# Generate huge matrix of binary traits specific to trees
tree_list <- generate_trees(num_trees)
binary_mat_list <- generate_disc_mat(trees)

# Separate traits into BM and WN
phylo_signal <- calculate_phylo_signal(binary_mat_list)
BM_mat_list <- select_BM_traits(binary_mat_list, phylo_signal)
WN_mat_list <- select_WN_traits(binary_mat_list, phylo_signal)

# Select genotypes and phenotypes
genotype_mat_list <- ???()
phenotype_mat_list <- ???()
genotype_AR_list <- ??()
phenotype_AR_list <- ?()

# Identify transition edges
genotype_sync_trans_list <- find_transition_edges(tree_list, genotype_mat_list, genotype_AR_list, "discrete")
phenotype_sync_trans_list <- find_transition_edges(tree_list, phenotype_mat_list, phenotype_AR_list, "discrete")
genotype_phyc_trans_list <- convert_to_phyc_trans(genotype_sync_trans_list)

# Calculate gamma
phyc_gamma_list <- calc_phyc_gamma(genotype_phyc_trans_list, phenotype_mat_list)
sync_gamma_list <- calc_sync_gamma(genotype_sync_trans_list, phenotype_sync_trans_list)

genotype_keeper_list <- subset_genotype(phyc_gamma_list, sync_gamma_list)
phenotype_keeper_list <- subset_phenotype(phyc_gamma_list, sync_gamma_list)

save_genotypes(genotype_keeper_list)
save_phenotypes(phenotype_keeper_list)

# End pseudocode
