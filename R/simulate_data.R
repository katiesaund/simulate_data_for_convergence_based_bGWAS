suppressWarnings(library(ape))
suppressWarnings(library(caper))
suppressWarnings(library(phytools))
suppressWarnings(library(optparse))
suppressWarnings(library(phylolm))
suppressWarnings(library(dplyr))

source("tree.R")
source("discrete_trait_lib.R")
source("transition_edges.R")


# pseudocode / outline

# Initialize variables / read in user input
num_trees <- 2 # change to user defined input
num_tips <- 30 # change to user defined input
num_phenos <- 2 # change to user defined input

# Generate huge matrix of binary traits specific to trees
tree_list <- generate_trees(num_trees, num_tips)
binary_AR_mat_list <- generate_disc_mat(tree_list)

# Separate traits into BM and WN
phylo_signal_list <- calculate_phylo_signal(tree_list, binary_AR_mat_list) 

# Select BM and WN phenotypes
BM_phenotype_names_list <- select_BM_traits(binary_AR_mat_list, phylo_signal_list, num_phenos) 
WN_phenotype_names_list <- select_WN_traits(binary_AR_mat_list, phylo_signal_list, num_phenos) # We're not getting enough WN, presumably because the traits are being created by evolving on the tree. I'll need to randomize them somehow.
phenotype_names_list <- combine_phenotype_names_lists(BM_phenotype_names_list, WN_phenotype_names_list)
# Select genotypes and phenotypes
genotype_AR_mat_list <- subsample_to_genotypes(binary_AR_mat_list,
                                            phylo_signal_list,
                                            lower_bound = -1.5, 
                                            upper_bound = 1.5, 
                                            num_genos = 100)
BM_phenotype_AR_mat_list <- subsample_to_phenotypes(binary_AR_mat_list, BM_phenotype_names_list)
WN_phenotype_AR_mat_list <- subsample_to_phenotypes(binary_AR_mat_list, WN_phenotype_names_list)

BM_phenotype_sync_trans_list <- find_transition_edges(tree_list, BM_phenotype_AR_mat_list, "discrete")
WN_phenotype_sync_trans_list <- find_transition_edges(tree_list, WN_phenotype_AR_mat_list, "discrete")
genotype_sync_trans_list <- find_transition_edges(tree_list, genotype_AR_mat_list, "discrete")
# Everything below is psuedocode

genotype_phyc_trans_list <- convert_to_phyc_trans(genotype_sync_trans_list)

# Calculate gamma
phyc_gamma_list <- calc_phyc_gamma(genotype_phyc_trans_list, phenotype_mat_list)
sync_gamma_list <- calc_sync_gamma(genotype_sync_trans_list, phenotype_sync_trans_list)

genotype_keeper_list <- subset_genotype(phyc_gamma_list, sync_gamma_list)
phenotype_keeper_list <- subset_phenotype(phyc_gamma_list, sync_gamma_list)

save_genotypes(genotype_keeper_list)
save_phenotypes(phenotype_keeper_list)

# End pseudocode
