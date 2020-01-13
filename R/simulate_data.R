suppressWarnings(library(ape))
suppressWarnings(library(caper))
suppressWarnings(library(phytools))
suppressWarnings(library(phylolm))

source("../../simulate_data_for_convergence_based_bGWAS/R/tree.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/discrete_trait_lib.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/transition_edges.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/gamma.R")

# Initialize variables / read in user input
num_trees <- 3 # change to user defined input
num_tips <- 100 # change to user defined input
num_phenos <- 2 # change to user defined input
num_start_trait <- 10000 # change to user defined input

# Generate huge matrix of binary traits specific to trees
tree_list <- generate_trees(num_trees, num_tips)
binary_AR_mat_list <- generate_disc_mat(tree_list, num_start_trait)

# Separate traits into BM and WN
phylo_signal_list <- calculate_phylo_signal(tree_list, binary_AR_mat_list) 

#temp stuff
phylo_signal_list[[1]][1:2] <- phylo_signal_list[[2]][1:2]  <- 0
phylo_signal_list[[1]][3:4] <- phylo_signal_list[[2]][3:4] <- 1
save(phylo_signal_list, file = "../data/phylo_signal_list.Rdata")
# end temp stuff

# Select BM and WN phenotypes
BM_phenotype_names_list <- select_BM_traits(binary_AR_mat_list, phylo_signal_list, num_phenos) 
WN_phenotype_names_list <- select_WN_traits(binary_AR_mat_list, phylo_signal_list, num_phenos) # We're not getting enough WN, presumably because the traits are being created by evolving on the tree. I'll need to randomize them somehow.
phenotype_names_list <- combine_phenotype_names_lists(BM_phenotype_names_list, WN_phenotype_names_list)
# Select genotypes and phenotypes
genotype_AR_mat_list <- select_geno_within_range(binary_AR_mat_list,
                                                 phylo_signal_list,
                                                 lower_bound = -1.5, 
                                                 upper_bound = 1.5, 
                                                 num_genos = 100)
BM_phenotype_AR_mat_list <- subsample_to_phenotypes(binary_AR_mat_list, BM_phenotype_names_list)
WN_phenotype_AR_mat_list <- subsample_to_phenotypes(binary_AR_mat_list, WN_phenotype_names_list)

BM_phenotype_sync_trans_list <- find_transition_edges(tree_list, BM_phenotype_AR_mat_list, "discrete")
WN_phenotype_sync_trans_list <- find_transition_edges(tree_list, WN_phenotype_AR_mat_list, "discrete")
genotype_sync_trans_list <- find_transition_edges(tree_list, genotype_AR_mat_list, "discrete")
genotype_phyc_trans_list <- convert_to_phyc_trans(genotype_AR_mat_list, genotype_sync_trans_list)

# Prep phenotype reconstruction by edges 
BM_pheno_recon_by_edge_list <- prep_pheno_recon_edges(BM_phenotype_AR_mat_list, 
                                                      tree_list)
WN_pheno_recon_by_edge_list <- prep_pheno_recon_edges(WN_phenotype_AR_mat_list, 
                                                      tree_list)

# Calculate gamma
BM_phyc_gamma_list <- calc_phyc_gamma_list(tree_list, genotype_phyc_trans_list, BM_pheno_recon_by_edge_list)
WN_phyc_gamma_list <- calc_phyc_gamma_list(tree_list, genotype_phyc_trans_list, WN_pheno_recon_by_edge_list)
BM_sync_gamma_list <- calc_sync_gamma_list(tree_list, genotype_sync_trans_list, BM_phenotype_sync_trans_list)
WN_sync_gamma_list <- calc_sync_gamma_list(tree_list, genotype_sync_trans_list, WN_phenotype_sync_trans_list)

# Save data
save_data(tree_list,
          genotype_AR_mat_list,
          genotype_phyc_trans_list, 
          genotype_sync_trans_list,
          BM_phenotype_AR_mat_list, 
          BM_pheno_recon_by_edge_list, 
          BM_phyc_gamma_list, 
          BM_sync_gamma_list,
          WN_phenotype_AR_mat_list, 
          WN_pheno_recon_by_edge_list, 
          WN_phyc_gamma_list, 
          WN_sync_gamma_list)

# Everything below is psuedocode
# genotype_keeper_list <- subset_genotype(phyc_gamma_list, sync_gamma_list)
# phenotype_keeper_list <- subset_phenotype(phyc_gamma_list, sync_gamma_list)
# save_genotypes(genotype_keeper_list)
# save_phenotypes(phenotype_keeper_list)

# End pseudocode
