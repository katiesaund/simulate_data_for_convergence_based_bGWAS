suppressWarnings(library(ape))
suppressWarnings(library(caper))
suppressWarnings(library(phytools))
suppressWarnings(library(phylolm))

source("../../simulate_data_for_convergence_based_bGWAS/R/tree.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/discrete_trait_lib.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/transition_edges.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/gamma.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/save_data.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/keep_interesting_genotypes.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/ancestral_reconstruction.R")

# Initialize variables / read in user input
num_trees <- 3 # change to user defined input
num_tips <- 100 # change to user defined input
num_phenos <- 2 # change to user defined input
num_start_trait <- 100 # change to user defined input
tree_edge_multiplier <- 100 # change to user defined input

# Generate huge matrix of binary traits specific to trees
tree_list <- generate_trees(num_trees, num_tips, tree_edge_multiplier)

binary_AR_mat_list <- generate_disc_mat(tree_list, num_start_trait)

binary_AR_mat_list <- add_WN(binary_AR_mat_list, tree_list) # Due to WN stuff the ancestral reconstructions are now wrong!

binary_AR_and_conf_mat <- ancestral_reconstruction(binary_AR_mat_list, tree_list)
binary_conf_mat_list <- binary_AR_and_conf_mat$conf_mat
binary_AR_mat_list <- binary_AR_and_conf_mat$AR_mat
# TODO - Add the confidence information to both gamma calc and the subsetting steps (select_geno_within_range, subsample_to_phenotypes)

# Separate traits into BM and WN
phylo_signal_list <- calculate_phylo_signal(tree_list, binary_AR_mat_list)

# #temp stuff
# phylo_signal_list[[1]][1:2] <- phylo_signal_list[[2]][1:2] <- 0
# phylo_signal_list[[1]][3:4] <- phylo_signal_list[[2]][3:4] <- 1
# #phylo_signal_list[[3]][1:2] <- 0
# # phylo_signal_list[[3]][3:4] <- 1
# # save(phylo_signal_list, file = "../data/phylo_signal_list.Rdata")
# # end temp stuff

# Select BM and WN phenotypes
BM_phenotype_names_list <- select_BM_traits(binary_AR_mat_list, phylo_signal_list, num_phenos)
WN_phenotype_names_list <- select_WN_traits(binary_AR_mat_list, phylo_signal_list, num_phenos)
phenotype_names_list <- combine_phenotype_names_lists(BM_phenotype_names_list, WN_phenotype_names_list)
# Select genotypes and phenotypes
genotype_AR_and_conf_mat_list <- select_geno_within_range(binary_AR_mat_list,
                                                          binary_conf_mat_list,
                                                          phylo_signal_list,
                                                          lower_bound = -1.5,
                                                          upper_bound = 1.5,
                                                          min_genos = 10)
genotype_AR_mat_list <- genotype_AR_and_conf_mat_list$AR_mat
genotype_conf_mat_list <- genotype_AR_and_conf_mat_list$conf_mat
BM_phenotype_AR_and_conf_mat_list <-
  subsample_to_phenotypes(binary_AR_mat_list,
                          binary_conf_mat_list,
                          BM_phenotype_names_list)
BM_phenotype_AR_mat_list <- BM_phenotype_AR_and_conf_mat_list$AR_mat
BM_phenotype_conf_mat_list <- BM_phenotype_AR_and_conf_mat_list$conf_mat
WN_phenotype_AR_and_conf_mat_list <-
  subsample_to_phenotypes(binary_AR_mat_list,
                          binary_conf_mat_list,
                          WN_phenotype_names_list)
WN_phenotype_AR_mat_list <- WN_phenotype_AR_and_conf_mat_list$AR_mat
WN_phenotype_conf_mat_list <- WN_phenotype_AR_and_conf_mat_list$conf_mat


# Order everything by edges instead of by tips then nodes ----
# BM pheno
BM_pheno_recon_by_edge_list <- prep_pheno_recon_edges(BM_phenotype_AR_mat_list,
                                                      tree_list)
BM_pheno_recon_conf_by_edge_list <- reorder_tip_and_node_to_edge_lists(BM_phenotype_conf_mat_list, tree_list)
BM_pheno_sync_trans_by_edge_list <- find_transition_edges(tree_list, BM_phenotype_AR_mat_list, "discrete")

# WN pheno
WN_pheno_recon_by_edge_list <- prep_pheno_recon_edges(WN_phenotype_AR_mat_list, tree_list)
WN_pheno_recon_conf_by_edge_list <- reorder_tip_and_node_to_edge_lists(WN_phenotype_conf_mat_list, tree_list)
WN_pheno_sync_trans_by_edge_list <- find_transition_edges(tree_list, WN_phenotype_AR_mat_list, "discrete")

# Geno
genotype_sync_trans_by_edge_list <- find_transition_edges(tree_list, genotype_AR_mat_list, "discrete")
genotype_phyc_trans_by_edge_list <- convert_to_phyc_trans(genotype_AR_mat_list, genotype_sync_trans_by_edge_list)
genotype_recon_conf_by_edge_list <- reorder_tip_and_node_to_edge_lists(genotype_conf_mat_list, tree_list)

# TODO ID high confidence phenotype recon and phenotype/genotype transitions by edges
BM_pheno_sync_trans_conf_by_edge_list <- assign_high_confidence_to_transition_edges()
WN_pheno_sync_trans_conf_by_edge_list
genotype_phyc_trans_conf_by_edge_list
genotype_sync_trans_conf_by_edge_list
# TODO incorporate high confidence information into the gamma calcuations next


# Calculate gamma
BM_phyc_gamma_list <- calc_phyc_gamma_list(tree_list, genotype_phyc_trans_by_edge_list, BM_pheno_recon_by_edge_list, phyc_geno_trans_BM_recon_conf_by_edge_list)
WN_phyc_gamma_list <- calc_phyc_gamma_list(tree_list, genotype_phyc_trans_by_edge_list, WN_pheno_recon_by_edge_list, phyc_geno_trans_WN_recon_conf_by_edge_list)
BM_sync_gamma_list <- calc_sync_gamma_list(tree_list, genotype_sync_trans_by_edge_list, BM_pheno_sync_trans_by_edge_list)
WN_sync_gamma_list <- calc_sync_gamma_list(tree_list, genotype_sync_trans_by_edge_list, WN_pheno_sync_trans_by_edge_list)

# Subset genotype matrix for each tree - phenotype combo
BM_phyc_genotype_keeper_list <- keep_good_genotypes(BM_phyc_gamma_list)
BM_sync_genotype_keeper_list <- keep_good_genotypes(BM_sync_gamma_list)
WN_phyc_genotype_keeper_list <- keep_good_genotypes(WN_phyc_gamma_list)
WN_sync_genotype_keeper_list <- keep_good_genotypes(WN_sync_gamma_list)

# Save data
save_data(tree_list,
          genotype_AR_mat_list,
          genotype_phyc_trans_by_edge_list,
          genotype_sync_trans_by_edge_list,
          BM_phenotype_AR_mat_list,
          BM_pheno_recon_by_edge_list,
          BM_phyc_gamma_list,
          BM_sync_gamma_list,
          WN_phenotype_AR_mat_list,
          WN_pheno_recon_by_edge_list,
          WN_phyc_gamma_list,
          WN_sync_gamma_list,
          BM_phyc_genotype_keeper_list,
          BM_sync_genotype_keeper_list,
          WN_phyc_genotype_keeper_list,
          WN_sync_genotype_keeper_list)
