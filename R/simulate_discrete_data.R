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
source("../../simulate_data_for_convergence_based_bGWAS/R/high_confidence.R")

# Initialize variables / read in user input
args <- commandArgs(trailingOnly = TRUE)
num_trees <- as.numeric(args[1]) 
num_phenos <- as.numeric(args[2]) 
num_tips <- as.numeric(args[3]) 
num_start_trait <- as.numeric(args[4])

# Defaults that shouldn't need user input
tree_edge_multiplier <- 100
bootstrap_threshold <- 70
bin_size <- 20

# Generate huge matrix of binary traits specific to trees
tree_list <- generate_trees(num_trees, num_tips, tree_edge_multiplier)

binary_AR_df_list <- generate_binary_df_list(tree_list, num_start_trait)
print("Finish generating first discrete matrix")

binary_AR_df_list <- add_WN(binary_AR_df_list, tree_list) # Due to WN stuff the ancestral reconstructions are now wrong!

binary_AR_and_conf_mat <- ancestral_reconstruction(binary_AR_df_list, tree_list) # So fix the ancestral reconstructions here
binary_conf_mat_list <- binary_AR_and_conf_mat$conf_mat
binary_AR_df_list <- binary_AR_and_conf_mat$AR_mat
print("Finish ancestral reconstruction")

# Separate traits into BM and WN
phylo_signal_list <- calculate_phylo_signal(tree_list, binary_AR_df_list)
print("Finish phylogenetic signal calculation")

# Select BM and WN phenotypes
BM_phenotype_names_list <- select_BM_traits(binary_AR_df_list, phylo_signal_list, num_phenos)
WN_phenotype_names_list <- select_WN_traits(binary_AR_df_list, phylo_signal_list, num_phenos)
phenotype_names_list <- combine_phenotype_names_lists(BM_phenotype_names_list, WN_phenotype_names_list)
# Select genotypes and phenotypes
genotype_AR_and_conf_mat_list <- select_geno_within_range(binary_AR_df_list,
                                                          binary_conf_mat_list,
                                                          phylo_signal_list,
                                                          lower_bound = -1.5,
                                                          upper_bound = 1.5,
                                                          min_genos = 10)
genotype_AR_mat_list <- genotype_AR_and_conf_mat_list$AR_mat
genotype_conf_mat_list <- genotype_AR_and_conf_mat_list$conf_mat
BM_phenotype_AR_and_conf_mat_list <-
  subsample_to_phenotypes(binary_AR_df_list,
                          binary_conf_mat_list,
                          BM_phenotype_names_list)
BM_phenotype_AR_mat_list <- BM_phenotype_AR_and_conf_mat_list$AR_mat
BM_phenotype_conf_mat_list <- BM_phenotype_AR_and_conf_mat_list$conf_mat
WN_phenotype_AR_and_conf_mat_list <-
  subsample_to_phenotypes(binary_AR_df_list,
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
print("Finish BM transition edges")

# WN pheno
WN_pheno_recon_by_edge_list <- prep_pheno_recon_edges(WN_phenotype_AR_mat_list, tree_list)
WN_pheno_recon_conf_by_edge_list <- reorder_tip_and_node_to_edge_lists(WN_phenotype_conf_mat_list, tree_list)
WN_pheno_sync_trans_by_edge_list <- find_transition_edges(tree_list, WN_phenotype_AR_mat_list, "discrete")
print("Finish WN transition edges")

# Geno
genotype_sync_trans_by_edge_list <- find_transition_edges(tree_list, genotype_AR_mat_list, "discrete")
genotype_phyc_trans_by_edge_list <- convert_to_phyc_trans(genotype_AR_mat_list, genotype_sync_trans_by_edge_list)
genotype_recon_by_edge_list <- reorder_tip_and_node_to_edge_lists(genotype_AR_mat_list, tree_list)
genotype_recon_conf_by_edge_list <- reorder_tip_and_node_to_edge_lists(genotype_conf_mat_list, tree_list)
print("Finish geno transition edges")

# Identify high confidence edges ----
print("start phyc bm hi conf")
phyc_geno_trans_BM_recon_hi_conf_obj_list <-
  prepare_high_confidence_objects_lists(genotype_phyc_trans_by_edge_list,
                                        tree_list,
                                        BM_pheno_recon_conf_by_edge_list,
                                        bootstrap_threshold,
                                        genotype_AR_mat_list,
                                        genotype_recon_conf_by_edge_list,
                                        genotype_recon_by_edge_list,
                                        snps_in_each_gene = NULL)
print("start phyc wn hi conf")
phyc_geno_trans_WN_recon_hi_conf_obj_list <-
  prepare_high_confidence_objects_lists(genotype_phyc_trans_by_edge_list,
                                        tree_list,
                                        WN_pheno_recon_conf_by_edge_list,
                                        bootstrap_threshold,
                                        genotype_AR_mat_list,
                                        genotype_recon_conf_by_edge_list,
                                        genotype_recon_by_edge_list,
                                        snps_in_each_gene = NULL)
print("start sync bm hi conf")
sync_geno_trans_BM_trans_hi_conf_obj_list <-
  prepare_high_confidence_objects_lists(genotype_sync_trans_by_edge_list,
                                        tree_list,
                                        BM_pheno_recon_conf_by_edge_list,
                                        bootstrap_threshold,
                                        genotype_AR_mat_list,
                                        genotype_recon_conf_by_edge_list,
                                        genotype_recon_by_edge_list,
                                        snps_in_each_gene = NULL)
print("start sync wn hi conf")
sync_geno_trans_WM_trans_hi_conf_obj_list <-
  prepare_high_confidence_objects_lists(genotype_sync_trans_by_edge_list,
                                        tree_list,
                                        WN_pheno_recon_conf_by_edge_list,
                                        bootstrap_threshold,
                                        genotype_AR_mat_list,
                                        genotype_recon_conf_by_edge_list,
                                        genotype_recon_by_edge_list,
                                        snps_in_each_gene = NULL)

print("Finish high confidence objects")

# Calculate gamma
BM_phyc_gamma_list <- calc_phyc_gamma_list(tree_list, BM_pheno_recon_by_edge_list, phyc_geno_trans_BM_recon_hi_conf_obj_list)
WN_phyc_gamma_list <- calc_phyc_gamma_list(tree_list, WN_pheno_recon_by_edge_list, phyc_geno_trans_WN_recon_hi_conf_obj_list)
BM_sync_gamma_list <- calc_sync_gamma_list(tree_list, BM_pheno_sync_trans_by_edge_list, sync_geno_trans_BM_trans_hi_conf_obj_list)
WN_sync_gamma_list <- calc_sync_gamma_list(tree_list, WN_pheno_sync_trans_by_edge_list, sync_geno_trans_WM_trans_hi_conf_obj_list)
print("Finish calculating gamma")

# Subset genotype matrix for each tree - phenotype combo
BM_phyc_genotype_keeper_list <- keep_good_genotypes(BM_phyc_gamma_list, bin_size)
BM_sync_genotype_keeper_list <- keep_good_genotypes(BM_sync_gamma_list, bin_size)
WN_phyc_genotype_keeper_list <- keep_good_genotypes(WN_phyc_gamma_list, bin_size)
WN_sync_genotype_keeper_list <- keep_good_genotypes(WN_sync_gamma_list, bin_size)

# Save data
print("Begin to save data")
save_data(tree_list = tree_list,
          genotype_AR_mat_list = genotype_AR_mat_list,
          genotype_phyc_trans_list = genotype_phyc_trans_by_edge_list,
          genotype_sync_trans_list = genotype_sync_trans_by_edge_list,
          BM_phenotype_AR_mat_list = BM_phenotype_AR_mat_list,
          BM_pheno_recon_by_edge_list = BM_pheno_recon_by_edge_list,
          BM_phyc_gamma_list = BM_phyc_gamma_list,
          BM_sync_gamma_list = BM_sync_gamma_list,
          WN_phenotype_AR_mat_list = WN_phenotype_AR_mat_list,
          WN_pheno_recon_by_edge_list = WN_pheno_recon_by_edge_list,
          WN_phyc_gamma_list = WN_phyc_gamma_list,
          WN_sync_gamma_list = WN_sync_gamma_list,
          BM_phyc_genotype_keeper_list = BM_phyc_genotype_keeper_list,
          BM_sync_genotype_keeper_list = BM_sync_genotype_keeper_list,
          WN_phyc_genotype_keeper_list = WN_phyc_genotype_keeper_list,
          WN_sync_genotype_keeper_list = WN_sync_genotype_keeper_list)
print("Finished")
