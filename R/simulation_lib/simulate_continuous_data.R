suppressWarnings(library(ape))
suppressWarnings(library(caper))
suppressWarnings(library(phytools))
suppressWarnings(library(scales))

# Initialize variables / read in user input ----
args <- commandArgs(trailingOnly = TRUE)
num_trees <- as.numeric(args[1]) 
num_phenos <- as.numeric(args[2]) 
num_tips <- as.numeric(args[3]) 
num_start_trait <- as.numeric(args[4])
rel_path <- args[5]

source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/tree.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/discrete_trait_lib.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/continuous_trait_lib.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/transition_edges.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/convergence.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/save_data.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/keep_interesting_genotypes.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/ancestral_reconstruction.R"))
source(paste0(rel_path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/high_confidence.R"))


# Defaults that shouldn't need user input ----
tree_edge_multiplier <- 100
bootstrap_threshold <- 70
bin_size <- 20

# Generate huge matrix of binary traits specific to trees ----
tree_list <- generate_trees(num_trees, num_tips, tree_edge_multiplier)
print("Finish making trees")

print("Start mkaing phenotypes")
# Make BM and WN continuous phenotypes and then ancestral reconstructions
cont_pheno_mat_list <- make_continuous_phenotypes(tree_list, num_phenos) # cont_pheno_WN_mat_list[[tree index]] (matrix, each column a phenotype)
cont_pheno_BM_mat_list <- cont_pheno_mat_list$cont_pheno_BM_mat_list
cont_pheno_WN_mat_list <- cont_pheno_mat_list$cont_pheno_WN_mat_list
# Ancestral Reconstruction of continous phenotypes
cont_pheno_BM_AR_and_conf_mat_list <- 
  ancestral_reconstruction(cont_pheno_BM_mat_list, 
                           tree_list, 
                           "continuous")

cont_pheno_WN_AR_and_conf_mat_list <- 
  ancestral_reconstruction(cont_pheno_WN_mat_list, 
                           tree_list, 
                           "continuous")

# Generated a mat_list, where each [[tree_num]] and each column is a phenotype. 
# It's ordered by tips and then nodes (rows).
# Separate traits into BM and WN

BM_phenotype_AR_mat_list <- cont_pheno_BM_AR_and_conf_mat_list$AR_mat
BM_phenotype_conf_mat_list <- cont_pheno_BM_AR_and_conf_mat_list$conf_mat
BM_phenotype_recon_edge_mat_list <- 
  cont_pheno_BM_AR_and_conf_mat_list$recon_edge_mat

WN_phenotype_AR_mat_list <- cont_pheno_WN_AR_and_conf_mat_list$AR_mat
WN_phenotype_conf_mat_list <- cont_pheno_WN_AR_and_conf_mat_list$conf_mat
WN_phenotype_recon_edge_mat_list <- 
  cont_pheno_WN_AR_and_conf_mat_list$recon_edge_mat
print("Finish pheno reconstructions")

# Order everything by edges instead of by tips then nodes ----
# BM pheno
BM_pheno_recon_by_edge_list <- prep_pheno_recon_edges(BM_phenotype_AR_mat_list,
                                                      tree_list)
BM_pheno_recon_conf_by_edge_list <- 
  reorder_tip_and_node_to_edge_lists(BM_phenotype_conf_mat_list, tree_list)
print("Finish BM edge ordering")

# WN pheno
WN_pheno_recon_by_edge_list <- 
  prep_pheno_recon_edges(WN_phenotype_AR_mat_list, tree_list)
WN_pheno_recon_conf_by_edge_list <-
  reorder_tip_and_node_to_edge_lists(WN_phenotype_conf_mat_list, tree_list)
print("Finish WN edges ordering")

print("Start making genotypes")
binary_AR_df_list <- generate_binary_df_list(tree_list, num_start_trait)
print("Finish generating first draft discrete matrix")

binary_AR_df_list <- add_geno_continuous(binary_AR_df_list, 
                                         tree_list, 
                                         num_phenos,
                                         cont_pheno_BM_mat_list,
                                         cont_pheno_WN_mat_list, 
                                         BM_phenotype_recon_edge_mat_list, 
                                         WN_phenotype_recon_edge_mat_list) 
# Due to new genotypes the ancestral reconstructions are now wrong and need to be recalculated in the next step.
print("Finish adding additional genotypes")

binary_AR_and_conf_mat <- 
  ancestral_reconstruction(binary_AR_df_list, tree_list, "discrete") # So fix the ancestral reconstructions here
binary_conf_mat_list <- binary_AR_and_conf_mat$conf_mat
binary_AR_df_list <- binary_AR_and_conf_mat$AR_mat
print("Finished genotype ancestral reconstruction")

phylo_signal_list <- calculate_phylo_signal(tree_list, binary_AR_df_list)
print("Finish phylogenetic signal calculation of genotypes")

# Select genotypes
genotype_AR_and_conf_mat_list <- select_geno_within_range(binary_AR_df_list,
                                                          binary_conf_mat_list,
                                                          phylo_signal_list,
                                                          lower_bound = -1.5,
                                                          upper_bound = 1.5,
                                                          min_genos = 10)
genotype_AR_mat_list <- genotype_AR_and_conf_mat_list$AR_mat
genotype_conf_mat_list <- genotype_AR_and_conf_mat_list$conf_mat
print("Finish subsetting genotypes to those with -1.5 < D < 1.5")

# Geno
genotype_cont_trans_by_edge_list <- 
  find_transition_edges(tree_list, genotype_AR_mat_list, "discrete")
genotype_recon_by_edge_list <- 
  reorder_tip_and_node_to_edge_lists(genotype_AR_mat_list, tree_list)
genotype_recon_conf_by_edge_list <- 
  reorder_tip_and_node_to_edge_lists(genotype_conf_mat_list, tree_list)
print("Finish identifying genotype transition edges")

# Identify high confidence edges ----
print("start cont bm hi conf")
cont_geno_trans_BM_trans_hi_conf_obj_list <-
  prepare_high_confidence_objects_lists(genotype_cont_trans_by_edge_list,
                                        tree_list,
                                        BM_pheno_recon_conf_by_edge_list,
                                        bootstrap_threshold,
                                        genotype_AR_mat_list,
                                        genotype_recon_conf_by_edge_list,
                                        genotype_recon_by_edge_list,
                                        snps_in_each_gene = NULL)
print("start cont wn hi conf")
cont_geno_trans_WM_trans_hi_conf_obj_list <-
  prepare_high_confidence_objects_lists(genotype_cont_trans_by_edge_list,
                                        tree_list,
                                        WN_pheno_recon_conf_by_edge_list,
                                        bootstrap_threshold,
                                        genotype_AR_mat_list,
                                        genotype_recon_conf_by_edge_list,
                                        genotype_recon_by_edge_list,
                                        snps_in_each_gene = NULL)

print("Finish high confidence objects")

# Calculate convergence

BM_cont_convergence_list <-
  calc_cont_convergence_list(tree_list, 
                       BM_phenotype_recon_edge_mat_list,
                       cont_geno_trans_BM_trans_hi_conf_obj_list)
WN_cont_convergence_list <- 
  calc_cont_convergence_list(tree_list, 
                       WN_phenotype_recon_edge_mat_list, 
                       cont_geno_trans_WM_trans_hi_conf_obj_list)
print("Finish calculating convergence")

# Subset genotype matrix for each tree - phenotype combo
BM_cont_genotype_keeper_list <- keep_good_genotypes(BM_cont_convergence_list, 
                                                    bin_size)
WN_cont_genotype_keeper_list <- keep_good_genotypes(WN_cont_convergence_list, 
                                                    bin_size)

# Save data
print("Begin to save data")
save_continuous_data(
  tree_list = tree_list,
  genotype_AR_mat_list = genotype_AR_mat_list,
  genotype_cont_trans_list = genotype_cont_trans_by_edge_list,
  BM_phenotype_AR_mat_list = BM_phenotype_AR_mat_list,
  BM_pheno_recon_by_edge_list = BM_pheno_recon_by_edge_list,
  BM_cont_convergence_list = BM_cont_convergence_list,
  WN_phenotype_AR_mat_list = WN_phenotype_AR_mat_list,
  WN_pheno_recon_by_edge_list = WN_pheno_recon_by_edge_list,
  WN_cont_convergence_list = WN_cont_convergence_list,
  BM_cont_genotype_keeper_list = BM_cont_genotype_keeper_list,
  WN_cont_genotype_keeper_list = WN_cont_genotype_keeper_list)
print("Finished")
