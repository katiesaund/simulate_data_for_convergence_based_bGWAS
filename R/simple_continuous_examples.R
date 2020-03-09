library(ape)
library(scales)

source("R/ancestral_reconstruction.R")
source("R/transition_edges.R")
source("R/tree.R")
source("R/continuous_trait_lib.R")
set.seed(1865)
tree <- ape::rcoal(20)
plot(tree)

pheno <- matrix(NA, nrow = ape::Ntip(tree), ncol = 1)
row.names(pheno) <- tree$tip.label
pheno[, 1] <- c(100, 0, 0, 100, 100, 0, 100, 100, 100, 100, 0, 100, 100, 0, 0, 100, 100, 100, 100, 0)

geno <- pheno
geno[geno == 100, ] <- 1


pheno_recon <- ancestral_reconstruction_by_ML(tree, pheno, 1, "continuous")
geno_recon <- ancestral_reconstruction_by_ML(tree, geno, 1, "discrete")

geno_trans_list <- identify_transition_edges(tree, geno, 1, geno_recon$node_anc_rec, "discrete")
  
pheon_recon_edge_mat <- convert_to_edge_mat(tree, pheno_recon$tip_and_node_recon)
  
geno_trans_index <- which(geno_trans_list$transition == 1)
geno_non_trans_index <- which(geno_trans_list$transition == 0)
trans_pheno_delta_edge <- calculate_phenotype_change_on_edge(geno_trans_index, pheon_recon_edge_mat)
non_trans_pheno_delta_edge <- calculate_phenotype_change_on_edge(geno_non_trans_index, pheon_recon_edge_mat)
all_pheno_delta_edge <- calculate_phenotype_change_on_edge(1:Nedge(tree), pheon_recon_edge_mat)

element_wise_geno_trans_sum <- sum(geno_trans_list$transition * all_pheno_delta_edge)
element_wise_geno_nontrans_sum <- sum((1 * !geno_trans_list$transition) * all_pheno_delta_edge)

hist(non_trans_pheno_delta_edge, col = "white", xlim = c(0, 80), breaks = 30)
hist(trans_pheno_delta_edge, col = "red", add = TRUE)

epsilon_b_scaled_to_num_geno_trans_edges <- function(geno_transition_vec, delta_pheno_vec) {
  # intersection / union, pheno scaled to the number of genotype transition edges
  sum_A <- sum(geno_transition_vec)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(0, sum_A))
  scaled_intersection <- sum(geno_transition_vec * scaled_pheno)
  sum_B <- sum(scaled_pheno)
  union <- sum_A + sum_B - scaled_intersection
  return(scaled_intersection / union)
}

epsilon_b_scaled_to_num_geno_trans_edges(geno_trans_list$transition, all_pheno_delta_edge) # 0.572

epsilon_b_scaled_to_one <- function(geno_transition_vec, delta_pheno_vec) {
  # intersection / union, pheno scaled to the number of genotype transition edges
  sum_A <- sum(geno_transition_vec)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(0, 1))
  scaled_intersection <- sum(geno_transition_vec * scaled_pheno)
  sum_B <- sum(scaled_pheno)
  union <- sum_A + sum_B - scaled_intersection
  return(scaled_intersection / union)
}

epsilon_b_scaled_to_one(geno_trans_list$transition, all_pheno_delta_edge) # 0.379

epsilon_using_median_values <- function(geno_transition_vec, delta_pheno_vec) {
  # Use the median of the nontrans edges
  non_trans_median_delta_pheno <- median(delta_pheno_vec[!geno_transition_vec])
  numerator <- sum(delta_pheno_vec[geno_transition_vec] > non_trans_median_delta_pheno)
  num_delta_pheno_greater_than_median <- sum(delta_pheno_vec > non_trans_median_delta_pheno)
  denominator <- sum(geno_transition_vec, num_delta_pheno_greater_than_median)
  print(numerator)
  print(denominator)
  return(numerator / denominator)
}
epsilon_using_median_values(geno_trans_list$transition, all_pheno_delta_edge) # 0.192

