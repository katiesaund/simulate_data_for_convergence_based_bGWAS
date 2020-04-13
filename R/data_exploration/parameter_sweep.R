# Use every threshold (based on phenotype value) in order to construct genotypes. Find highest possible epsilon. 

library(ape)
library(scales)
library(phytools)
source("R/ancestral_reconstruction.R")
source("R/transition_edges.R")
source("R/tree.R")
source("R/continuous_trait_lib.R")


# LIB
get_pheno_delta_geno_recon <- function(tr, ph, ge) {
  pheno_recon <- ancestral_reconstruction_by_ML(tr, ph, 1, "continuous")
  geno_recon <- ancestral_reconstruction_by_ML(tr, ge, 1, "discrete")
  geno_trans_list <- identify_transition_edges(tr, ge, 1, geno_recon$node_anc_rec, "discrete")
  pheon_recon_edge_mat <- convert_to_edge_mat(tr, pheno_recon$tip_and_node_recon)
  geno_trans_index <- which(geno_trans_list$transition == 1)
  geno_non_trans_index <- which(geno_trans_list$transition == 0)
  trans_pheno_delta_edge <- calculate_phenotype_change_on_edge(geno_trans_index, pheon_recon_edge_mat)
  non_trans_pheno_delta_edge <- calculate_phenotype_change_on_edge(geno_non_trans_index, pheon_recon_edge_mat)
  all_pheno_delta_edge <- calculate_phenotype_change_on_edge(1:Nedge(tr), pheon_recon_edge_mat)
  element_wise_geno_trans_sum <- sum(geno_trans_list$transition * all_pheno_delta_edge)
  element_wise_geno_nontrans_sum <- sum((1 * !geno_trans_list$transition) * all_pheno_delta_edge)
  
  geno_recon_edge_color <- reorder_tip_and_node_to_edge(geno_recon$tip_and_node_recon, tr)
  geno_recon_edge_color[geno_recon_edge_color == 1] <- "red"
  geno_recon_edge_color[geno_recon_edge_color == 0] <- "black"
  
  geno_trans_edge_color <- geno_trans_list$transition
  geno_trans_edge_color[geno_trans_edge_color == 1] <- "cyan"
  geno_trans_edge_color[geno_trans_edge_color == 0] <- "black"
  return(list("geno_trans_vec" = geno_trans_list, 
              "delta_pheno_vec" = all_pheno_delta_edge, 
              "trans_pheno_delta_edge" = trans_pheno_delta_edge, 
              "non_trans_pheno_delta_edge" = non_trans_pheno_delta_edge, 
              "geno_trans_index" = geno_trans_index, 
              "geno_non_trans_index" = geno_non_trans_index, 
              "element_wise_geno_trans_sum" = element_wise_geno_trans_sum, 
              "element_wise_geno_nontrans_sum" = element_wise_geno_nontrans_sum,
              "pheno_node_recon" = pheno_recon$node_anc_rec, 
              "geno_trans_edge_color" = geno_trans_edge_color, 
              "geno_recon_edge_color" = geno_recon_edge_color))
}


make_WN_phenotype <- function(pheno_mat) {
  pheno_vec <- pheno_mat[, 1]
  lambda_has_high_signal <- TRUE
  while (lambda_has_high_signal) {
    jumbled_pheno <- sample(unname(pheno_vec), size = length(pheno_vec), replace = FALSE)
    jumbled_lambda <- phytools::phylosig(tree = treeA,
                                         x = jumbled_pheno,
                                         method = "lambda")
    lambda_has_high_signal <- jumbled_lambda$lambda < -0.05 & jumbled_lambda$lambda > 0.05
  }
  pheno_mat[, 1] <- jumbled_pheno
  return(pheno_mat)
}

epsilon_b_scaled_to_one <- function(geno_transition_vec, delta_pheno_vec) {
  # intersection / union, pheno scaled to the number of genotype transition edges
  sum_A <- sum(geno_transition_vec)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(0, 1))
  scaled_intersection <- sum(geno_transition_vec * scaled_pheno)
  sum_B <- sum(scaled_pheno)
  union <- sum_A + sum_B - scaled_intersection
  return(scaled_intersection / union)
}

plot_geno_non_and_trans_hist <- function(ge_tr, ge_no_tr, title) {
  max_x <- max(c(ge_tr, ge_no_tr))
  min_x <- min(c(ge_tr, ge_no_tr))
  # break_seq <- seq(from = min_x, to = max_x, by = 1/20)
  
  par(mar = c(2, 2, 2, 2))
  hist(ge_no_tr, col = rgb(0, 0, 0, 0.5), xlim = c(min_x, max_x) ,main = title)
  hist(ge_tr, col = rgb(1, 0, 0, 0.5), add = TRUE)
}

# Set up
set.seed(1865)
treeA <- ape::rcoal(100)
treeA <- midpoint.root(treeA)
genoA <- matrix(NA, nrow = ape::Ntip(treeA), ncol = 1)
genoA[, 1] <- c(1, 1, 0, 0, 1,
                1, 1, 1, 1, 1,
                0, 0, 1, 1, 0,
                1, 1, 1, 1, 0)
row.names(genoA) <- treeA$tip.label 
phenoA <- genoA
set.seed(1)
phenoA[, 1] <- ape::rTraitCont(treeA, model = "BM")

WN <- FALSE
if (WN) {
  phenoA <- make_WN_phenotype(phenoA)
}
dotTree(treeA, phenoA, length = 10, ftype = "i")

# Create genotype matrix with 1/0 determined by phenotype rank
geno_mat <- matrix(rank(phenoA), nrow = Ntip(treeA), ncol = 2 * Ntip(treeA))
row.names(geno_mat) <- treeA$tip.label
colnames(geno_mat) <- c(paste0("pheno_rank_below_", 1:Ntip(treeA)), paste0("pheno_rank_above_", 1:Ntip(treeA)))

for (i in 1:Ntip(treeA)) {
  geno_mat[, i] <- as.numeric(geno_mat[, i] < i)
}

for (i in (Ntip(treeA) + 1):(2 * Ntip(treeA))) {
  geno_mat[, i] <- as.numeric(geno_mat[, i] > (i - Ntip(treeA)))
}

geno_mat <- geno_mat[, colSums(geno_mat) > 2] # remove genotypes where ancestral reconstruction doesn't work well

# Calculations
geno_out <- list()

for (i in 1:ncol(geno_mat)) {
  temp_geno <- geno_mat[, i, drop = FALSE]
  geno_out[[i]] <- get_pheno_delta_geno_recon(treeA, phenoA, temp_geno)
}



# Make a table of epsilon score
num_data_sets <- ncol(geno_mat)
epsilon_mat <- matrix(0, nrow = num_data_sets, ncol = 1)
colnames(epsilon_mat) <- c("b_scaled_to_one")
row.names(epsilon_mat) <- colnames(geno_mat)

for (i in 1:ncol(geno_mat)) {
  epsilon_mat[i, 1] <- epsilon_b_scaled_to_one(geno_out[[i]]$geno_trans_vec$transition, geno_out[[i]]$delta_pheno_vec) 
}

hist(epsilon_mat[, 1])
  
# Which ones had the highest epsilon?
row.names(epsilon_mat[epsilon_mat[,1] > 0.12, , drop = FALSE])

# Function to generate a vector of colors
whiteToRed = colorRampPalette(c("white", "firebrick"))
# Generate a character vector of colors with each shade of red 
heatMapCols = whiteToRed(100)
pheatmap::pheatmap(epsilon_mat, 
                   color = heatMapCols, 
                   main = "epsilons",
                   cluster_cols = FALSE,
                   cluster_rows = FALSE, 
                   scale = "none", 
                   show_rownames = TRUE)

# The intermediate values below 32 - 57, above 31-56 -- but skipping a lot of them.