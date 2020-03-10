library(ape)
library(scales)
library(phytools)
source("R/ancestral_reconstruction.R")
source("R/transition_edges.R")
source("R/tree.R")
source("R/continuous_trait_lib.R")


# Tree A
set.seed(1)
treeA <- ape::rcoal(20)
treeA <- midpoint.root(treeA)
phenoA <- matrix(NA, nrow = ape::Ntip(treeA), ncol = 1)
phenoA[, 1] <- c(1, 0, 0, 0, 0, 
                 0, 0, 0, 0, 1, 
                 0, 0, 0, 0, 1, 
                 0, 0, 0, 0, 0)
genoA <- phenoA
row.names(phenoA) <- treeA$tip.label 
dotTree(treeA, phenoA, length = 10, ftype = "i")


# Tree B
set.seed(1865)
treeB <- ape::rcoal(20)
treeB <- midpoint.root(treeB)
phenoB <- matrix(NA, nrow = ape::Ntip(treeB), ncol = 1)
phenoB[, 1] <- c(1, 1, 0, 0, 1,
                 1, 1, 1, 1, 1,
                 0, 0, 1, 1, 0,
                 1, 1, 1, 1, 0)
row.names(phenoB) <- treeB$tip.label 
genoB <- phenoB
dotTree(treeB, phenoB, length = 10, ftype = "i")

# Tree C
set.seed(2)
treeC <- ape::rcoal(20)
treeC <- midpoint.root(treeC)
phenoC <- matrix(NA, nrow = ape::Ntip(treeC), ncol = 1)
phenoC[, 1] <- c(1, 1, 1, 1, 0,
                 0, 0, 1, 1, 1, 
                 1, 1, 1, 1, 1, 
                 0, 0, 1, 1, 1)
row.names(phenoC) <- treeC$tip.label 
genoC <- phenoC
dotTree(treeC, phenoC, length = 10, ftype = "i")

# Tree D
set.seed(10)
treeD <- ape::rcoal(20)
treeD <- midpoint.root(treeD)
phenoD <- matrix(NA, nrow = ape::Ntip(treeD), ncol = 1)
phenoD[, 1] <- c(1, 1, 1, 1, 0,
                 0, 0, 1, 1, 1, 
                 1, 1, 1, 1, 1, 
                 0, 0, 1, 1, 1)
row.names(phenoD) <- treeD$tip.label 
genoD <- phenoD

phenoD[, 1] <- c(0.5, 0.75, 0.8, 1, 0.01,
                 0.02, 0.05, 1, 0.95, 0.9, 
                 1, 0.90, 0.97, 0.98, 0.96, 
                 0.1, 0.05, 0.8, .75, 0.85)
dotTree(treeD, phenoD, length = 10, ftype = "i")
# LIB
epsilon_b_scaled_to_num_geno_trans_edges <- function(geno_transition_vec, delta_pheno_vec) {
  # intersection / union, pheno scaled to the number of genotype transition edges
  sum_A <- sum(geno_transition_vec)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(0, sum_A))
  scaled_intersection <- sum(geno_transition_vec * scaled_pheno)
  sum_B <- sum(scaled_pheno)
  union <- sum_A + sum_B - scaled_intersection
  return(scaled_intersection / union)
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

plot_geno_non_and_trans_hist <- function(ge_tr, ge_no_tr) {
  max_x <- max(c(ge_tr, ge_no_tr))
  min_x <- min(c(ge_tr, ge_no_tr))
  par(mar = c(2, 2, 2, 2))
  hist(ge_no_tr, col = rgb(0, 0, 0, 0.5), xlim = c(min_x, max_x))
  hist(ge_tr, col = rgb(1, 0, 0, 0.5), add = TRUE)
}

# Calculations
# A 
a_out <- get_pheno_delta_geno_recon(treeA, phenoA, genoA)
par(mfrow = c(2, 2))
dotTree(treeA,  phenoA, length = 10, ftype = "i")
edgelabels(round(a_out$delta_pheno_vec, 1))
graphics::plot(treeA,
               edge.width = 5,
               edge.color = a_out$geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
graphics::plot(treeA,
               edge.width = 5,
               edge.color = a_out$geno_trans_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
par(mfrow = c(1, 1))
plot_geno_non_and_trans_hist(a_out$trans_pheno_delta_edge, a_out$non_trans_pheno_delta_edge)
phytools::contMap(treeA,
                  phenoA[, 1, drop = TRUE],
                  method = "user",
                  anc.states = a_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(a_out$delta_pheno_vec, 1))
epsilon_b_scaled_to_one(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec) # 0.318
epsilon_b_scaled_to_num_geno_trans_edges(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec) # 0.412
epsilon_using_median_values(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec) # 0.130

# B 
b_out <- get_pheno_delta_geno_recon(treeB, phenoB, genoB)
par(mfrow = c(2, 2))
dotTree(treeB,  phenoB, length = 10, ftype = "i")
edgelabels(round(b_out$delta_pheno_vec, 1))
graphics::plot(treeB,
               edge.width = 5,
               edge.color = b_out$geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
graphics::plot(treeB,
               edge.width = 5,
               edge.color = b_out$geno_trans_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
par(mfrow = c(1, 1))
plot_geno_non_and_trans_hist(b_out$trans_pheno_delta_edge, b_out$non_trans_pheno_delta_edge)
phytools::contMap(treeB,
                  phenoB[, 1, drop = TRUE],
                  method = "user",
                  anc.states = b_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(b_out$delta_pheno_vec, 1))
epsilon_b_scaled_to_one(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec) # 0.351
epsilon_b_scaled_to_num_geno_trans_edges(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec) # 0.494
epsilon_using_median_values(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec) # 0.167

# C
c_out <- get_pheno_delta_geno_recon(treeC, phenoC, genoC)
par(mfrow = c(2, 2))
dotTree(treeC,  phenoC, length = 10, ftype = "i")
edgelabels(round(c_out$delta_pheno_vec, 1))
graphics::plot(treeC,
               edge.width = 5,
               edge.color = c_out$geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
graphics::plot(treeC,
               edge.width = 5,
               edge.color = c_out$geno_trans_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
par(mfrow = c(1, 1))
plot_geno_non_and_trans_hist(c_out$trans_pheno_delta_edge, c_out$non_trans_pheno_delta_edge)
phytools::contMap(treeC,
                  phenoC[, 1, drop = TRUE],
                  method = "user",
                  anc.states = c_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(c_out$delta_pheno_vec, 1))
epsilon_b_scaled_to_one(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec) # 0.332
epsilon_b_scaled_to_num_geno_trans_edges(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec) # 0.544
epsilon_using_median_values(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec) # 0.23

# D
d_out <- get_pheno_delta_geno_recon(treeD, phenoD, genoD)
par(mfrow = c(2, 2))
dotTree(treeD,  phenoD, length = 10, ftype = "i")
edgelabels(round(d_out$delta_pheno_vec, 1))
graphics::plot(treeD,
               edge.width = 5,
               edge.color = d_out$geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
graphics::plot(treeD,
               edge.width = 5,
               edge.color = d_out$geno_trans_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0)
par(mfrow = c(1, 1))
plot_geno_non_and_trans_hist(d_out$trans_pheno_delta_edge, d_out$non_trans_pheno_delta_edge)
phytools::contMap(treeD,
                  phenoD[, 1, drop = TRUE],
                  method = "user",
                  anc.states = d_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(d_out$delta_pheno_vec, 1))
epsilon_b_scaled_to_one(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec) # 0.239
epsilon_b_scaled_to_num_geno_trans_edges(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec) # 0.400
epsilon_using_median_values(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec) # 0

# # let's dig into the continuous results ancestral reconstruction confidence issue
# set.seed(1)
# reconstruction <- ape::ace(pheno[, 1, drop = TRUE],
#                            tree,
#                            model = "BM",
#                            type = "continuous",
#                            method = "REML", # pic?
#                            marginal = FALSE)
