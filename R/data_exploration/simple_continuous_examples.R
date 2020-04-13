library(ape)
library(scales)
library(phytools)
source("R/simulation_lib/ancestral_reconstruction.R")
source("R/simulation_lib/transition_edges.R")
source("R/simulation_lib/tree.R")
source("R/simulation_lib/continuous_trait_lib.R")

# Tree A
set.seed(1865)
treeA <- ape::rcoal(20)
treeA <- midpoint.root(treeA)
phenoA <- matrix(NA, nrow = ape::Ntip(treeA), ncol = 1)
phenoA[, 1] <- c(1, 1, 0, 0, 1,
                 1, 1, 1, 1, 1,
                 0, 0, 1, 1, 0,
                 1, 1, 1, 1, 0)
row.names(phenoA) <- treeA$tip.label 
genoA <- phenoA
phenoA <- 100 * phenoA
dotTree(treeA, phenoA, length = 10, ftype = "i")

# Tree B
set.seed(1)
treeB <- ape::rcoal(20)
treeB <- midpoint.root(treeB)
phenoB <- matrix(NA, nrow = ape::Ntip(treeB), ncol = 1)
phenoB[, 1] <- c(1, 0, 0, 0, 0, 
                 0, 0, 0, 0, 1, 
                 0, 0, 0, 0, 1, 
                 0, 0, 0, 0, 0)
genoB <- phenoB
row.names(phenoB) <- treeB$tip.label 
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

epsilon_twice_scaled_to_one <- function(geno_transition_vec, delta_pheno_vec) {
  # intersection / union, pheno scaled to the number of genotype transition edges
  sum_A <- sum(geno_transition_vec)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(0, 1))
  scaled_intersection <- sum(geno_transition_vec * scaled_pheno)
  sum_B <- sum(scaled_pheno)
  denom <- sum_A + sum_B 
  return(2 * (scaled_intersection / denom))
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

epsilon_bin_areas <- function(geno_transition_vec, delta_pheno_vec) {
  num_bins <- 20
  min <- 0
  max <- 1
  zero_to_one_bins <- seq(from = min, to = max, by = 1 / num_bins)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(min, max))
  
  geno_trans_bin_list <- geno_non_trans_bin_list <- 
    intersection <- union <- rep(0, num_bins)
  for (i in 1:num_bins) {
    geno_trans_bin_list[i] <- 
      sum(zero_to_one_bins[i] <= scaled_pheno & 
            scaled_pheno <= zero_to_one_bins[i + 1] &
            geno_transition_vec == 1)
    geno_non_trans_bin_list[i] <- 
      sum(zero_to_one_bins[i] <= scaled_pheno &
            scaled_pheno <= zero_to_one_bins[i + 1] & 
            geno_transition_vec == 0)
    
    intersection[i] <-
      min(geno_trans_bin_list[i], geno_non_trans_bin_list[i])
    union[i] <- 
      geno_trans_bin_list[i] + geno_non_trans_bin_list[i] - intersection[i]
  }
  return(1 - (sum(intersection) / sum(union)))
}

epsilon_bin_use_med <- function(geno_transition_vec, delta_pheno_vec) {
  num_bins <- 20
  min <- 0
  max <- 1
  zero_to_one_bins <- seq(from = min, to = max, by = 1 / num_bins)
  scaled_pheno <- rescale(delta_pheno_vec, to = c(min, max))
  
  geno_trans_bin_list <- geno_non_trans_bin_list <- 
    intersection <- union <- rep(0, num_bins)
  for (i in 1:num_bins) {
    geno_trans_bin_list[i] <- 
      sum(zero_to_one_bins[i] <= scaled_pheno & 
            scaled_pheno <= zero_to_one_bins[i + 1] &
            geno_transition_vec == 1)
    geno_non_trans_bin_list[i] <- 
      sum(zero_to_one_bins[i] <= scaled_pheno &
            scaled_pheno <= zero_to_one_bins[i + 1] & 
            geno_transition_vec == 0)
    
    intersection[i] <-
      min(geno_trans_bin_list[i], geno_non_trans_bin_list[i])
    union[i] <- 
      geno_trans_bin_list[i] + geno_non_trans_bin_list[i] - intersection[i]
  }
  
  non_trans_med_delta <- median(delta_pheno_vec[geno_transition_vec == 0])
  trans_med_delta <- median(delta_pheno_vec[geno_transition_vec == 1])
  
  output <- 1 - (sum(intersection) / sum(union))
  if (non_trans_med_delta > trans_med_delta) {
    output <- (sum(intersection) / sum(union))
  }
  return(output)
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

plot_geno_non_and_trans_hist <- function(ge_tr, ge_no_tr, title) {
  max_x <- max(c(ge_tr, ge_no_tr))
  min_x <- min(c(ge_tr, ge_no_tr))
  # break_seq <- seq(from = min_x, to = max_x, by = 1/20)
  
  par(mar = c(2, 2, 2, 2))
  hist(ge_no_tr, col = rgb(0, 0, 0, 0.5), xlim = c(min_x, max_x) ,main = title)
  hist(ge_tr, col = rgb(1, 0, 0, 0.5), add = TRUE)
}

# Calculations
a_out <- get_pheno_delta_geno_recon(treeA, phenoA, genoA)
b_out <- get_pheno_delta_geno_recon(treeB, phenoB, genoB)
c_out <- get_pheno_delta_geno_recon(treeC, phenoC, genoC)
d_out <- get_pheno_delta_geno_recon(treeD, phenoD, genoD)

# Plots
# A 
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
plot_geno_non_and_trans_hist(a_out$trans_pheno_delta_edge, a_out$non_trans_pheno_delta_edge, "A")
phytools::contMap(treeA,
                  phenoA[, 1, drop = TRUE],
                  method = "user",
                  anc.states = a_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(a_out$delta_pheno_vec, 1))
ep_one_a <- epsilon_b_scaled_to_one(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec) # 0.318
ep_num_a <- epsilon_b_scaled_to_num_geno_trans_edges(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec) # 0.412
ep_nt_med_a <- epsilon_using_median_values(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec) # 0.130
ep_bin_a <- epsilon_bin_areas(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec)
ep_bin_med_a <- epsilon_bin_use_med(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec)
ep_twice_a <- epsilon_twice_scaled_to_one(a_out$geno_trans_vec$transition, a_out$delta_pheno_vec)
# B 
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
plot_geno_non_and_trans_hist(b_out$trans_pheno_delta_edge, b_out$non_trans_pheno_delta_edge, "B")
phytools::contMap(treeB,
                  phenoB[, 1, drop = TRUE],
                  method = "user",
                  anc.states = b_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(b_out$delta_pheno_vec, 1))
ep_one_b <- epsilon_b_scaled_to_one(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec) # 0.351
ep_num_b <- epsilon_b_scaled_to_num_geno_trans_edges(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec) # 0.494
ep_nt_med_b <- epsilon_using_median_values(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec) # 0.167
ep_bin_b <- epsilon_bin_areas(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec)
ep_bin_med_b <- epsilon_bin_use_med(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec)
ep_twice_b <- epsilon_twice_scaled_to_one(b_out$geno_trans_vec$transition, b_out$delta_pheno_vec)

# C

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
plot_geno_non_and_trans_hist(c_out$trans_pheno_delta_edge, c_out$non_trans_pheno_delta_edge ,"C")
phytools::contMap(treeC,
                  phenoC[, 1, drop = TRUE],
                  method = "user",
                  anc.states = c_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(c_out$delta_pheno_vec, 1))
ep_one_c <- epsilon_b_scaled_to_one(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec) # 0.332
ep_num_c <- epsilon_b_scaled_to_num_geno_trans_edges(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec) # 0.544
ep_nt_med_c <- epsilon_using_median_values(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec) # 0.23
ep_bin_c <- epsilon_bin_areas(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec)
ep_bin_med_c <- epsilon_bin_use_med(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec)
ep_twice_c <- epsilon_twice_scaled_to_one(c_out$geno_trans_vec$transition, c_out$delta_pheno_vec)

# D
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
plot_geno_non_and_trans_hist(d_out$trans_pheno_delta_edge, d_out$non_trans_pheno_delta_edge, "D")
phytools::contMap(treeD,
                  phenoD[, 1, drop = TRUE],
                  method = "user",
                  anc.states = d_out$pheno_node_recon,
                  plot = TRUE)
edgelabels(round(d_out$delta_pheno_vec, 1))
ep_one_d <- epsilon_b_scaled_to_one(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec) # 0.239
ep_num_d <- epsilon_b_scaled_to_num_geno_trans_edges(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec) # 0.400
ep_nt_med_d <- epsilon_using_median_values(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec) # 0
ep_bin_d <- epsilon_bin_areas(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec)
ep_bin_med_d <- epsilon_bin_use_med(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec)
ep_twice_d <- epsilon_twice_scaled_to_one(d_out$geno_trans_vec$transition, d_out$delta_pheno_vec)

# Fake E
e_geno_trans_vec <- c(rep(1, 10), rep(0, Nedge(treeA) - 10))
e_delta_pheno_vec <- c(seq(from = 0.0, to = 0.8, by = 0.08), 
                       seq(from = 0.0, to = 0.8, by = 0.03))
e_trans_delta <- e_delta_pheno_vec[e_geno_trans_vec == 1]
e_non_trans_delta <- e_delta_pheno_vec[e_geno_trans_vec == 0]
ep_one_e <- epsilon_b_scaled_to_one(e_geno_trans_vec, e_delta_pheno_vec) # 0.186
ep_num_e <- epsilon_b_scaled_to_num_geno_trans_edges(e_geno_trans_vec, e_delta_pheno_vec) # 0.297
ep_nt_med_e <- epsilon_using_median_values(e_geno_trans_vec, e_delta_pheno_vec) # 0
ep_bin_e <- epsilon_bin_areas(e_geno_trans_vec, e_delta_pheno_vec) # 0.571
ep_bin_med_e <- epsilon_bin_use_med(e_geno_trans_vec, e_delta_pheno_vec)
ep_twice_e <- epsilon_twice_scaled_to_one(e_geno_trans_vec, e_delta_pheno_vec)


# FAKE F
f_geno_trans_vec <- c(rep(1, 10), rep(0, Nedge(treeA) - 10))
f_delta_pheno_vec <- c(seq(from = 0.0, to = 0.4, by = 0.04), 
                       seq(from = 0.4, to = 0.8, by = 0.015))
f_trans_delta <- f_delta_pheno_vec[f_geno_trans_vec == 1]
f_non_trans_delta <- f_delta_pheno_vec[f_geno_trans_vec == 0]
ep_one_f <- epsilon_b_scaled_to_one(f_geno_trans_vec, f_delta_pheno_vec) # 0.074
ep_num_f <- epsilon_b_scaled_to_num_geno_trans_edges(f_geno_trans_vec, f_delta_pheno_vec) # 0.104
ep_nt_med_f <- epsilon_using_median_values(f_geno_trans_vec, f_delta_pheno_vec) # 0
ep_bin_f <- epsilon_bin_areas(f_geno_trans_vec, f_delta_pheno_vec) # 1
ep_bin_med_f <- epsilon_bin_use_med(f_geno_trans_vec, f_delta_pheno_vec)
ep_twice_f <- epsilon_twice_scaled_to_one(f_geno_trans_vec, f_delta_pheno_vec)


# FAKE AA
aa_geno_trans_vec <- c(rep(1, 9), rep(0, Nedge(treeA) - 9))
aa_delta_pheno_vec <- c(seq(from = 0.6, to = 1.0, by = 0.05), 
                       seq(from = 0.00, to = 0.4, by = 0.014))
aa_trans_delta <- aa_delta_pheno_vec[aa_geno_trans_vec == 1]
aa_non_trans_delta <- aa_delta_pheno_vec[aa_geno_trans_vec == 0]
ep_one_aa <- epsilon_b_scaled_to_one(aa_geno_trans_vec, aa_delta_pheno_vec) # 0.490
ep_num_aa <- epsilon_b_scaled_to_num_geno_trans_edges(aa_geno_trans_vec, aa_delta_pheno_vec) # 1.077
ep_nt_med_aa <- epsilon_using_median_values(aa_geno_trans_vec, aa_delta_pheno_vec) # 0.281
ep_bin_aa <- epsilon_bin_areas(aa_geno_trans_vec, aa_delta_pheno_vec) # 1.0
ep_bin_med_aa <- epsilon_bin_use_med(aa_geno_trans_vec, aa_delta_pheno_vec) # 1.0
ep_twice_aa <- epsilon_twice_scaled_to_one(aa_geno_trans_vec, aa_delta_pheno_vec)

# Plotted from what should be highest epsilon to lowest epsilon
par(mfrow = c(3, 3))
plot_geno_non_and_trans_hist(aa_trans_delta, aa_non_trans_delta, "AA")
plot_geno_non_and_trans_hist(a_out$trans_pheno_delta_edge, a_out$non_trans_pheno_delta_edge, "A")
plot_geno_non_and_trans_hist(b_out$trans_pheno_delta_edge, b_out$non_trans_pheno_delta_edge, "B")
plot_geno_non_and_trans_hist(c_out$trans_pheno_delta_edge, c_out$non_trans_pheno_delta_edge, "C")
plot_geno_non_and_trans_hist(d_out$trans_pheno_delta_edge, d_out$non_trans_pheno_delta_edge, "D")
plot_geno_non_and_trans_hist(e_trans_delta, e_non_trans_delta, "E")
plot_geno_non_and_trans_hist(f_trans_delta, f_non_trans_delta, "F")

# Make a table of epsilon score
num_data_sets <- 7
epsilon_mat <- matrix(0, nrow = num_data_sets, ncol = 6)
colnames(epsilon_mat) <- c("b_scaled_to_one", 
                           "scaled_to_one_TWICE",
                           "b_scaled_to_num_trans_edge", 
                           "nontrans_med", 
                           "bin", 
                           "bin_med")
row.names(epsilon_mat) <- c("aa", "a", "b", "c", "d", "e", "f")
epsilon_mat[1, ] <- c(ep_one_aa, ep_twice_aa, ep_num_aa, ep_nt_med_aa, ep_bin_aa, ep_bin_med_aa)
epsilon_mat[2, ] <- c(ep_one_a, ep_twice_a, ep_num_a, ep_nt_med_a, ep_bin_a, ep_bin_med_a)
epsilon_mat[3, ] <- c(ep_one_b, ep_twice_b, ep_num_b, ep_nt_med_b, ep_bin_b, ep_bin_med_b)
epsilon_mat[4, ] <- c(ep_one_c, ep_twice_c, ep_num_c, ep_nt_med_c, ep_bin_c, ep_bin_med_c)
epsilon_mat[5, ] <- c(ep_one_d, ep_twice_d, ep_num_d, ep_nt_med_d, ep_bin_d, ep_bin_med_d)
epsilon_mat[6, ] <- c(ep_one_e, ep_twice_e, ep_num_e, ep_nt_med_e, ep_bin_e, ep_bin_med_e)
epsilon_mat[7, ] <- c(ep_one_f, ep_twice_f, ep_num_f, ep_nt_med_f, ep_bin_f, ep_bin_med_f)

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

# Conclusion: epsilon with phenotype delta scaled to one seems like it has the best properties. 
# epsilon = sum(genotype_trasition * scaled_delta_pheno) / (sum(genotype_transition) + sum(scaled_delta_pheno) - numerator)
# max(sum(genotype_transition)) == Nedge(tree) (every edge is a transition -- which, admittedly is nearly impossible)
# min(sum(geotype_transition)) == 0 (no edge is a transition)

# max(sum(scaled_delta_pheno)) == Nedge(tree) (all delta pheno the same on each tree edge, so all get set to 1)
# min(sum(scaled_delta_pheno)) == 0 (no changes on any edge) -- but that would never happen b/c we're not using uniform phenotypes

# max(numerator) == Nedge (elementwise multiplication of 1 x 1 x Nedge(tree))
# min(numerator) == 0 (elementwise multiplication of 0 x 0 x Nedge(tree)))

# max(denominator) == 2 x Nedge (Nedge + Nedge - 0)
# min(denominator) == 0? technically 0 + 0 - 0, but this would never happen. More realistically :  0 + sum(scaled_delta_pheno)) - 0