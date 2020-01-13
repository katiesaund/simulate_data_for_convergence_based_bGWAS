generate_disc_mat <- function(tree_list){
  # Generate a huge discrete, binary matrix. 
  num_genotypes <- 1000 # TODO User defined input? Just some very large number?
  num_sim_trees <- length(tree_list)
  geno_mat_list <- geno_AR_mat_list <- rep(list(NULL), num_sim_trees)
  for (i in 1:num_sim_trees) {
    set.seed(1)
    
    # Create a matrix: 
    # Rows = tips then nodes
    # Columns = individual, simulated genotypes
    geno_tip_and_AR_mat <-
      replicate(num_genotypes,
                ape::rTraitDisc(tree_list[[i]],
                                ancestor = TRUE,
                                root.value = sample(x = c(1,2), size = 1, replace = FALSE, prob = c(0.5, 0.5)), # Root.value refers to which factor (1 or 2) of the two states (k=2, states are either 0 or 1).
                                type = "ER",
                                k = 2,
                                states = c(0,1)))
    geno_mat <- geno_tip_and_AR_mat[1:ape::Ntip(tree_list[[i]]), , drop = FALSE]
    storage.mode(geno_mat) <- "numeric"
    storage.mode(geno_tip_and_AR_mat) <- "numeric"
    cols_to_keep <- colSums(geno_mat) > 1 & colSums(geno_mat) < nrow(geno_mat) - 1
    geno_mat <- geno_mat[, cols_to_keep, drop = FALSE]
    geno_tip_and_AR_mat <- geno_tip_and_AR_mat[, cols_to_keep, drop = FALSE]
    geno_mat <- as.data.frame(geno_mat)
    geno_tip_and_AR_mat <- as.data.frame(geno_tip_and_AR_mat)
    colnames(geno_tip_and_AR_mat) <- colnames(geno_mat) <- paste0("sim", 1:ncol(geno_mat))
    geno_mat_list[[i]] <- geno_mat
    geno_AR_mat_list[[i]] <- geno_tip_and_AR_mat
  }
  return(geno_AR_mat_list)
}

# Phylogenetic signal functions -------
calculate_phylo_signal <- function(tree_list, binary_AR_mat_list){
  num_tree <- length(tree_list)
  d_stat_list <- rep(list(), num_tree)
  for (i in 1:num_tree) {
    bin_no_AR_mat <- binary_AR_mat_list[[i]][1:ape::Ntip(tree_list[[i]]), , drop = FALSE]
    d_stat_list[[i]] <- estimate_d_stat(bin_no_AR_mat, tree_list[[i]])
  }
  return(d_stat_list)
}

estimate_d_stat <- function(trait_df, tree){
  # Calculate the D statistic for the discrete trait we made
  tree$node.label <- NULL
  num_genotype <- ncol(trait_df)
  d_stat_vec <- rep(NA, num_genotype)
  for (i in 1:num_genotype) {
    temp_trait <- trait_df[, i, drop = FALSE]
    temp_trait <- convert_trait_vec_to_df(temp_trait, tree)
    compar_data_obj <- 
      caper::comparative.data(data = temp_trait,
                              phy = tree,
                              names.col = "ID")
    d_stat_results <- phylo.d2(data = compar_data_obj,
                               phy = tree,
                               binvar = "trait")
    d_stat_vec[i] <- d_stat_results$DEstimate
  }
  return(d_stat_vec)
}

convert_trait_vec_to_df <- function(trait_vec, tree){
  trait_df <- as.data.frame(trait_vec)
  trait_df <- cbind(tree$tip.label, trait_df)
  colnames(trait_df) <- c("ID", "trait")
  row.names(trait_df) <- NULL
  return(trait_df)
}

#' phylo.d2
#' @details This is a modification of the phylo.d function originally a part of
#' the caper R pacakge. The authors a different package, sensiPhy, rewrote the
#' phylo.d function to work with string variable name inputs rather than the
#' symbol of the variable name. I copied the function from:
#' https://github.com/paternogbc/sensiPhy/issues/186, code written by Caterina
#' Penone.
phylo.d2 <- function(data,
                     phy,
                     names.col,
                     binvar,
                     permut = 1000,
                     rnd.bias = NULL)
{
  if (!missing(data)) {
    if (!inherits(data, "comparative.data")) {
      if (missing(names.col)) {
        stop("names column is missing")
      }
      names.col <- deparse(substitute(names.col))
      data <- caicStyleArgs(data = data, phy = phy, names.col = names.col)
    }
  }
  binvar <- binvar
  bininds <- match(binvar, names(data$data))
  if (is.na(bininds))
    (stop("'", binvar, "' is not a variable in data."))
  ds <- data$data[, bininds]
  if (any(is.na(ds)))
    stop("'", binvar, "' contains missing values.")
  if (is.character(ds))
    ds <- as.factor(ds)
  if (length(unique(ds)) > 2)
    stop("'", binvar, "' contains more than two states.")
  if (length(unique(ds)) < 2)
    stop("'", binvar, "' only contains a single state.")
  propStates <- unclass(table(ds))
  propState1 <- propStates[1]/sum(propStates)
  names(dimnames(propStates)) <- binvar
  if (is.factor(ds))
    ds <- as.numeric(ds)
  if (!is.numeric(permut))
    (stop("'", permut, "' is not numeric."))
  if (!is.null(rnd.bias)) {
    rnd.bias <- deparse(substitute(rnd.bias))
    rnd.ind <- match(rnd.bias, names(data$data))
    if (is.na(rnd.ind))
      (stop("'", rnd.bias, "' is not a variable in data."))
    rnd.bias <- data$data[, rnd.bias]
  }
  el <- data$phy$edge.length
  elTip <- data$phy$edge[, 2] <= length(data$phy$tip.label)
  if (any(el[elTip] == 0))
    stop("Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate")
  if (any(el[!elTip] == 0))
    stop("Phylogeny contains zero length internal branches. Use di2multi.")
  ds.ran <- replicate(permut, sample(ds, prob = rnd.bias))
  if (is.null(data$vcv)) {
    vcv <- VCV.array(data$phy)
  }
  else {
    vcv <- data$vcv
  }
  ds.phy <- rmvnorm(permut, sigma = unclass(vcv))
  ds.phy <- as.data.frame(t(ds.phy))
  ds.phy.thresh <- apply(ds.phy, 2, quantile, propState1)
  ds.phy <- sweep(ds.phy, 2, ds.phy.thresh, "<")
  ds.phy <- as.numeric(ds.phy)
  dim(ds.phy) <- dim(ds.ran)
  ds.ran <- cbind(Obs = ds, ds.ran)
  ds.phy <- cbind(Obs = ds, ds.phy)
  dimnames(ds.ran) <- dimnames(ds.phy) <- list(data$phy$tip.label,
                                               c("Obs", paste("V", 1:permut, sep = "")))
  phy <- reorder(data$phy, "pruningwise")
  ds.ran.cc <- contrCalc(vals = ds.ran, phy = phy, ref.var = "V1",
                         picMethod = "phylo.d", crunch.brlen = 0)
  ds.phy.cc <- contrCalc(vals = ds.phy, phy = phy, ref.var = "V1",
                         picMethod = "phylo.d", crunch.brlen = 0)
  ransocc <- colSums(ds.ran.cc$contrMat)
  physocc <- colSums(ds.phy.cc$contrMat)
  if (round(ransocc[1], digits = 6) != round(physocc[1], digits = 6))
    stop("Problem with character change calculation in phylo.d")
  obssocc <- ransocc[1]
  ransocc <- ransocc[-1]
  physocc <- physocc[-1]
  soccratio <- (obssocc - mean(physocc))/(mean(ransocc) -
                                            mean(physocc))
  soccpval1 <- sum(ransocc < obssocc)/permut
  soccpval0 <- sum(physocc > obssocc)/permut
  dvals <- list(DEstimate = soccratio,
                Pval1 = soccpval1,
                Pval0 = soccpval0,
                Parameters = list(Observed = obssocc,
                                  MeanRandom = mean(ransocc),
                                  MeanBrownian = mean(physocc)),
                StatesTable = propStates,
                Permutations = list(random = ransocc, brownian = physocc),
                NodalVals = list(observed = ds.ran.cc$nodVal[, 1, drop = FALSE],
                                 random = ds.ran.cc$nodVal[, -1, drop = FALSE],
                                 brownian = ds.phy.cc$nodVal[, -1, drop = FALSE]),
                binvar = binvar,
                data = data,
                nPermut = permut,
                rnd.bias = rnd.bias)
  class(dvals) <- "phylo.d"
  return(dvals)
}
# ----

select_BM_traits <- function(binary_mat_list, 
                             phylo_signal_list, 
                             num_to_pick = NULL) {
  num_trees <- length(binary_mat_list)
  BM_trait_names_list <- rep(list(), num_trees)
  for (i in 1:num_trees) {
    BM_trait_logical <- phylo_signal_list[[i]] < 0.05 & phylo_signal_list[[i]] > -0.05
    BM_colnames <- colnames(binary_mat_list[[i]][, BM_trait_logical, drop = FALSE])
    if (!is.null(num_to_pick)) {
      if (num_to_pick <= length(BM_colnames)){
       BM_colnames <- BM_colnames[1:num_to_pick]
      } else {
      stop("Not enough BM traits")
      }
    }
    BM_trait_names_list[[i]] <- BM_colnames
  } 
  return(BM_trait_names_list)
}

select_WN_traits <- function(binary_mat_list, 
                             phylo_signal_list, 
                             num_to_pick = NULL) {
  num_trees <- length(binary_mat_list)
  WN_trait_names_list <- rep(list(), num_trees)
  for (i in 1:num_trees) {
    WN_trait_logical <- phylo_signal_list[[i]] < 1.05 & phylo_signal_list[[i]] > 0.95
    WN_colnames <- colnames(binary_mat_list[[i]][, WN_trait_logical, drop = FALSE])
    if (!is.null(num_to_pick)) {
      if (num_to_pick <= length(WN_colnames)){
        WN_colnames <- WN_colnames[1:num_to_pick]
      } else {
        stop("Not enough WN traits")
      }
    }
    WN_trait_names_list[[i]] <- WN_colnames
  } 
  return(WN_trait_names_list)
}

subsample_to_phenotypes <- function(binary_AR_mat_list, 
                                    phenotype_names_list){
  num_mat <- length(binary_AR_mat_list)
  temp_mat_list <- rep(list(), num_mat)
  for (i in 1:num_mat){
    temp_mat_list[[i]] <- binary_AR_mat_list[[i]][, colnames(binary_AR_mat_list[[i]]) %in% phenotype_names_list[[i]], drop = FALSE]
  }
  return(temp_mat_list)
}

subsample_to_genotypes <- function(binary_AR_mat_list, 
                                   phylo_signal_list,
                                   lower_bound = NULL, 
                                   upper_bound = NULL, 
                                   num_genos = NULL){
  num_mat <- length(binary_AR_mat_list)
  geno_mat_list <- binary_AR_mat_list
  for (i in 1:num_mat) {
    if (!is.null(lower_bound)) {
      greater_than_lower_log <- phylo_signal_list[[i]] > lower_bound
      phylo_signal_list[[i]] <- phylo_signal_list[[i]][greater_than_lower_log]
      geno_mat_list[[i]] <- geno_mat_list[[i]][, greater_than_lower_log, drop = FALSE]
    }
    if (!is.null(upper_bound)) {
      less_than_upper_log <- phylo_signal_list[[i]] < upper_bound
      phylo_signal_list <- phylo_signal_list[[i]][less_than_upper_log]
      geno_mat_list[[i]] <- geno_mat_list[[i]][, less_than_upper_log, drop = FALSE]
    }
    if (!is.null(num_genos)) {
      if (ncol(geno_mat_list[[i]]) < num_genos) {
        stop(paste0("Geno mat", i, "doesn't have enough genotypes"))
      }
    }
  }
  return(geno_mat_list)
}

combine_phenotype_names_lists <- function(list_1, list_2){
  num_names <- length(list_1)
  name_list <- rep(list(), num_names)
  for (i in 1:num_names) {
    name_list[[i]] <- unique(c(list_1[[i]], list_2[[i]]))
  }
  return(name_list)
}

## OLD STUFF BELOW

save_tree_specific_mat <- function(sim_geno_mat,
                                   tree,
                                   index){
  if (nrow(sim_geno_mat) != Ntip(tree)) {
    stop("genotype matrix should have 1 row per tree tip")
  }

  row.names(sim_geno_mat) <- tree$tip.label

  write.table(sim_geno_mat,
              sep = "\t",
              row.names = TRUE,
              file = paste0("../data/",
                            "simulated_genotype_",
                            index,
                            ".tsv"))
  write.tree(tree,
             file = paste0("../data/",
                           "simulated_tree_",
                           index,
                           ".tree"))
} # end save_tree_specific_mat()



# make_discrete_phenotypes <- function(tree_list, num_pheno){
#   # Generate two discrete traits that follows a given tree, one BM, one WN
#   # A D statistic == 0 -> BM (high phylogenetic signal)
#   # A D statistic == 1 -> Random (low phylogenetic signal)
# 
#   num_sim_trees <- length(tree_list)
#   set.seed(1)
#   for (j in 1:num_pheno) {
#     for (i in 1:num_sim_trees) {
#       # Generate a trait with high phylogenetic signal (D stat ~0) ---------------
#       discrete_BM_trait <- generate_BM_trait(tree_list[[i]])
# 
#       # Generate a trait with low phylogenetic signal (D stat => 1) --------------
#       discrete_WN_trait <- generate_WN_trait(tree_list[[i]])
# 
#       # Save the discrete traits -----------------------------------------------
#       row.names(discrete_BM_trait) <- discrete_BM_trait[, 1]
#       discrete_BM_trait <- discrete_BM_trait[, 2, drop = FALSE]
#       row.names(discrete_WN_trait) <- discrete_WN_trait[, 1]
#       discrete_WN_trait <- discrete_WN_trait[, 2, drop = FALSE]
#       colnames(discrete_WN_trait) <- colnames(discrete_WN_trait) <- "pheno"
#       write.table(discrete_BM_trait,
#                   sep = "\t",
#                   file = paste0("../data/",
#                                 "discrete_pheno_BM_tree_",
#                                 i,
#                                 "_pheno_",
#                                 j,
#                                 ".tsv"))
#       write.table(discrete_WN_trait,
#                   sep = "\t",
#                   file = paste0("../data/",
#                                 "discrete_pheno_WN_tree_",
#                                 i,
#                                 "_pheno_",
#                                 j,
#                                 ".tsv"))
#     }
#   }
# }


# make_genotype_mat <- function(tree_list, num_genotypes){
#   # Generate a matrix of discrete genotypes that each follows BM.
#   # A D statistic == 0 -> BM (high phylogenetic signal)
#   # A D statistic == 1 -> Random (low phylogenetic signal)
#   num_sim_trees <- length(tree_list)
#   set.seed(1)
#   for (i in 1:num_sim_trees) {
#     set.seed(1)
#     geno_mat <-
#       replicate(num_genotypes,
#                 ape::rTraitDisc(tree_list[[i]],
#                                 root.value = sample(x = c(1,2), size = 1, replace = FALSE, prob = c(0.5, 0.5)), # Root.value refers to which factor (1 or 2) of the two states (k=2, states are either 0 or 1).
#                                 type = "ER",
#                                 k = 2,
#                                 states = c(0,1)))
#     storage.mode(geno_mat) <- "numeric"
#     geno_mat <- geno_mat[, colSums(geno_mat) > 1]
#     geno_mat <- geno_mat[, colSums(geno_mat) < nrow(geno_mat) - 1]
#     geno_mat <- subset_to_N_state_flips(geno_mat, 2)
#     geno_mat <- as.data.frame(geno_mat)
#     colnames(geno_mat) <- paste0("sim", 1:ncol(geno_mat))
#     save_tree_specific_mat(geno_mat, tree_list[[i]], i)
#   }
#   return(NULL)
# }



#' #' Title
#' #'
#' #' @param tree
#' #'
#' #' @return discrete_WN_trait Data.frame. 2 columns. Colnames: ID & trait.
#' #'   1st column = tree tip labels. 2nd column = binary (0/1). No row names.
#' generate_WN_trait <- function(tree){
#'   d_stat_lack_low_signal <- TRUE
#'   while (d_stat_lack_low_signal) { # While the D-stat is "BM", do the following:
#'     continuous_BM_trait <- phytools::fastBM(tree = tree,
#'                                             bounds = c(0, 1),
#'                                             sig2 = 100)
#'     discrete_WN_trait <- round(continuous_BM_trait, 0)
#'     discrete_WN_trait <- convert_trait_vec_to_df(discrete_WN_trait,
#'                                                  tree)
#'     current_d_stat <- estimate_d_stat(discrete_WN_trait, tree)
#'     d_stat_lack_low_signal <- current_d_stat < 0.95 | current_d_stat > 1.05
#'   }
#'   return(discrete_WN_trait)
#' }

#' #' Title
#' #'
#' #' @param tree
#' #'
#' #' @return discrete_BM_trait Data.frame. 2 columns. Colnames: ID & trait.
#' #'   1st column = tree tip labels. 2nd column = binary (0/1). No row names.
#' generate_BM_trait <- function(tree){
#'   d_stat_lack_high_signal <- TRUE
#'   while (d_stat_lack_high_signal) { # While the D-stat is "random", do the following:
#' 
#'     # Generate a temporary continuous trait to use
#'     indep_var <- phylolm::rTrait(n = 1, phy = tree)
#'     indep_var <- cbind(rep(1, length(indep_var)), indep_var)
#' 
#'     # Simulate a binary trait (0/1) along the current tree
#'     num_traits <- 1
#'     while (num_traits == 1) {
#'       discrete_BM_trait <-
#'         # Since n = 1 in phylolm::rbinTrait, the output is a numeric vector of
#'         #    0-1 values with names from the tip labels in the tree.
#'         phylolm::rbinTrait(n = 1, # num replicates
#'                            phy = tree, # tree
#'                            alpha = 1, # phylogenetic correlation
#'                            beta = c(-1, 0.5), # coeffients for the logistic
#'                            #    regression model
#'                            X = indep_var) # X = Design matrix, one row per tree
#'       #    tip
#' 
#'       discrete_BM_trait <- convert_trait_vec_to_df(discrete_BM_trait, tree)
#'       num_traits <- length(unique(discrete_BM_trait$trait))
#'       # Can't move on to next section until the discrete trait has two traits
#'       #     (num_traits == 2) because phylo.d2 requires two traits
#'     }
#'     current_d_stat <- estimate_d_stat(discrete_BM_trait, tree)
#'     d_stat_lack_high_signal <- current_d_stat > 0.05 | current_d_stat < -0.05
#' 
#'     
#'     # So if the D statistic is random it will restart the while loop,
#'     #    if the D statistic fits BM the code will move on
#'   }
#'   return(discrete_BM_trait)
#' }