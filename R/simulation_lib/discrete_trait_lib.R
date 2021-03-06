#' Generate a binary data.frame for each  tree
#'
#' @description The function generates a unique data.frame for each tree. The 
#'   function attempts to create a number of traits equal to num_genotypes, but
#'   the application of several filters may reduce this number. The filters 
#'   remove traits that are too ubiquitous or too rare and duplicated traits.
#'
#' @param tree_list List of phylogenetic trees. 
#' @param num_genotypes Integer. The number of discrete traits to initially 
#'   produce. 
#'
#' @return geno_AR_df_list A list with length equal to the number of trees in
#'   the input tree_list. The list entries are dataframes. Each dataframe has 
#'   nrow = number of discrete traits that pass the quality control steps. Ncol 
#'   = number of tips in the corresponding tree. 
#' @export
#'
generate_binary_df_list <- function(tree_list, num_genotypes){
  # Generate a huge discrete, binary matrix.
  num_sim_trees <- length(tree_list)
  geno_AR_df_list <- rep(list(NULL), num_sim_trees)
  set.seed(1)
  for (i in 1:num_sim_trees) {
    # Create a matrix:
    # Rows = tips then nodes
    # Columns = individual, simulated genotypes
    geno_tip_and_AR_mat <-
      replicate(num_genotypes,
                ape::rTraitDisc(tree_list[[i]],
                                ancestor = TRUE,
                                root.value = sample(x = c(1,2),
                                                    size = 1, 
                                                    replace = FALSE,
                                                    prob = c(0.5, 0.5)), 
                                # Root.value refers to which factor (1 or 2) of 
                                # the two states (k=2, states are either 0 or 1).
                                type = "ER",
                                k = 2,
                                states = c(0,1)))
    geno_mat <- geno_tip_and_AR_mat[1:ape::Ntip(tree_list[[i]]), , drop = FALSE]
    storage.mode(geno_mat) <- "numeric"
    storage.mode(geno_tip_and_AR_mat) <- "numeric"
    cols_to_keep <- 
      colSums(geno_mat) > 1 & colSums(geno_mat) < nrow(geno_mat) - 1
    geno_tip_and_AR_mat <- geno_tip_and_AR_mat[, cols_to_keep, drop = FALSE]
    geno_tip_and_AR_mat <- unique(geno_tip_and_AR_mat, MARGIN = 2)
    geno_tip_and_AR_mat <- as.data.frame(geno_tip_and_AR_mat)
    colnames(geno_tip_and_AR_mat) <- paste0("sim", 1:ncol(geno_tip_and_AR_mat))
    geno_AR_df_list[[i]] <- geno_tip_and_AR_mat
  }
  return(geno_AR_df_list)
}

# Phylogenetic signal functions -------

#' Wrapper function for estimate_d_stat() to work on lists
calculate_phylo_signal <- function(tree_list, binary_AR_mat_list){
  num_tree <- length(tree_list)
  d_stat_list <- rep(list(), num_tree)
  for (i in 1:num_tree) {
    bin_no_AR_mat <- 
      binary_AR_mat_list[[i]][1:ape::Ntip(tree_list[[i]]), , drop = FALSE]
    d_stat_list[[i]] <- estimate_d_stat(bin_no_AR_mat, tree_list[[i]])
  }
  return(d_stat_list)
}

#' Calculate the D statistic for the discrete trait we made
estimate_d_stat <- function(trait_df, tree){
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

#' Select the traits that follow Brownian motion
select_BM_traits <- function(binary_mat_list,
                             phylo_signal_list,
                             num_to_pick = NULL) {
  num_trees <- length(binary_mat_list)
  BM_trait_names_list <- rep(list(), num_trees)
  for (i in 1:num_trees) {
    BM_trait_logical <- phylo_signal_list[[i]] < 0.05 & phylo_signal_list[[i]] > -0.05
    BM_colnames <- colnames(binary_mat_list[[i]][, BM_trait_logical, drop = FALSE])
    if (!is.null(num_to_pick)) {
      if (num_to_pick <= length(BM_colnames)) {
       BM_colnames <- BM_colnames[1:num_to_pick]
      } else {
      stop("Not enough BM traits")
      }
    }
    BM_trait_names_list[[i]] <- BM_colnames
  }
  return(BM_trait_names_list)
}

#' Select the traits that follow white noise
select_WN_traits <- function(binary_mat_list,
                             phylo_signal_list,
                             num_to_pick = NULL) {
  num_trees <- length(binary_mat_list)
  WN_trait_names_list <- rep(list(), num_trees)
  for (i in 1:num_trees) {
    WN_trait_logical <- 
      phylo_signal_list[[i]] < 1.05 & phylo_signal_list[[i]] > 0.95
    WN_colnames <- 
      colnames(binary_mat_list[[i]][, WN_trait_logical, drop = FALSE])
    if (!is.null(num_to_pick)) {
      if (num_to_pick <= length(WN_colnames)) {
        WN_colnames <- WN_colnames[1:num_to_pick]
      } else {
        stop("Not enough WN traits")
      }
    }
    WN_trait_names_list[[i]] <- WN_colnames
  }
  return(WN_trait_names_list)
}

#' Subset to just those phenotypes in the provided list: phenotype_names_list
subsample_to_phenotypes <- function(binary_AR_mat_list,
                                    binary_conf_mat_list,
                                    phenotype_names_list){
  num_mat <- length(binary_AR_mat_list)
  AR_mat_list <- conf_mat_list <- rep(list(), num_mat)
  for (i in 1:num_mat) {
    AR_mat_list[[i]] <- 
      binary_AR_mat_list[[i]][, colnames(binary_AR_mat_list[[i]]) %in% phenotype_names_list[[i]], 
                              drop = FALSE]
    conf_mat_list[[i]] <- 
      binary_conf_mat_list[[i]][, colnames(binary_conf_mat_list[[i]]) %in% phenotype_names_list[[i]], 
                                drop = FALSE]

  }
  return(list("AR_mat" = AR_mat_list,
              "conf_mat" = conf_mat_list))
}

#' Only keep genotypes whose phylogenetic signal is within the given range
select_geno_within_range <- function(binary_AR_mat_list,
                                     binary_AR_conf_list,
                                     phylo_signal_list,
                                     lower_bound = NULL,
                                     upper_bound = NULL,
                                     min_genos = 1){
  num_mat <- length(binary_AR_mat_list)
  geno_mat_list <- binary_AR_mat_list
  for (i in 1:num_mat) {
    if (!is.null(lower_bound)) {
      greater_than_lower_log <- phylo_signal_list[[i]] > lower_bound
      phylo_signal_list[[i]] <- phylo_signal_list[[i]][greater_than_lower_log]
      geno_mat_list[[i]] <- 
        geno_mat_list[[i]][, greater_than_lower_log, drop = FALSE]
      binary_AR_conf_list[[i]] <- 
        binary_AR_conf_list[[i]][, greater_than_lower_log, drop = FALSE]
    }
    if (!is.null(upper_bound)) {
      less_than_upper_log <- phylo_signal_list[[i]] < upper_bound
      phylo_signal_list[[i]] <- phylo_signal_list[[i]][less_than_upper_log]
      geno_mat_list[[i]] <-
        geno_mat_list[[i]][, less_than_upper_log, drop = FALSE]
      binary_AR_conf_list[[i]] <-
        binary_AR_conf_list[[i]][, less_than_upper_log, drop = FALSE]

    }
    if (!is.null(min_genos)) {
      if (ncol(geno_mat_list[[i]]) < min_genos) {
        stop(paste0("Geno mat", i, "doesn't have enough genotypes"))
      }
    }
  }

  return(list("conf_mat" = binary_AR_conf_list,
              "AR_mat" = geno_mat_list))
}

#' Add phenotypes that are modeled by white noise to the origin set that were
#'   modeled with Brownian motion 
add_WN <- function(binary_AR_mat_list, tree_list) {
  num_trees <- length(tree_list)
  for (i in 1:num_trees) {
    num_trait <- ncol(binary_AR_mat_list[[i]])
    num_tip <- ape::Ntip(tree_list[[i]])
    num_to_add <- round(num_trait / 4, 0)
    random_col_to_add <- flipped_to_add_10p <- flipped_to_add_25p <- 
      flipped_to_add_40p <- 
      binary_AR_mat_list[[i]][, 1:num_to_add, drop = FALSE]
    for (j in 1:ncol(random_col_to_add)) {
      
      # select some columns to work with
      current_col <- random_col_to_add[1:num_tip, j]
      
      # Completely jumble the selected columns
      jumbled_tips <- sample(current_col,
                             size = length(current_col),
                             replace = FALSE)
      random_col_to_add[1:num_tip, j] <- jumbled_tips
      
      # Flip 10% of the bits
      flipped_tip_10p <- flip_some_tips(current_col, flip_freq = 0.10)
      flipped_to_add_10p[1:num_tip, j] <- flipped_tip_10p
      
      # Flip 25% of the bits
      flipped_tip_25p <- flip_some_tips(current_col, flip_freq = 0.25)
      flipped_to_add_25p[1:num_tip, j] <- flipped_tip_25p
      
      # Flip 40% of the bits
      flipped_tip_40p <- flip_some_tips(current_col, flip_freq = 0.40)
      flipped_to_add_40p[1:num_tip, j] <- flipped_tip_40p
    }
    binary_AR_mat_list[[i]] <- cbind(binary_AR_mat_list[[i]],
                                     random_col_to_add,
                                     flipped_to_add_10p, 
                                     flipped_to_add_25p, 
                                     flipped_to_add_40p)

    # Remove traits found in 0, 1, N, or N-1 samples
    temp_only_tips <- binary_AR_mat_list[[i]][1:num_tip, , drop = FALSE]
    cols_to_keep <-
      colSums(temp_only_tips) > 1 & colSums(temp_only_tips) < (nrow(temp_only_tips) - 1)
    binary_AR_mat_list[[i]] <- binary_AR_mat_list[[i]][, cols_to_keep, drop = FALSE]
    
    # Remove duplicate columns (genotypes)
    binary_AR_mat_list[[i]] <- as.data.frame(t(unique(t(binary_AR_mat_list[[i]])))) 
    
    # Write genotype names (colnames)
    colnames(binary_AR_mat_list[[i]]) <-
      paste0("sim", 1:ncol(binary_AR_mat_list[[i]]))
  }
  return(binary_AR_mat_list)
}

#' Flip the bits of some tips given a user provided frequency
#' flip_freq should be between 0 and 1
flip_some_tips <- function(numeric_vec, flip_freq = 0.10) {
  len <- length(numeric_vec)
  index <- 1:len
  indices_to_flip <- sample(index, round(len * flip_freq, 0), replace = FALSE)
  for (i in indices_to_flip) {
    numeric_vec[i] <- as.numeric(!numeric_vec[i])
  }
  return(numeric_vec)
}
