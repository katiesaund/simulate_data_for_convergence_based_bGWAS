#' Wrapper function for ancestral_reconstruction_by_ML() to run on lists
ancestral_reconstruction <- function(AR_mat_list, tree_list, disc_or_cont){
  num_trees <- length(tree_list)
  geno_recon_and_conf_list <- rep(list(0), num_trees)
  AR_mat <- conf_mat <- rep(list(matrix(0)), num_trees)
  recon_edge_mat <- rep(list(), num_trees)

  for (i in 1:num_trees) {
    temp_tree <- tree_list[[i]]
    num_tips <- ape::Ntip(temp_tree)
    tip_mat <- AR_mat_list[[i]][1:num_tips, , drop = FALSE]
    num_geno <- ncol(tip_mat)
    geno_recon_and_conf <- rep(list(0), num_geno)

    num_row <- ape::Ntip(temp_tree) + ape::Nnode(temp_tree)
    AR_mat[[i]] <- conf_mat[[i]] <- matrix(0, nrow = num_row, ncol = num_geno)
    recon_edge_mat[[i]] <- rep(list(), num_geno)

    for (j in 1:num_geno) {
      geno_recon_and_conf[[j]] <- ancestral_reconstruction_by_ML(temp_tree, tip_mat, j, disc_or_cont)
      AR_mat[[i]][, j] <- geno_recon_and_conf[[j]]$tip_and_node_recon
      conf_mat[[i]][, j] <- geno_recon_and_conf[[j]]$tip_and_node_rec_conf
      recon_edge_mat[[i]][[j]] <- geno_recon_and_conf[[j]]$recon_edge_mat
      storage.mode(recon_edge_mat[[i]][[j]]) <- "numeric"
    }
    colnames(AR_mat[[i]]) <- colnames(conf_mat[[i]]) <- colnames(AR_mat_list[[i]])
    
    if (nrow(AR_mat_list[[i]]) == nrow(conf_mat[[i]])) {
      row.names(AR_mat[[i]]) <- row.names(conf_mat[[i]]) <- row.names(AR_mat_list[[i]])
    } else {
      row.names(AR_mat[[i]]) <- row.names(conf_mat[[i]]) <- c(row.names(AR_mat_list[[i]]), paste0("node_", 1:ape::Nnode(tree_list[[i]])))
    }
    
    geno_recon_and_conf_list[[i]] <- geno_recon_and_conf
  }

  return(list("conf_mat" = conf_mat,
              "AR_mat" = AR_mat, 
              "recon_edge_mat" = recon_edge_mat))
}

#' ancestral_reconstruction_by_ML
#'
#' @description Compute ancestral reconstruction from a continuous phenotype or
#' discrete genotype or phenotype.
#'
#' @details Note on the marginal setting in the function ape::ace(). Marginal =
#'   FALSE means that the marginal is in fact calculated. Do not set marginal to
#'   TRUE, because this only does the condition (only downstream edges
#'   considered). Ace never calculates the joint distribution.
#'
#' @param tr Phylo. Rooted phylogenetic tree.
#' @param mat Matrix. Either the phenotype matrix or the genotype matrix. Dim:
#'   nrow = Ntip(tr) x ncol = {1 if phenotype or number of genotypes}.
#' @param num Numeric. Indicating the row of the matrix from which the
#'   ancestral reconstruction is to be built.
#' @param disc_cont Character. Either "discrete" or "continuous."
#'
#' @return results a list with multiple outputs:
#'   * node_anc_rec Vector. Ancestral reconstruction for each node. Length =
#'       Nnode(tr).
#'   * tip_and_node_rec_conf A vector with a 1 for each tip with the node
#'       ancestral reconstruction confidence appended. Length = Nnode(tr) +
#'       Ntip(tr).
#'   * recon_edge_mat Matrix. Ancestral reconstruction configured into the edge
#'        matrix. Dim: nrow = Nedge(tr) x ncol = 2. 1st column == parent node.
#'        2nd column == child node.
#'   * tip_and_node_recon A vector with the tip values followed by the node
#'       ancestral reconstruction. Length == Ntip(tr) + Nnode(tr).
#'
#' @noRd
#'
ancestral_reconstruction_by_ML <- function(tr, mat, num, disc_cont){
  # Function -------------------------------------------------------------------
  # Compute ancestral reconstruction
  ML_significance_threshold <- .875 # ML literature suggests that a ratio of 7:1
                                    # suggests a high confidence ancestral
                                    # reconstruction per node .875/.125 = 7.

  if (disc_cont == "continuous") {
    # This is only for a continuous phenotype
    recon_method <- "pic" # To resolve this bug:
      # Error in nlm(function(p) dev.BM(p), p = c(1, rep(mean(x), nb.node)), 
      # hessian = TRUE) : missing value in parameter 
    
    # RECONSTRUCTION
    cont_results <- continuous_ancestral_reconstruction(tr,
                                                        mat,
                                                        num,
                                                        disc_cont,
                                                        recon_method)
    ML_anc_rec <- cont_results$ML_anc_rec
    tip_and_node_recon <- cont_results$tip_and_node_recon

    # CONFIDENCE IN RECONSTRUCTION
    tip_and_node_anc_rec_conf <-
      continuous_get_recon_confidence(tip_and_node_recon)

  } else {
    # This is always the choice for genotypes and discrete phenotypes
    recon_method <- "ML" # ML == Maximum Likelihood.
    
    # RECONSTRUCTION
    discrete_results <-
      discrete_ancestral_reconstruction(tr, mat, num, disc_cont, recon_method)
    ML_anc_rec <- discrete_results$ML_anc_rec
    tip_and_node_recon <- discrete_results$tip_and_node_recon

    # CONFIDENCE IN RECONSTRUCTION
    tip_and_node_anc_rec_conf <-
      discrete_get_recon_confidence(discrete_results$reconstruction,
                                    tr,
                                    ML_significance_threshold)
  }
  reconstruction_as_edge_mat <- convert_to_edge_mat(tr, tip_and_node_recon)

  # Return outputs -------------------------------------------------------------
  results <- list("node_anc_rec" = ML_anc_rec,
                  "tip_and_node_rec_conf" = tip_and_node_anc_rec_conf,
                  "recon_edge_mat" = reconstruction_as_edge_mat,
                  "tip_and_node_recon" = tip_and_node_recon)
  return(results)
}

#' continuous_ancestral_reconstruction
#'
#' @description Perform ancestral state reconstruction using ape::ace() function
#'   on continuous phenotype.
#'
#' @param tr Phylo.
#' @param mat Matrix. Continuous phenotype. Dim: nrow = Ntip(tr) x ncol = 1.
#' @param num  Number. Column index for matrix.
#' @param disc_cont Character string. Must be "continuous."
#' @param recon_method Character string. Either "ML", "REML", "pic", or "GLS."
#'
#' @return List of two outputs:
#'  * ML_anc_rec. Reconstruction values for each node.
#'  * tip_and_node_recon. Vector. Phenotype values for each tip followed by the
#'      ancestral reconstruction at each node. Length = Ntip(tr) + Nnode(tr).
#'
#' @noRd
#'
continuous_ancestral_reconstruction <- function(tr,
                                                mat,
                                                num,
                                                disc_cont,
                                                recon_method){
  # Function -------------------------------------------------------------------
  set.seed(1)
  reconstruction <- ape::ace(mat[, num, drop = TRUE],
                        tr,
                        model = "BM",
                        type = disc_cont,
                        method = recon_method,
                        marginal = FALSE)

  # Vector containing reconstruction data for all internal nodes N+ where tips
  # are 1-N.
  ML_anc_rec <- reconstruction$ace
  tip_and_node_recon <- c(mat[, num, drop = TRUE], ML_anc_rec)
  names(tip_and_node_recon) <- c(1:sum(ape::Ntip(tr), ape::Nnode(tr)))

  # Return output --------------------------------------------------------------
  return(list("ML_anc_rec" = ML_anc_rec,
              "tip_and_node_recon" = tip_and_node_recon))
}

#' continuous_get_recon_confidence
#'
#' @description Given a vector of a continuous trait reconstruction, generate a
#'   dummy confidence vector. Ancestral reconstruction for continuous values
#'   only gives back a 95% CI. We can't use any of this information to decide
#'   which nodes are low confidence so treat all reconstructed values as high
#'   confidence, which is stored as a value 1.
#'
#' @param recon_vector Numeric vector. Values of the reconstruction. Length ==
#'   Ntip(tr) + Nnode(tr).
#'
#' @return tip_and_node_anc_rec_conf: Vector of all 1s, indicating 'high'
#'   confidence. Length == Ntip(tr) + Nnode(tr). Tips folowed by nodes.
#'
#' @noRd
#'
continuous_get_recon_confidence <- function(recon_vector){
  # Function -------------------------------------------------------------------
  tip_and_node_anc_rec_conf <- rep(1, length(recon_vector))

  # Return output --------------------------------------------------------------
  return(tip_and_node_anc_rec_conf)
}

#' convert_to_edge_mat
#'
#' @description Convert the reconstruction to be in the same format as
#'   tree$edge, where each row is an edge. This format will make it much easier
#'   to calculate the phenotype change on each edge for continuous phenotypes.
#'
#' @param tr Phylo.
#' @param tip_and_node_reconstruction Vector. Ordered by tips then by nodes.
#'   Length = Ntip(tr) = Nnode(tr). Observed values at each tip and ancestral
#'   reconstruction values at each node.
#'
#' @return reconstruction_as_edge_mat. Matrix. Dim = nrow = Ntip(tr) x ncol = 2.
#'   1st column = parent node. 2nd column = child node. Ancestral reconstruction
#'   values for each node.
#'
#' @noRd
convert_to_edge_mat <- function(tr, tip_and_node_reconstruction){
  # Function -------------------------------------------------------------------
  reconstruction_as_edge_mat <- tr$edge
  for (k in 1:nrow(reconstruction_as_edge_mat)) {
    reconstruction_as_edge_mat[k, 1] <-
      tip_and_node_reconstruction[tr$edge[k, 1]]
    reconstruction_as_edge_mat[k, 2] <-
      tip_and_node_reconstruction[tr$edge[k, 2]]
  }

  # Return output --------------------------------------------------------------
  return(reconstruction_as_edge_mat)
}

#' Perform ancestral reconstruction of discrete data
#'
#' @description Perform ancestral state reconstruction using ape::ace() function
#'   on discrete phenotype or genotype. The best model to describe the
#'   probabilities of state changes is decided by the subfunction:
#'   build_better_reconstruction(), which chooses either 'ARD', all rates
#'    different, or ER', equal rates and returns the better of the two
#'    reconstructions.
#'
#' @param tr Phylo.
#' @param mat Matrix. Genotype or discrete phenotype matrix.
#' @param num Number. Column index for matrix.
#' @param disc_cont Character. Should always be "discrete."
#' @param recon_method Character string. Either "ML", "REML", "pic", or "GLS."
#'
#' @return List with multiple outputs:
#'   * ML_anc_rec: Vector. Reconstruction (numeric) at node values only. Length
#'       = Nnode(tr).
#'   * tip_and_node_recon: A vector with the tip values followed by the node
#'       ancestral reconstruction. Length = Nnode(tr) + Ntip(tr).
#'   * reconstruction: Ape::ace() object. Contains multiple pieces of
#'       information. Class = "ace." Type = "list".
#'       $loglik. Number. Log likelihood of the ancestral reconstruction.
#'       $rates. Number.
#'       $se. Number. Standard error.
#'       $index.matrix. Matrix.
#'       $lik.anc. Matrix. Likelihood for each state at each tree node.
#'       $call. Record of the ace() call. Class = "call."
#'
#' @noRd
discrete_ancestral_reconstruction <- function(tr,
                                              mat,
                                              num,
                                              disc_cont,
                                              recon_method){
  # Check inputs ---------------------------------------------------------------
  if (disc_cont != "discrete") {
    stop("Discrete ancestral reconstruction only")
  }

  # Function -------------------------------------------------------------------
  set.seed(1)
  reconstruction <-
    build_better_reconstruction(mat, tr, disc_cont, num, recon_method)

  # Extract the mostly likely character state using which.max
  ML_anc_rec <-
    as.numeric(colnames(reconstruction$lik.anc)[apply(reconstruction$lik.anc,
                                                      1,
                                                      which.max)])
  names(ML_anc_rec) <- c( (ape::Ntip(tr) + 1):(ape::Ntip(tr) + ape::Nnode(tr)))
  tip_and_node_recon <- c(mat[, num, drop = TRUE], ML_anc_rec)
  names(tip_and_node_recon) <- c(1:sum(ape::Ntip(tr), ape::Nnode(tr)))

  # Return outputs -------------------------------------------------------------
  return(list("tip_and_node_recon" = tip_and_node_recon,
              "ML_anc_rec" = ML_anc_rec,
              "reconstruction" = reconstruction))
}

#' discrete_get_recon_confidence
#'
#' @description Given the ancestral reconstruction identify which nodes are high
#'   confidence based on the maximum likelihood of the reconstruction and the
#'   confidence treshold (0.875; a ratio of about 7:1; based on ML literature).
#'   Confidence at the tips is assigned 1 because these values were measured and
#'   therefore have total confidence in the value.
#'
#' @param recon Ape::ace() object. Contains multiple pieces of information.
#'   Class = "ace." Type = "list." Ancestral reconstruction.
#'   * $loglik. Number. Log likelihood of the ancestral reconstruction.
#'   * $rates. Number.
#'   * $se. Number. Standard error.
#'   * $index.matrix. Matrix.
#'   * $lik.anc. Matrix. Likelihood for each state at each tree node.
#'   * $call. Record of the ace() call. Class = "call."
#' @param tr Phylo.
#' @param ML_cutoff Number between 0 and 1. If the maximum likelihood of the
#'   reconstruction is above this number it's considered a high confidence
#'   reconstruction, otherwise it's a low confidence reconstruction.
#'
#' @noRd
#'
discrete_get_recon_confidence <- function(recon, tr, ML_cutoff){
  # Function -------------------------------------------------------------------
  # Get the highest confidence value at each node
  anc_rec_confidence <- apply(recon$lik.anc, 1, max)

  # Count all tips as high confidence
  tip_and_node_anc_rec_conf <- c(rep(1, ape::Ntip(tr)), anc_rec_confidence)

  # Count anything lower than threshold as low confidence
  tip_and_node_anc_rec_conf <-
    discretize_conf_with_cutoff(tip_and_node_anc_rec_conf, ML_cutoff)

  # Return output --------------------------------------------------------------
  return(tip_and_node_anc_rec_conf)
}

#' Reconstruct discrete data and return reconstruction built with the better
#' model
#'
#' @description Given a discrete genotype or phenotype, build an ancestral
#'   reconstruction using an equal rates model ("ER") and then repeat with an
#'   all rates different model ("ARD"). Choose the best model based on p-value
#'   from likelihood ratio test and difference in AIC. Return the better
#'   reconstruction.
#'
#' @param mat Matrix. Genotype or phenotype binary matrix. (Only discrete
#'   phenotype). Dim: genotype: ncol = ntip(tree) x nrow = number of genotypes.
#'   phenotype: ncol = ntip(tree) x nrow = 1.
#' @param tr Phylo.
#' @param disc_cont String. "discrete" or "continuous."
#' @param num Number. Column index for matrix.
#' @param recon_method String. Either "ML", "REML", "pic", or "GLS".
#'
#' @return Reconstruction of class `ace` from `ape`.
#'
#' @noRd
#'
build_better_reconstruction <- function(mat, tr, disc_cont, num, recon_method){
  # Check input ----------------------------------------------------------------
  if (disc_cont != "discrete") {
    stop("Only pick recon model for discrete. Continuous must use BM.")
  }
  if (num > ncol(mat)) {
    stop("Index must be 1 <= index <= ncol(matrix)")
  }

  # Function -------------------------------------------------------------------
  # Note, SYMreconstruction removed because SYM == ER for binary inputs.
  # Use this function to choose the best model for reconstruction.

  # Cutoffs for comparing the ER and ARD:
  alpha <- 0.05 # For likelihood test
  significant_difference_in_AIC <- 2

  # Reference for model testing:
  # https://www.r-phylo.org/wiki/HowTo/Ancestral_State_Reconstruction &
  # http://blog.phytools.org/2017/07/comparing-fitted-discrete-character.html
  # Test ER vs ARD
  set.seed(1)
  ERreconstruction  <- ape::ace(mat[, num, drop = TRUE],
                           tr,
                           type = disc_cont,
                           method = recon_method,
                           marginal = FALSE,
                           model = "ER")

  # Some ARD models don't work well with the data and given a warning message
  # like:  "In sqrt(diag(solve(h))) : NaNs produced".  To ensure the ER model is
  # preferred in this case use the following warning catching:
  error_msg <- "ARD_bad_fit"
  set.seed(1)
  ARDreconstruction <- tryCatch(ape::ace(mat[, num, drop = TRUE],
                                    tr,
                                    type = disc_cont,
                                    method = recon_method,
                                    marginal = FALSE,
                                    model = "ARD"),
                                warning = function(x) {
    error_msg
    }
  )
  # If ARD gave a warning, pick ER
  best_model <- "ER"
  if (length(ARDreconstruction) == 1) {
    best_model <- "ER"
  } else {
    # Pick best of ER and ARD
    p_ER_ARD <-
      1 - stats::pchisq(
        2 * abs(ERreconstruction$loglik - ARDreconstruction$loglik), 1)
    if (p_ER_ARD < alpha &
        stats::AIC(ERreconstruction) >
        (significant_difference_in_AIC + stats::AIC(ARDreconstruction))) {
      best_model <- "ARD"
    }
  }

  # Given best model, return best ancestral reconstruction
  best_reconstruction <- ERreconstruction
  if (best_model == "ARD") {
    best_reconstruction <- ARDreconstruction
  }

  # Return output --------------------------------------------------------------
  return(best_reconstruction)
} # end build_better_reconstruction()

#' discretize_conf_with_cutoff
#'
#' @description Given a vector with values that describe confidence, binarize
#'  the vector a accoriding to a cutoff value.
#' @param confidence_vector Numeric vector.
#' @param threshold Number.
#'
#' @return Confidence vector. Binary vector.
#'
#' @noRd
#'
discretize_conf_with_cutoff <- function(confidence_vector, threshold){
  # Function -------------------------------------------------------------------
  confidence_vector[confidence_vector  < threshold] <- 0
  confidence_vector[confidence_vector >= threshold] <- 1

  # return output --------------------------------------------------------------
  return(confidence_vector)
} # end discretize_conf_with_cutoff()