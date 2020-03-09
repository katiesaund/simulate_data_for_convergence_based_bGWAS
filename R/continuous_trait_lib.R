make_continuous_phenotypes <- function(tree_list, num_pheno){
  num_trees <- length(tree_list)
  pheno_mat_list_BM <- pheno_mat_list_WN <- rep(list(NULL), num_trees)
  set.seed(1)
  for (i in 1:num_trees) {
    pheno_mat_list_BM[[i]] <- pheno_mat_list_WN[[i]] <- 
      matrix(NA, nrow = ape::Ntip(tree_list[[i]]), ncol = num_pheno)
    colnames(pheno_mat_list_BM[[i]]) <- paste0("BM_pheno_", 1:num_pheno)
    colnames(pheno_mat_list_WN[[i]]) <- paste0("WN_pheno_", 1:num_pheno)
    
    for (j in 1:num_pheno) {
      lamdba_not_close_to_1 <- TRUE
      
      while (lamdba_not_close_to_1) {
        continuous_BM_pheno <- phytools::fastBM(tree = tree_list[[i]])
        # Check that lambda is close to 1 for BM phenotype
        BM_lambda <- phytools::phylosig(tree = tree_list[[i]],
                                        x = continuous_BM_pheno,
                                        method = "lambda")
        lamdba_not_close_to_1 <- BM_lambda$lambda < 0.95 & BM_lambda$lambda > 1.05
      }

      # When lambda is close to zero, the phylogenetic signal is low (White Noise)
      lambda_has_high_signal <- TRUE
      while (lambda_has_high_signal) {
        jumbled_pheno <- sample(unname(continuous_BM_pheno), 
                                size = ape::Ntip(tree_list[[i]]), 
                                replace = FALSE)
        names(jumbled_pheno) <- tree_list[[i]]$tip.label
        jumbled_lambda <- phytools::phylosig(tree = tree_list[[i]],
                                             x = jumbled_pheno,
                                             method = "lambda")
        lambda_has_high_signal <- jumbled_lambda$lambda < -0.05 & jumbled_lambda$lambda > 0.05
      }
      cont_low_signal_pheno <- jumbled_pheno
      
      continuous_BM_pheno <- 
        matrix(c(names(continuous_BM_pheno), continuous_BM_pheno), ncol = 2)  
      row.names(pheno_mat_list_BM[[i]]) <- continuous_BM_pheno[, 1]
      continuous_BM_pheno <- continuous_BM_pheno[, 2, drop = FALSE]
      
      cont_low_signal_pheno <- matrix(c(names(cont_low_signal_pheno), cont_low_signal_pheno), ncol = 2) 
      row.names(pheno_mat_list_WN[[i]]) <- cont_low_signal_pheno[, 1]
      cont_low_signal_pheno <- cont_low_signal_pheno[, 2, drop = FALSE]
      
      pheno_mat_list_BM[[i]][, j] <- continuous_BM_pheno
      pheno_mat_list_WN[[i]][, j] <- cont_low_signal_pheno
    }
  }
  results <- list("cont_pheno_BM_mat_list" = pheno_mat_list_BM, 
                  "cont_pheno_WN_mat_list" = pheno_mat_list_WN)
  return(results)
}

#' Calculate phenotype change per tree edge
#'
#' @description Quantify absoluate value of phenotype change on each tree edge.
#'
#' @param edge_list Numeric vector. Each number is the index of the tree edge to
#'   be used
#' @param phenotype_by_edges Mat. Dimensions: Nedge x 2 matrix. Entries are the
#'   phenotype value at the node, where row is the edge, 1st column is the
#'   parent node and 2nd column is the child node.
#'
#' @return Numeric vector.   Length = length(edge_list).
#' @noRd
calculate_phenotype_change_on_edge <- function(edge_list, phenotype_by_edges){
  # Check input ----------------------------------------------------------------
  if (max(edge_list) > nrow(phenotype_by_edges)) {
    stop("Cannot calculate phenotype change on edges.")
  }
  # check_dimensions(phenotype_by_edges,
  #                  exact_rows = NULL,
  #                  min_rows = max(edge_list),
  #                  exact_cols = NULL,
  #                  min_cols = 2)
  if (!is.vector(edge_list)) {
    stop("edge_list must be a vector of indices of edges")
  }
  # check_is_number(edge_list[1])
  
  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)) {
    delta[j] <-
      abs(phenotype_by_edges[edge_list[j], 1] -
            phenotype_by_edges[edge_list[j], 2])
  }
  
  # Check and return output ----------------------------------------------------
  if (sum(delta < 0) > 0) {
    stop("Delta phenotype should always be recorded as an absolute value.")
  }
  
  return(delta)
}


add_geno_continuous <- function(binary_AR_mat_list, tree_list, num_pheno, cont_pheno_BM_mat_list, cont_pheno_WN_mat_list) {
  num_trees <- length(tree_list)
  for (i in 1:num_trees) {
    num_trait <- ncol(binary_AR_mat_list[[i]])
    num_tip <- ape::Ntip(tree_list[[i]])
    num_node <- ape::Nnode(tree_list[[i]])
    num_to_add <- round(num_trait / 4, 0)
    if (num_to_add < 10) {
      num_to_add <- 10
    }
    num_to_add_per_pheno <- round(num_trait / 10, 0)
    if (num_to_add_per_pheno < 10) {
      num_to_add_per_pheno <- 10
    }
    
    random_col_to_add <- flipped_to_add_10p <- flipped_to_add_25p <- 
      flipped_to_add_40p <- 
      binary_AR_mat_list[[i]][, 1:num_to_add, drop = FALSE]
    for (j in 1:ncol(random_col_to_add)) {
      current_col <- random_col_to_add[1:num_tip, j]
      jumbled_tips <- sample(current_col,
                             size = length(current_col),
                             replace = FALSE)
      random_col_to_add[1:num_tip, j] <- jumbled_tips
      flipped_tip_10p <- flip_some_tips(current_col, flip_freq = 0.10)
      flipped_to_add_10p[1:num_tip, j] <- flipped_tip_10p
      
      flipped_tip_25p <- flip_some_tips(current_col, flip_freq = 0.25)
      flipped_to_add_25p[1:num_tip, j] <- flipped_tip_25p
      
      flipped_tip_40p <- flip_some_tips(current_col, flip_freq = 0.40)
      flipped_to_add_40p[1:num_tip, j] <- flipped_tip_40p
    }
    
    genos_to_add_bc_pheno <- matrix(NA, nrow = num_tip, ncol = 0)
    for (k in 1:num_pheno) {
      # BM
      temp_BM_matrix <- cont_pheno_BM_mat_list[[i]][, k, drop = FALSE]
      geno_like_BM_pheno <- matrix(rep(as.numeric(t(temp_BM_matrix)), each = num_to_add_per_pheno), nrow = nrow(temp_BM_matrix), byrow = TRUE)
      temp_BM_median <- median(as.numeric(temp_BM_matrix[, 1, drop = TRUE]))
      
      # WN
      temp_WN_matrix <- cont_pheno_WN_mat_list[[i]][, k, drop = FALSE]
      geno_like_WN_pheno <- matrix(rep(as.numeric(t(temp_WN_matrix)), each = num_to_add_per_pheno), nrow = nrow(temp_WN_matrix), byrow = TRUE)
      temp_WN_median <- median(as.numeric(temp_WN_matrix[, 1, drop = TRUE]))
      
      geno_like_BM_pheno <- 1 * (geno_like_BM_pheno > temp_BM_median) # 1 times matrix converts it from logical to numeric
      geno_like_WN_pheno <- 1 * (geno_like_WN_pheno > temp_WN_median)
      
      flip_sequence <- seq(from = 0, to = 0.99999, by = 1 / num_to_add_per_pheno)
      for (m in 2:num_to_add_per_pheno) {
        geno_like_BM_pheno[, m] <- flip_some_tips(geno_like_BM_pheno[, m], flip_freq = flip_sequence[m])
        geno_like_WN_pheno[, m] <- flip_some_tips(geno_like_WN_pheno[, m], flip_freq = flip_sequence[m])
      }
      genos_to_add_bc_pheno <- cbind(geno_like_BM_pheno, geno_like_WN_pheno)
    }

    # Add fake ancestral reconstructions
    genos_to_add_bc_pheno <- rbind(genos_to_add_bc_pheno, matrix(0, nrow = num_node, ncol = ncol(genos_to_add_bc_pheno)))
    
    binary_AR_mat_list[[i]] <- cbind(binary_AR_mat_list[[i]],
                                     random_col_to_add,
                                     flipped_to_add_10p, 
                                     flipped_to_add_25p, 
                                     flipped_to_add_40p, 
                                     genos_to_add_bc_pheno)
    
    temp_only_tips <- binary_AR_mat_list[[i]][1:num_tip, , drop = FALSE]
    cols_to_keep <- colSums(temp_only_tips) > 1 & colSums(temp_only_tips) < (nrow(temp_only_tips) - 1)
    binary_AR_mat_list[[i]] <- binary_AR_mat_list[[i]][, cols_to_keep, drop = FALSE]
    # Remove duplicate columns (genotypes)
    binary_AR_mat_list[[i]] <- as.data.frame(t(unique(t(binary_AR_mat_list[[i]])))) 
    
    # Write genotype names (colnames)
    colnames(binary_AR_mat_list[[i]]) <-
      paste0("sim", 1:ncol(binary_AR_mat_list[[i]]))
  }
  return(binary_AR_mat_list)
}

#' Calculate phenotype change per tree edge
#'
#' @description Quantify absoluate value of phenotype change on each tree edge.
#'
#' @param edge_list Numeric vector. Each number is the index of the tree edge to
#'   be used
#' @param phenotype_by_edges Mat. Dimensions: Nedge x 2 matrix. Entries are the
#'   phenotype value at the node, where row is the edge, 1st column is the
#'   parent node and 2nd column is the child node.
#'
#' @return Numeric vector.   Length = length(edge_list).
#' @noRd
calculate_phenotype_change_on_edge <- function(edge_list, phenotype_by_edges){
  # Check input ----------------------------------------------------------------
  if (max(edge_list) > nrow(phenotype_by_edges)) {
    stop("Cannot calculate phenotype change on edges.")
  }
  if (!is.vector(edge_list)) {
    stop("edge_list must be a vector of indices of edges")
  }
  
  # Function -------------------------------------------------------------------
  delta <- rep(NA, length(unique(edge_list)))
  for (j in 1:length(edge_list)) {
    delta[j] <-
      abs(phenotype_by_edges[edge_list[j], 1] -
            phenotype_by_edges[edge_list[j], 2])
  }
  
  # Check and return output ----------------------------------------------------
  if (sum(delta < 0) > 0) {
    stop("Delta phenotype should always be recorded as an absolute value.")
  }
  
  return(delta)
}
                  