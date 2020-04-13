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


add_geno_continuous <- function(binary_AR_mat_list,
                                tree_list, 
                                num_pheno, 
                                cont_pheno_BM_mat_list, 
                                cont_pheno_WN_mat_list, 
                                BM_phenotype_recon_edge_mat_list, 
                                WN_phenotype_recon_edge_mat_list) {
  num_trees <- length(tree_list)
  for (i in 1:num_trees) {
    num_trait <- ncol(binary_AR_mat_list[[i]])
    num_tip <- ape::Ntip(tree_list[[i]])
    num_node <- ape::Nnode(tree_list[[i]])
    num_edge <- ape::Nedge(tree_list[[i]])
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
    
    # 1. Given a simulated genotype, jumble up the tips to various degrees
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
    
    # 2. Make genotypes based on phenotype ranking
    rank_based_genos <- matrix(NA, nrow = num_tip, ncol = 0)
    for (k in 1:num_pheno) {
      # BM
      temp_BM_matrix <- cont_pheno_BM_mat_list[[i]][, k, drop = FALSE]
      geno_like_BM_pheno <- make_rank_based_geno_mat(temp_BM_matrix, tree_list[[i]])
        
      # WN
      temp_WN_matrix <- cont_pheno_WN_mat_list[[i]][, k, drop = FALSE]
      geno_like_WN_pheno <- make_rank_based_geno_mat(temp_WN_matrix, tree_list[[i]])
      
      rank_based_genos <- cbind(geno_like_BM_pheno, geno_like_WN_pheno)
    }

    # Add fake ancestral reconstructions to genotypes based on phenotype ranking
    rank_based_genos <- rbind(rank_based_genos, 
                                   matrix(0, nrow = num_node, ncol = ncol(rank_based_genos)))
    
    # 3. Add genotypes loosely built based on which edges the phenotype changes lot
    #    A) Given delta phenotype for each edge, identify the edges with high delta (3rd quartile)
    #    B) Starting from a root value of 0 (artibrary) move down each clade and flip all child values when reaching a high delta edge
    #    C) return genotype at just the tips once this procedure finishes
    #    D) it won't be perfect b/c this isn't a reverse engineering of the ancesral reconstruction process
    
    pheno_delta_based_genos <- matrix(NA, nrow = num_tip, ncol = 0)
    for (k in 1:num_pheno) {
      # Calculate absolute value of child - parent node value
      pheno_delta_BM <-
        calculate_phenotype_change_on_edge(1:num_edge, 
                                           BM_phenotype_recon_edge_mat_list[[i]][[k]])
      
      pheno_delta_WN <-
        calculate_phenotype_change_on_edge(1:num_edge, 
                                           WN_phenotype_recon_edge_mat_list[[i]][[k]])
      
      BM_based_geno <- construct_geno_from_pheno(tree_list[[i]], pheno_delta_BM)
      WN_based_geno <- construct_geno_from_pheno(tree_list[[i]], pheno_delta_WN)
      
      pheno_delta_based_genos <- cbind(BM_based_geno, WN_based_geno)
    }
    
    # Add fake ancestral reconstructions to genotypes based on phenotype deltas
    pheno_delta_based_genos <- rbind(pheno_delta_based_genos, 
                              matrix(0, nrow = num_node, ncol = ncol(pheno_delta_based_genos)))
    
    binary_AR_mat_list[[i]] <- cbind(binary_AR_mat_list[[i]],
                                     random_col_to_add,
                                     flipped_to_add_10p, 
                                     flipped_to_add_25p, 
                                     flipped_to_add_40p, 
                                     rank_based_genos, 
                                     pheno_delta_based_genos)
    
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
                  
# For constructing a genotype from root to tip, based on phenotype changes
construct_geno_from_pheno <- function(tree, delta_pheno_vec) {
  tips <- ape::Ntip(tree)
  num_edge <- ape::Nedge(tree)
  
  if (length(delta_pheno_vec) != num_edge) { stop("delta pheno vector needs entry for each edge")}
  
  third_quartile_value <- summary(delta_pheno_vec)[5] # number
  hi_delta_edge_log <- delta_pheno_vec > third_quartile_value # vector of T/F length == nedge(Tree)
  fake_geno_edge_mat <- matrix(0, nrow = nrow(tree$edge), ncol = 2)
  for (i in 1:nrow(tree$edge)) {
    if (hi_delta_edge_log[i]) {
      child_node <- tree$edge[i, 2]
      child_node_and_tips <- phytools::getDescendants(tree, child_node)
      child_node_and_tips <- c(child_node, child_node_and_tips)
      new_value <- as.numeric(!fake_geno_edge_mat[i, 1])
      fake_geno_edge_mat[tree$edge %in% child_node_and_tips] <- new_value
    }
  }
  
  # Now get tip values for the tree in tree tip order
  geno_at_tips <- matrix(NA, nrow = tips, ncol = 1)
  row.names(geno_at_tips) <- tree$tip.label
  for (i in 1:tips) {
    tip_id <-  strsplit(tree$tip.label, "")[[i]][2] # if tip label encoded as "t1" for example
    tip_value <- fake_geno_edge_mat[, 2][tree$edge[, 2] == tip_id]
    geno_at_tips[i, 1] <- tip_value
  }
  return(geno_at_tips)
}

# make a rank based genotype (sweep all possible values divisions based on rank)
make_rank_based_geno_mat <- function(pheno_vec, tree) {
  num_tips <- ape::Ntip(tree)
  geno_mat <- matrix(rank(pheno_vec), nrow = num_tips, ncol = 2 * num_tips)
  for (i in 1:num_tips) {
    geno_mat[, i] <- as.numeric(geno_mat[, i] < i)
  }
  
  for (i in (num_tips + 1):(2 * num_tips)) {
    geno_mat[, i] <- as.numeric(geno_mat[, i] > (i - num_tips))
  }
  return(geno_mat)
}