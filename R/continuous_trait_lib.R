check_if_phenotype_normal <- function(pheno){
  result <- stats::shapiro.test(unlist(pheno))
  alpha <- 0.05
  is_normal <- TRUE
  if (result$p < alpha) {
    is_normal <- FALSE
  }
  return(is_normal)
}

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
      lamdba_not_close_to_1 <- pheno_not_normal <- TRUE
      
      while (pheno_not_normal) {
        while (lamdba_not_close_to_1) {
          continuous_BM_pheno <- phytools::fastBM(tree = tree_list[[i]])
          # Check that lambda is close to 1 for BM phenotype
          BM_lambda <- phytools::phylosig(tree = tree_list[[i]],
                                          x = continuous_BM_pheno,
                                          method = "lambda")
          lamdba_not_close_to_1 <- BM_lambda$lambda < 0.95 & BM_lambda$lambda > 1.05
        }
        pheno_not_normal <- !(check_if_phenotype_normal(continuous_BM_pheno))
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
    # write.table(continuous_BM_pheno, 
    #             sep = "\t", 
    #             file = paste0("../data/", 
    #                           "continuous_pheno_BM_tree_", 
    #                           i, 
    #                           "_pheno_", 
    #                           j, 
    #                           ".tsv"))
    # write.table(cont_low_signal_pheno, 
    #             sep = "\t", 
    #             file = paste0("../data/", 
    #                           "continuous_pheno_WN_tree_", 
    #                           i, 
    #                           "_pheno_",
    #                           j,
    #                           ".tsv"))
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

                  