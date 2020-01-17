generate_trees <- function(num_trees, num_tips, edge_multiplier){
  # Generate a series of trees with a set number of tips
  simulated_trees <- rep(list(0), num_trees)
  set.seed(1)
  for (i in 1:num_trees) {
    # Use rcoal instead of rtree to get realistic branch lengths.
    simulated_trees[[i]] <- ape::rcoal(n = num_tips, rooted = TRUE)
    simulated_trees[[i]]$node.label <- rep(100, Nnode(simulated_trees[[i]]))
    simulated_trees[[i]]$edge.length <-
      edge_multiplier * simulated_trees[[i]]$edge.length
  }
  return(simulated_trees)
}

#' reorder_tip_and_node_to_edge
#'
#' @description Reorder a vector organzied by tips then nodes into a vector
#'  organized by tree edges.
#'
#'  @details This function grabs the value of each child from every node and
#'   stores them in the order of the tree's edge matrix. This effectively
#'   drops the value of the tree's root because the root is never a child and
#'   therefore not in the child column of the edge matrix. Each node/tip that
#'   is not the root is a child exactly once.
#'
#' @param tips_and_node_vector Numeric vector. Length = Ntip(tr) + Nnode(tr).
#' @param tr Phylo.
#'
#' @return ordered_by_edges. Numeric vector. Length = Nedge(tr).
#' @noRd
reorder_tip_and_node_to_edge <- function(tips_and_node_vector, tr){
  # Check input ----------------------------------------------------------------
  # check_tree_is_valid(tr)
  # check_for_root_and_bootstrap(tr)
  # check_equal(length(tips_and_node_vector), sum(ape::Ntip(tr), ape::Nnode(tr)))
  # check_is_number(tips_and_node_vector[1])

  # Function -------------------------------------------------------------------
  ordered_by_edges <- rep(NA, ape::Nedge(tr))
  for (i in 1:ape::Nedge(tr)) {
    ordered_by_edges[i] <- tips_and_node_vector[tr$edge[i, 2]]
  }

  # Return output --------------------------------------------------------------
  return(ordered_by_edges)
}# end reorder_tip_and_node_to_edge()

prep_pheno_recon_edges <- function(phenotype_AR_mat_list, tree_list) {
  num_trees <- length(tree_list)
  pheno_recon_by_edges_list <- list()
  for (i in 1:num_trees) {
    num_pheno <- ncol(phenotype_AR_mat_list[[i]])
    temp_recon_by_edges <- list()
    for (j in 1:num_pheno) {
      temp_recon_by_edges[[j]] <- reorder_tip_and_node_to_edge(phenotype_AR_mat_list[[i]][, j, drop = TRUE], tree_list[[i]])
    }
    pheno_recon_by_edges_list[[i]] <- temp_recon_by_edges
  }
  return(pheno_recon_by_edges_list)
}


