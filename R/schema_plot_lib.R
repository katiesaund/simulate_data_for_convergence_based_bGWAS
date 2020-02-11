library(ape)
library(phytools)
#' Check that the given node if in the tree
#'
#' @description Test if a node value is plausibly contained within the tree.
#'
#' @param node_val Integer. Index of node.
#' @param tr Phylo.
#'
#' @noRd
check_node_is_in_tree <- function(node_val, tr){
  # Check input & function -----------------------------------------------------
  if (node_val > ape::Nnode(tr) + ape::Ntip(tr)) {
    stop("Node number is too high; not found in tree.")
  }
  if (node_val < 1 | node_val %% 1 != 0) {
    stop("Node must be positive integer")
  }
}

#' Test if a node is also a tip
#'
#' @description Test if a node is a tree tip. An internal node should return
#'   false.
#'
#' @param node_num Integer. Index of node.
#' @param tr Phylo.
#'
#' @return Logical. TRUE OR FALSE.
#' @noRd
is_tip <- function(node_num, tr){
  # Check input ----------------------------------------------------------------
  if (node_num < 1 || node_num %% 1 != 0) {
    stop("Node number must be a positive integer")
  }
  check_node_is_in_tree(node_num, tr)
  
  # Function & return output ---------------------------------------------------
  return(node_num <= ape::Ntip(tr))
}

reorder_tip_and_node_to_edge <- function(tips_and_node_vector, tr){
  # Function -------------------------------------------------------------------
  ordered_by_edges <- rep(NA, ape::Nedge(tr))
  for (i in 1:ape::Nedge(tr)) {
    ordered_by_edges[i] <- tips_and_node_vector[tr$edge[i, 2]]
  }
  
  # Return output --------------------------------------------------------------
  return(ordered_by_edges)
}

#' Identify transition edges on a tree
#'
#' @description Given a reconstruction identify which edges on tree are
#'   transitions. A transition edge is one in which the parent node and child
#'   node differ.
#'
#' @param tr Phylo.
#' @param vec Vector of the trait
#' @param node_recon Numeric vector. Length = Nnode(tr). Ancestral
#'   reconstruction value for each node.
#' @param disc_cont Character string. Either "discrete" or "continuous".
#'
#' @return List.
#'   \describe{
#'    \item{transition}{Continuous: NA, because this is meaningless for
#'    continuous data as all edges will be transitions. Discrete: Numeric vector
#'    of 0 or 1. 0 indicates parent and child node are identical. 1 indicates
#'    parent and child node differ.}
#'    \item{trans_dir}{Numeric Vector. Continuous and discrete: gives the
#'    direction of the transition. When parent < child or parent_0_child_1 value
#'    is +1. When parent > child or parent_1_child_0 value is -1.}
#'   }
#' @noRd
identify_transition_edges <- function(tr, vec, node_recon, disc_cont){
  # FUNCTION -------------------------------------------------------------------
  transition <- transition_direction <-
    parent_node <- child_node <- integer(ape::Nedge(tr))
  older <- 1 # older node is 1st column in tr$edge
  younger <- 2 # younger node is 2nd column in tr$edge
  parent_0_child_1 <- 1
  parent_1_child_0 <- -1
  parent_equals_child <- 0
  both_parent_and_child_are_one <- 2
  
  
  for (i in 1:ape::Nedge(tr)) {
    if (is_tip(tr$edge[i, older], tr)) {
      stop("tree invalid")
    }
    parent_node[i] <- node_recon[tr$edge[i, older] - ape::Ntip(tr)]
    if (is_tip(tr$edge[i, younger], tr)) {
      # child is a tip
      child_node[i]  <- vec[tr$edge[i, younger]]
    } else {
      # child is internal nodes
      child_node[i]  <- node_recon[tr$edge[i, younger] - ape::Ntip(tr)]
    }
    
    transition[i] <- sum(parent_node[i] + child_node[i])
    # transition[i] is either 0, 1, or 2 for discrete traits
    # transition[i] is not to be used when the trait is continuous because all,
    # or very nearly all edges are transition edges.
    
    if (parent_node[i] > child_node[i]) {
      transition_direction[i] <- parent_1_child_0
    } else if (parent_node[i] < child_node[i]) {
      transition_direction[i] <- parent_0_child_1
    }
    
    if (disc_cont == "discrete") {
      transition[transition == both_parent_and_child_are_one] <-
        parent_equals_child # parent_node == child_node, then no transition.
    } else {
      transition <- NA # Not used when the trait is continuous.
    }
  }
  # Check and return output ----------------------------------------------------
  results <- list("transition" = transition, "trans_dir" = transition_direction)
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


#' Convert a vector to an edge matrix
#'
#' @description Convert the reconstruction to be in the same format as
#'   tree$edge, where each row is an edge. This format will make it much easier
#'   to calculate the phenotype change on each edge for continuous phenotypes.
#'
#' @param tr Phylo.
#' @param tip_and_node_reconstruction Vector. Ordered by tips then by nodes.
#'   Length = Ntip(tr) + Nnode(tr). Observed values at each tip and ancestral
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
