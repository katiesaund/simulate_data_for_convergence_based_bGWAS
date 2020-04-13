#' prep_geno_trans_for_phyc
#'
#' @description Discrete testing requires two different definitions of genotype
#'  transition, one for PhyC and one for Synchronous Test. This function
#'  converts geno_trans$transition from the object created for synchronous test
#'  to the version required for the PhyC test.
#'
#'  @details For PhyC prepare genotype transition as below: keep only WT ->
#'   mutant transitions (0 -> 1)
#'
#' @param geno Matrix. Nrow = Ntip(Tree). Ncol = number of genotypes. Binary.
#' @param genotype_transition List of lists. Each sublist has two vectors.
#'  * $transition. Length  == Nedge(tree). 0/1
#'  * $trans_dir. -1/0/1. Length == Nedge(tree).
#'
#' @return genotype_transition. List with $transition and $trans_dir.
prep_geno_trans_for_phyc <- function(geno, genotype_transition){
  # Function -------------------------------------------------------------------
  for (k in 1:ncol(geno)) {
    parent_WT_child_mutant <- 1 # 1 implies parent < child
    parent_mutant_child_WT <- -1 # -1 implies parent > child
    no_transition <- 0 # 0 implies parent == child
    genotype_transition[[k]]$transition <-
      as.numeric(genotype_transition[[k]]$trans_dir == parent_WT_child_mutant)
    genotype_transition[[k]]$trans_dir[genotype_transition[[k]]$trans_dir
                                       == parent_mutant_child_WT] <-
      no_transition # erase transitions in the opposite direction
  }

  # Return output --------------------------------------------------------------
  return(genotype_transition)
} # end prep_geno_trans_for_phyc()

#' is_tip
#'
#' @description Test if a node is a tree tip. An internal node should return
#'  false.
#'
#' @param node_num Integer. Index of node.
#' @param tr Phylo.
#'
#' @return Logical. TRUE OR FALSE.
#' @noRd
is_tip <- function(node_num, tr){
  # Check input ----------------------------------------------------------------
  # check_tree_is_valid(tr)
  # check_is_number(node_num)
  if (node_num < 1 || node_num %% 1 != 0) {
    stop("Node number must be a positive integer")
  }
  #check_node_is_in_tree(node_num, tr)

  # Function & return output ---------------------------------------------------
  return(node_num <= ape::Ntip(tr))
} # end is_tip()

#' identify_transition_edges
#'
#' @description Given a reconstruction identify which edges on tree are
#'  transitions. A transition edge is one in which the parent node and child
#'  node differ.
#'
#' @param tr Phylo.
#' @param mat Matrix. Phenotype (either continuous or binary) or a binary
#'   genotype. Dim: nrow = Ntip(tr) x ncol = {1 if phenotype or number of
#'   genotypes}.
#' @param num Integer. Index of current genotype (column number in genotype matrix).
#' @param node_recon Numeric vector. Either pheno_recon_and_conf$node_anc_rec or
#'  geno_recon_and_conf[[k]]$node_anc_rec. Length = Nnode(tr). Ancestral
#'  reconstruction value for each node.
#' @param disc_cont Character string. Either "discrete" or "continuous".
#'
#' @return List.
#'  * $transition: Continuous: NA, because this is meaningless for continuous data
#'    as all edges will be transitions. Discrete: Numeric vector of 0 or 1. 0
#'    indicates parent and child node are identical. 1 indicates parent and
#'    child node differ.
#'  * $trans_dir: Numeric Vector. Continuous and discrete: gives the direction
#'    of the transition. When parent < child or parent_0_child_1 value is +1.
#'    When parent > child or parent_1_child_0 value is -1.
#' @noRd
identify_transition_edges <- function(tr, mat, num, node_recon, disc_cont){
  # Check input ----------------------------------------------------------------
  # check_for_root_and_bootstrap(tr)
  # check_tree_is_valid(tr)
  # check_is_number(num)
  # check_str_is_discrete_or_continuous(disc_cont)

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
      child_node[i]  <- mat[, num][tr$edge[i, younger]]
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

find_transition_edges <- function(tree_list,
                                  mat_list,
                                  cont_disc_str){
  num_tree <- length(tree_list)
  nested_trans_list <- rep(list(), num_tree)
  for (i in 1:num_tree) {
    num_tip <- ape::Ntip(tree_list[[i]])
    num_col <- ncol(mat_list[[i]])
    num_row <- nrow(mat_list[[i]])
    just_tips_mat <- mat_list[[i]][1:num_tip, , drop = FALSE]
    node_recon <- mat_list[[i]][(num_tip + 1):num_row, , drop = FALSE]
    temp_trans_list <- rep(list(), num_col)
    for (j in 1:num_col) {
      temp_node_recon <- node_recon[, j, drop = TRUE]
      temp_trans_list[[j]] <- identify_transition_edges(tree_list[[i]], just_tips_mat, j, temp_node_recon, cont_disc_str)
    }
    nested_trans_list[[i]] <- temp_trans_list
  }
  return(nested_trans_list)
}

convert_to_phyc_trans <- function(genotype_AR_mat_list,
                                  genotype_sync_trans_list){
  num_mat <- length(genotype_AR_mat_list)
  phyc_trans_list <- rep(list(), num_mat)
  for (i in 1:num_mat) {
    current_trans_list <- genotype_sync_trans_list[[i]]
    num_tip <- ape::Ntip(tree_list[[i]])
    just_tips_mat <- genotype_AR_mat_list[[i]][1:num_tip, , drop = FALSE]
    phyc_trans_list[[i]] <- prep_geno_trans_for_phyc(just_tips_mat, current_trans_list)
  }
  return(phyc_trans_list)
}

#' keep_two_plus_hi_conf_tran_ed
#'
#' @description Since we're looking for convergence of transitions we need a
#'  second quality control step where we remove genotypes that have only 1
#'  transition-edge or where the transition edges are identical!
#'
#' @param genotype_transition List of multiple vectors ($transition and
#'  $trans_dir). Length(list) = number of genotypes. Length(vector) = Nedge(tr).
#' @param genotype_confidence List of vectors. Length(list) = number of
#'  genotypes. Length(vector) = Nedge(tr). Binary.
#'
#' @return at_least_two_hi_conf_trans_ed Logical vector.
#'  Length == length(genotype_transition) == length(genotype_confidence).
#' @noRd
keep_two_plus_hi_conf_tran_ed <- function(genotype_transition,
                                          genotype_confidence){
  # Check inputs ---------------------------------------------------------------
  if (!is.vector(genotype_transition[[1]]$transition)) {
    stop("Input must be a numeric vector")
  }

  # Function -------------------------------------------------------------------
  at_least_two_hi_conf_trans_ed <-
    rep(FALSE, length(genotype_transition))
  for (p in 1:length(genotype_transition)) {
    if (sum(genotype_transition[[p]]$transition *
            genotype_confidence[[p]]) > 1) {
      at_least_two_hi_conf_trans_ed[p] <- TRUE
    }
  }

  # Check and return output ----------------------------------------------------
  return(at_least_two_hi_conf_trans_ed)
}