#' assign_high_confidence_to_transition_edges
#'
#' @description Identify all edges for which the edge is high confidence and a
#'  transition edge.
#'
#' @param tr Phylo.
#' @param all_confidence_by_edge List of vectors. Each vector is binary.
#'  Length(list) == number of genomes.
#' @param geno_trans_by_edge List of vectors. Each vector is binary.
#'  Length(list) == number of genomes.
#' @param geno Matrix. Binary.
#'
#' @return edge_confident_and_trans_edge. List of vector. Each vector is binary.
#'  Length(list) == number of genomes.
#' @noRd
#'
assign_high_confidence_to_transition_edges <- function(all_confidence_by_edge,
                                                       geno_trans_by_edge,
                                                       geno){
  # Function -------------------------------------------------------------------
  edge_confident_and_trans_edge <- rep(list(NULL), ncol(geno))
  for (k in 1:ncol(geno)) {
    edge_confident_and_trans_edge[[k]] <-
      as.numeric( (all_confidence_by_edge[[k]] +
                    geno_trans_by_edge[[k]]$transition) == 2)
  }

  # Return output --------------------------------------------------------------
  return(edge_confident_and_trans_edge)
} 

#' This is a wrapper function for prepare_high_confidence_objects()
#' @noRd
prepare_high_confidence_objects_lists <- function(genotype_trans_by_edge_mat_list,
                                                  tree_list,
                                                  phenotype_conf_by_edge_mat_list,
                                                  boot_threshold,
                                                  geno_AR_mat_list,
                                                  geno_conf_edge_mat_list,
                                                  geno_recon_edge_mat_list,
                                                  snps_in_each_gene = NULL) {
  num_tree <- length(tree_list)
  overall_results <- rep(list(NULL), num_tree)
  for (i in 1:num_tree) {
    num_pheno <- ncol(phenotype_conf_by_edge_mat_list[[i]])
    mini_results <- rep(list(NULL), num_pheno)
    tree <- tree_list[[i]]
    num_tip <- ape::Ntip(tree)
    genotype_transition <- genotype_trans_by_edge_mat_list[[i]]
    geno_mat <- geno_AR_mat_list[[i]][1:num_tip, , drop = FALSE]

    # convert matrices to lists of vectors
    geno_conf_edge <- geno_conf_edge_mat_list[[i]]
    geno_conf_edge <- tapply(geno_conf_edge, 
                             rep(1:ncol(geno_conf_edge),
                                 each = nrow(geno_conf_edge)), 
                             function(i) i)

    geno_recon_edge <- geno_recon_edge_mat_list[[i]]
    geno_recon_edge <- tapply(geno_recon_edge,
                              rep(1:ncol(geno_recon_edge), 
                                  each = nrow(geno_recon_edge)), 
                              function(i) i)

    for (j in 1:num_pheno) {
      temp_pheno_conf_by_edge_vec <-
        phenotype_conf_by_edge_mat_list[[i]][, j, drop = TRUE]
      mini_results[[j]] <- 
        prepare_high_confidence_objects(genotype_transition,
                                        tree,
                                        temp_pheno_conf_by_edge_vec,
                                        boot_threshold,
                                        geno_mat,
                                        geno_conf_edge,
                                        geno_recon_edge,
                                        snps_in_each_gene)
      
    }
    overall_results[[i]] <- mini_results
  }
  return(overall_results)
}

#' prepare_high_confidence_objects
#'
#' @description Identify high confidence edges (considering: tree bootstrap
#'  values, phenotype reconstruction, tree edge lengths, and ancestral
#'  reconstruction of genotype).
#'
#' @param genotype_transition List of lists. Number of lists = number of
#'  genotypes. Each list is made of a $transition and $trans_dir list.
#'  Length(transition) == Length(trans_dir) == Nedge(tree)
#' @param tr Phylo.
#' @param pheno_conf_ordered_by_edges List of confidence values. Binary.
#'  Length(list) == Nedge(tr).
#' @param boot_threshold Numeric. Between 0 and 1.
#' @param geno Matrix. Binary. Nrow = Ntip(tree). Ncol = Number of genotypes.
#' @param geno_conf_edge List of lists. Binary. Number of lists = number of
#'  genotypes. Length(each individual list) == Nedge(Tree)
#' @param geno_recon_edge List of lists. Binary. Number of lists = number of
#'  genotypes. Length(each individual list) == Nedge(Tree)
#' @param snps_in_each_gene Either Null or named table where names are genotypes
#'  and the values are number of not-yet-grouped-genotypes that go into the
#'  grouped genotype.
#'
#' @return List of objects.
#'  * dropped_genotypes Character vector. Names of the genotypes not being kept.
#'  * hi_confidence_transition_edge. List.
#'  * genotype. Matrix.
#'  * snps_per_gene. Either Null or named table where names are genotypes and
#'       the Values are number of not-yet-grouped-genotypes that go into the
#'       grouped genotype.
#'  * genotype_transition Object with two lists: $trans_dir and $transition.
#'  * geno_recon_edge. List of lists. Binary. Number of lists = number of
#'      genotypes. Length(each individual list) == Nedge(tree).
#'  * high_conf_ordered_by_edges. List.
#'  * num_high_conf_trans_edges. Numeric vector. Count of number of
#'      high confidence transitions per genotype. Vector is named with genotype
#'      names.
#'  * tr_and_pheno_hi_conf. Vector. Binary. Length = Nedge(tree).
#'
#' @noRd
#'
prepare_high_confidence_objects <- function(genotype_transition,
                                            tr,
                                            pheno_conf_ordered_by_edges,
                                            boot_threshold,
                                            geno,
                                            geno_conf_edge,
                                            geno_recon_edge,
                                            snps_in_each_gene){
  # Check input ----------------------------------------------------------------
  if (!is.null(snps_in_each_gene)) {
    if (length(snps_in_each_gene) < 1) {
      stop("There should be more grouped genotypes")
    }
  }

  # Function -------------------------------------------------------------------
  tree_conf <- get_bootstrap_confidence(tr, boot_threshold)
  tree_conf_ordered_by_edges <- reorder_tip_and_node_to_edge(tree_conf, tr)
  short_edges <- identify_short_edges(tr)

  high_confidence_edges <-
    pheno_conf_ordered_by_edges + tree_conf_ordered_by_edges + short_edges == 3

  # check_equal(length(high_confidence_edges), ape::Nedge(tr))
  all_high_confidence_edges <- rep(list(0), ncol(geno))

  # ADD IN GENO RECONSTRUCTION CONFIDENCE
  for (k in 1:ncol(geno)) {
    all_high_confidence_edges[[k]] <-
      as.numeric(geno_conf_edge[[k]] + high_confidence_edges == 2)
  }
  only_high_conf_geno_trans <-
    assign_high_confidence_to_transition_edges(all_high_confidence_edges,
                                               genotype_transition,
                                               geno)
  for (i in 1:ncol(geno)) {
    genotype_transition[[i]]$transition <- only_high_conf_geno_trans[[i]]
    genotype_transition[[i]]$trans_dir <-
      only_high_conf_geno_trans[[i]] * genotype_transition[[i]]$trans_dir
  }
  names(only_high_conf_geno_trans) <- colnames(geno)

  # KEEP ONLY GENOTYPES WITH AT LEAST TWO HIGH CONFIDENCE TRANSITION EDGES
  geno_to_keep <- keep_two_plus_hi_conf_tran_ed(genotype_transition,
                                                all_high_confidence_edges)
  genotype_transition <- genotype_transition[geno_to_keep]
  geno_recon_edge <- geno_recon_edge[geno_to_keep]
  high_conf_ordered_by_edges <- all_high_confidence_edges[geno_to_keep]
  dropped_genotypes <- get_dropped_genotypes(geno, geno_to_keep)
  geno <- geno[, geno_to_keep, drop = FALSE]
  snps_in_each_gene <-
    snps_in_each_gene[names(snps_in_each_gene) %in% colnames(geno)]

  names(genotype_transition) <- names(geno_recon_edge) <- colnames(geno)

  # Check output ---------------------------------------------------------------
  if (length(genotype_transition) == 0) {
    stop("No genotypes to test because all genotypes failed quality control")
  }

  # Return output --------------------------------------------------------------
  results <-
    list("dropped_genotypes" = dropped_genotypes,
          "hi_confidence_transition_edge" = only_high_conf_geno_trans,
          "genotype" = geno,
          "snps_per_gene" = snps_in_each_gene,
          "genotype_transition" = genotype_transition,
          "geno_recon_edge" = geno_recon_edge,
          "high_conf_ordered_by_edges" = high_conf_ordered_by_edges,
          "tr_and_pheno_hi_conf" = high_confidence_edges)
  return(results)
} # end prepare_high_confidence_objects()

#' Wrapper function for reorder_tip_and_node_to_edge() to work on lists
reorder_tip_and_node_to_edge_lists <- function(geno_conf_mat_list, tree_list) {
  num_mat <- length(geno_conf_mat_list)
  geno_conf_by_edges_mat_list <- rep(list(0), num_mat)
  for (i in 1:num_mat) {
    tree <- tree_list[[i]]
    num_edge <- ape::Nedge(tree)
    num_geno <- ncol(geno_conf_mat_list[[i]])
    geno_conf_by_edges_mat_list[[i]] <- 
      matrix(NA, nrow = num_edge, ncol = num_geno)
    for (k in 1:num_geno) {
      temp_column <- 
        reorder_tip_and_node_to_edge(unname(geno_conf_mat_list[[i]][, k, drop = TRUE]),
                                     tree)
      geno_conf_by_edges_mat_list[[i]][, k] <- temp_column
    }
    colnames(geno_conf_by_edges_mat_list[[i]]) <- colnames(geno_conf_mat_list[[i]])
    row.names(geno_conf_by_edges_mat_list[[i]]) <- 
      paste0("edge", 1:nrow(geno_conf_by_edges_mat_list[[i]]))
  }
  return(geno_conf_by_edges_mat_list)
}

#' Get names of genotypes dropped from analysis
#'
#' @description  Given a matrix with both desirable genotypes (keepers) and
#'  undesirable genotypes and a numeric index of the keepers, get the names of
#'  the undesirable genotypes.
#'
#' @param geno Matrix. Numeric, binary genotype matrix. Columns = genotypes.
#'   Rows = samples.
#' @param keepers Logical vector. Length == length(genotype_transition)  ==
#'  length(genotype_confidence) == Number of genotypes. True indicates the
#'  genotype has at least 2 high confidence genotype transition edges. False
#'  indicates genotype lacks at least 2 high confidence genotype transition
#'  edges.
#'
#' @return dropped_genotype_names. Character vector. Names of the non-keepers
#'  (genotypes not to be processed downstream).
#'
#' @noRd
#'
get_dropped_genotypes <- function(geno, keepers){
  # Check input ----------------------------------------------------------------
  if (ncol(geno) != length(keepers)) {
    stop("Keepers must have an entry for each genotype.")
  }
  if (sum(keepers) == 0) {
    stop("There are no genotypes left to test.")
  }

  # Function -------------------------------------------------------------------
  dropped_genotype_names <- colnames(geno)[!keepers]

  # Return output --------------------------------------------------------------
  return(dropped_genotype_names)
} # end get_dropped_genotypes()

list_to_matrix_columns <- function(lst) {
  num_col <- length(lst)
  num_row <- length(lst[[1]])
  mat <- matrix(NA, ncol = num_col, nrow = num_row)
  for (i in 1:num_col) {
    mat[, i] <- lst[[i]]
  }
  return(mat)
}

convert_conf_obj_to_conf_mat <- function(conf_obj) {
  num_tree <- length(conf_obj)
  big_conf_mat_list <- rep(list(NULL), num_tree)
  for (i in 1:num_tree) {
    num_pheno <- length(conf_obj[[i]])
    conf_mat <- rep(list(NULL), num_pheno)
    for (j in 1:num_pheno) {
      conf_mat[[j]] <-
        list_to_matrix_columns(conf_obj[[i]][[j]]$high_conf_ordered_by_edges)
    }
    big_conf_mat_list[[i]] <- conf_mat
  }
  return(big_conf_mat_list)
}
