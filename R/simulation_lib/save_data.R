# Save all of the data necessary for running hogwash on simulated data
save_data <- function(tree_list,
                      genotype_AR_mat_list,
                      genotype_phyc_trans_list,
                      genotype_sync_trans_list,
                      BM_phenotype_AR_mat_list,
                      BM_pheno_recon_by_edge_list,
                      BM_phyc_convergence_list,
                      BM_sync_convergence_list,
                      WN_phenotype_AR_mat_list,
                      WN_pheno_recon_by_edge_list,
                      WN_phyc_convergence_list,
                      WN_sync_convergence_list,
                      BM_phyc_genotype_keeper_list,
                      BM_sync_genotype_keeper_list,
                      WN_phyc_genotype_keeper_list,
                      WN_sync_genotype_keeper_list) {

  num_trees <- length(tree_list)
  summary_col_names <- c("phenotype_phylogenetic_signal",
                         "tree_id",
                         "phenotype_id",
                         "test",
                         "intersection",
                         "num_hi_conf_edges",
                         "pheno_beta",
                         "geno_beta",
                         "epsilon", 
                         "genotype")
  num_cols <- length(summary_col_names)
  summary_df <- data.frame(matrix(NA, nrow = 1 , ncol = num_cols))
  colnames(summary_df) <- summary_col_names
  for (i in 1:num_trees) {
    current_tree <- tree_list[[i]]
    num_tips <- ape::Ntip(current_tree)
    num_pheno <- ncol(BM_phenotype_AR_mat_list[[i]])


    ape::write.tree(current_tree,
                    file = paste0("../data/",
                                  "simulated_discrete_tree_",
                                  i,
                                  ".tree"))

    for (j in 1:num_pheno) {
      # Save genotypes ----
      BM_keepers <- sort(unique(c(BM_phyc_genotype_keeper_list[[i]][[j]],
                                  BM_sync_genotype_keeper_list[[i]][[j]])))
      WN_keepers <- sort(unique(c(WN_phyc_genotype_keeper_list[[i]][[j]],
                                  WN_sync_genotype_keeper_list[[i]][[j]])))

      BM_geno_mat <- 
        genotype_AR_mat_list[[i]][1:num_tips, BM_keepers, drop = FALSE]
      WN_geno_mat <- 
        genotype_AR_mat_list[[i]][1:num_tips, WN_keepers, drop = FALSE]

      write.table(BM_geno_mat,
                  sep = "\t",
                  row.names = TRUE,
                  file = paste0("../data/",
                                "simulated_genotype_for_discrete_pheno_BM_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))

      write.table(WN_geno_mat,
                  sep = "\t",
                  row.names = TRUE,
                  file = paste0("../data/",
                                "simulated_genotype_for_discrete_pheno_WN_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))

      # Save phenotypes ----
      temp_disc_BM_trait <-
        BM_phenotype_AR_mat_list[[i]][1:num_tips, j, drop = FALSE]
      write.table(temp_disc_BM_trait,
                  sep = "\t",
                  file = paste0("../data/",
                                "simulated_discrete_pheno_BM_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))

      temp_disc_WN_trait <-
        WN_phenotype_AR_mat_list[[i]][1:num_tips, j, drop = FALSE]
      write.table(temp_disc_WN_trait,
                  sep = "\t",
                  file = paste0("../data/",
                                "simulated_discrete_pheno_WN_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))

      # Save convergence results ----
      BM_phyc_convergence <- BM_phyc_convergence_list[[i]][[j]]
      save(BM_phyc_convergence, file = paste0("../data/",
                                        "simulated_discrete_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_phyc_convergence.RData"))
      BM_sync_convergence <- BM_sync_convergence_list[[i]][[j]]
      save(BM_sync_convergence, file = paste0("../data/",
                                        "simulated_discrete_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_sync_convergence.RData"))

      WN_phyc_convergence <- WN_phyc_convergence_list[[i]][[j]]
      save(WN_phyc_convergence, file = paste0("../data/",
                                        "simulated_discrete_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_phyc_convergence.RData"))
      WN_sync_convergence <- WN_sync_convergence_list[[i]][[j]]
      save(WN_sync_convergence, file = paste0("../data/",
                                        "simulated_discrete_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_sync_convergence.RData"))

      BM_phyc_convergence_df <- 
        as.data.frame(matrix(NA, 
                             nrow = length(BM_phyc_convergence$intersection),
                             ncol = num_cols))
      BM_sync_convergence_df <- 
        as.data.frame(matrix(NA, 
                             nrow = length(BM_sync_convergence$intersection), 
                             ncol = num_cols))
      WN_phyc_convergence_df <- 
        as.data.frame(matrix(NA, 
                             nrow = length(WN_phyc_convergence$intersection), 
                             ncol = num_cols))
      WN_sync_convergence_df <-
        as.data.frame(matrix(NA, 
                             nrow = length(WN_sync_convergence$intersection), 
                             ncol = num_cols))

      colnames(BM_phyc_convergence_df) <- colnames(BM_sync_convergence_df) <-
        colnames(WN_phyc_convergence_df) <- colnames(WN_sync_convergence_df) <-
        summary_col_names

      BM_phyc_convergence_df$phenotype_phylogenetic_signal <- 
        BM_sync_convergence_df$phenotype_phylogenetic_signal <- "BM"
      WN_phyc_convergence_df$phenotype_phylogenetic_signal <-
        WN_sync_convergence_df$phenotype_phylogenetic_signal <- "WN"

      BM_phyc_convergence_df$test <- WN_phyc_convergence_df$test <- "phyc"
      BM_sync_convergence_df$test <- WN_sync_convergence_df$test <- "sync"

      BM_phyc_convergence_df$tree_id <- i
      BM_phyc_convergence_df$phenotype_id <- j
      BM_phyc_convergence_df$intersection <- BM_phyc_convergence$intersection
      BM_phyc_convergence_df$num_hi_conf_edges <- BM_phyc_convergence$num_hi_conf_edges
      BM_phyc_convergence_df$pheno_beta <- BM_phyc_convergence$pheno_beta
      BM_phyc_convergence_df$geno_beta <- BM_phyc_convergence$geno_beta
      BM_phyc_convergence_df$epsilon <- BM_phyc_convergence$epsilon
      BM_phyc_convergence_df$genotype <- BM_phyc_convergence$genotype
      
      BM_sync_convergence_df$tree_id <- i
      BM_sync_convergence_df$phenotype_id <- j
      BM_sync_convergence_df$intersection <- BM_sync_convergence$intersection
      BM_sync_convergence_df$num_hi_conf_edges <- BM_sync_convergence$num_hi_conf_edges
      BM_sync_convergence_df$pheno_beta <- BM_sync_convergence$pheno_beta
      BM_sync_convergence_df$geno_beta <- BM_sync_convergence$geno_beta
      BM_sync_convergence_df$epsilon <- BM_sync_convergence$epsilon
      BM_sync_convergence_df$genotype <- BM_sync_convergence$genotype
      
      WN_phyc_convergence_df$tree_id <- i
      WN_phyc_convergence_df$phenotype_id <- j
      WN_phyc_convergence_df$intersection <- WN_phyc_convergence$intersection
      WN_phyc_convergence_df$num_hi_conf_edges <- WN_phyc_convergence$num_hi_conf_edges
      WN_phyc_convergence_df$pheno_beta <- WN_phyc_convergence$pheno_beta
      WN_phyc_convergence_df$geno_beta <- WN_phyc_convergence$geno_beta
      WN_phyc_convergence_df$epsilon <- WN_phyc_convergence$epsilon
      WN_phyc_convergence_df$genotype <- WN_phyc_convergence$genotype
      
      WN_sync_convergence_df$tree_id <- i
      WN_sync_convergence_df$phenotype_id <- j
      WN_sync_convergence_df$intersection <- WN_sync_convergence$intersection
      WN_sync_convergence_df$num_hi_conf_edges <- WN_sync_convergence$num_hi_conf_edges
      WN_sync_convergence_df$pheno_beta <- WN_sync_convergence$pheno_beta
      WN_sync_convergence_df$geno_beta <- WN_sync_convergence$geno_beta
      WN_sync_convergence_df$epsilon <- WN_sync_convergence$epsilon
      WN_sync_convergence_df$genotype <- WN_sync_convergence$genotype
      
      summary_df <- rbind(summary_df,
                          BM_phyc_convergence_df, 
                          BM_sync_convergence_df, 
                          WN_phyc_convergence_df,
                          WN_sync_convergence_df)
    }
  }

  summary_df <- summary_df[2:nrow(summary_df), , drop = FALSE] # remove the NA row
  save(summary_df, file = "../data/simulated_discrete_convergence_summary.RData")
  write.table(summary_df,
              sep = "\t",
              file = "../data/simulated_discrete_convergence_summary.tsv", 
              row.names = FALSE)
}

save_continuous_data <- function(tree_list,
                                 genotype_AR_mat_list,
                                 genotype_cont_trans_list,
                                 BM_phenotype_AR_mat_list,
                                 BM_pheno_recon_by_edge_list,
                                 BM_cont_convergence_list,
                                 WN_phenotype_AR_mat_list,
                                 WN_pheno_recon_by_edge_list,
                                 WN_cont_convergence_list,
                                 BM_cont_genotype_keeper_list,
                                 WN_cont_genotype_keeper_list) {
  
  num_trees <- length(tree_list)
  summary_col_names <- c("phenotype_phylogenetic_signal",
                         "tree_id",
                         "phenotype_id",
                         "test",
                         "intersection",
                         "num_hi_conf_edges",
                         "pheno_beta",
                         "geno_beta",
                         "epsilon", 
                         "genotype")
  num_cols <- length(summary_col_names)
  summary_df <- data.frame(matrix(NA, nrow = 1 , ncol = num_cols))
  colnames(summary_df) <- summary_col_names
  for (i in 1:num_trees) {
    current_tree <- tree_list[[i]]
    num_tips <- ape::Ntip(current_tree)
    num_pheno <- ncol(BM_phenotype_AR_mat_list[[i]])
    ape::write.tree(current_tree,
                    file = paste0("../data/",
                                  "simulated_continuous_tree_",
                                  i,
                                  ".tree"))
    
    for (j in 1:num_pheno) {
      # Save genotypes ----
      BM_keepers <- BM_cont_genotype_keeper_list[[i]][[j]]
      WN_keepers <- WN_cont_genotype_keeper_list[[i]][[j]]
      
      BM_geno_mat <- genotype_AR_mat_list[[i]][1:num_tips,
                                               BM_keepers,
                                               drop = FALSE]
      WN_geno_mat <- genotype_AR_mat_list[[i]][1:num_tips,
                                               WN_keepers,
                                               drop = FALSE]
      
      write.table(BM_geno_mat,
                  sep = "\t",
                  row.names = TRUE,
                  file = paste0("../data/",
                                "simulated_genotype_for_continuous_pheno_BM_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))
      
      write.table(WN_geno_mat,
                  sep = "\t",
                  row.names = TRUE,
                  file = paste0("../data/",
                                "simulated_genotype_for_continuous_pheno_WN_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))
      
      # Save phenotypes ----
      temp_cont_BM_trait <-
        BM_phenotype_AR_mat_list[[i]][1:num_tips, j, drop = FALSE]
      write.table(temp_cont_BM_trait,
                  sep = "\t",
                  file = paste0("../data/",
                                "simulated_continuous_pheno_BM_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))
      
      temp_cont_WN_trait <-
        WN_phenotype_AR_mat_list[[i]][1:num_tips, j, drop = FALSE]
      write.table(temp_cont_WN_trait,
                  sep = "\t",
                  file = paste0("../data/",
                                "simulated_continuous_pheno_WN_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))
      
      # Save convergence results ----
      BM_cont_convergence <- BM_cont_convergence_list[[i]][[j]]
      save(BM_cont_convergence, file = paste0("../data/",
                                        "simulated_continuous_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_cont_convergence.RData"))
      

      WN_cont_convergence <- WN_cont_convergence_list[[i]][[j]]
      save(WN_cont_convergence, file = paste0("../data/",
                                        "simulated_continuous_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_cont_convergence.RData"))
      
      BM_cont_convergence_df <- 
        as.data.frame(matrix(NA,
                             nrow = length(BM_cont_convergence$intersection), 
                             ncol = num_cols))
      WN_cont_convergence_df <-
        as.data.frame(matrix(NA, 
                             nrow = length(WN_cont_convergence$intersection),
                             ncol = num_cols))
      
      colnames(BM_cont_convergence_df) <- 
        colnames(WN_cont_convergence_df) <- 
        summary_col_names
      
      BM_cont_convergence_df$test <- WN_cont_convergence_df$test <- "cont"
      
      BM_cont_convergence_df$tree_id <- i
      BM_cont_convergence_df$phenotype_id <- j
      BM_cont_convergence_df$intersection <- BM_cont_convergence$intersection
      BM_cont_convergence_df$num_hi_conf_edges <- BM_cont_convergence$num_hi_conf_edges
      BM_cont_convergence_df$pheno_beta <- BM_cont_convergence$pheno_beta
      BM_cont_convergence_df$geno_beta <- BM_cont_convergence$geno_beta
      BM_cont_convergence_df$epsilon <- BM_cont_convergence$epsilon
      BM_cont_convergence_df$genotype <- BM_cont_convergence$genotype
      
      WN_cont_convergence_df$tree_id <- i
      WN_cont_convergence_df$phenotype_id <- j
      WN_cont_convergence_df$intersection <- WN_cont_convergence$intersection
      WN_cont_convergence_df$num_hi_conf_edges <- WN_cont_convergence$num_hi_conf_edges
      WN_cont_convergence_df$pheno_beta <- WN_cont_convergence$pheno_beta
      WN_cont_convergence_df$geno_beta <- WN_cont_convergence$geno_beta
      WN_cont_convergence_df$epsilon <- WN_cont_convergence$epsilon
      WN_cont_convergence_df$genotype <- WN_cont_convergence$genotype
      
      summary_df <- rbind(summary_df, BM_cont_convergence_df, WN_cont_convergence_df)
    }
  }
  
  summary_df <- summary_df[2:nrow(summary_df), , drop = FALSE] # remove the NA row
  save(summary_df, file = "../data/simulated_continuous_convergence_summary.RData")
  write.table(summary_df,
              sep = "\t",
              file = "../data/simulated_continuous_convergence_summary.tsv", 
              row.names = FALSE)
}