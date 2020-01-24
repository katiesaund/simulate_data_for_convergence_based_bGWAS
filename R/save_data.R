save_data <- function(tree_list,
                      genotype_AR_mat_list,
                      genotype_phyc_trans_list,
                      genotype_sync_trans_list,
                      BM_phenotype_AR_mat_list,
                      BM_pheno_recon_by_edge_list,
                      BM_phyc_gamma_list,
                      BM_sync_gamma_list,
                      WN_phenotype_AR_mat_list,
                      WN_pheno_recon_by_edge_list,
                      WN_phyc_gamma_list,
                      WN_sync_gamma_list,
                      BM_phyc_genotype_keeper_list,
                      BM_sync_genotype_keeper_list,
                      WN_phyc_genotype_keeper_list,
                      WN_sync_genotype_keeper_list) {

  num_trees <- length(tree_list)
  summary_col_names <- c("phenotype_phylogenetic_signal",
                         "tree_id",
                         "phenotype_id",
                         "test",
                         "gamma_avg",
                         "gamma_percent",
                         "gamma_count",
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


    write.tree(current_tree,
               file = paste0("../data/",
                             "simulated_tree_",
                             i,
                             ".tree"))

    for (j in 1:num_pheno) {
      # Save genotypes ----
      BM_keepers <- sort(unique(c(BM_phyc_genotype_keeper_list[[i]][[j]],
                                  BM_sync_genotype_keeper_list[[i]][[j]])))
      WN_keepers <- sort(unique(c(WN_phyc_genotype_keeper_list[[i]][[j]],
                                  WN_sync_genotype_keeper_list[[i]][[j]])))

      BM_geno_mat <- genotype_AR_mat_list[[i]][1:num_tips, BM_keepers, drop = FALSE]
      WN_geno_mat <- genotype_AR_mat_list[[i]][1:num_tips, WN_keepers, drop = FALSE]

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

      # Save gamma results ----
      BM_phyc_gamma <- BM_phyc_gamma_list[[i]][[j]]
      save(BM_phyc_gamma, file = paste0("../data/",
                                        "simulated_discrete_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_phyc_gamma.RData"))
      BM_sync_gamma <- BM_sync_gamma_list[[i]][[j]]
      save(BM_sync_gamma, file = paste0("../data/",
                                        "simulated_discrete_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_sync_gamma.RData"))

      WN_phyc_gamma <- WN_phyc_gamma_list[[i]][[j]]
      save(WN_phyc_gamma, file = paste0("../data/",
                                        "simulated_discrete_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_phyc_gamma.RData"))
      WN_sync_gamma <- WN_sync_gamma_list[[i]][[j]]
      save(WN_sync_gamma, file = paste0("../data/",
                                        "simulated_discrete_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_sync_gamma.RData"))

      BM_phyc_gamma_df <- as.data.frame(matrix(NA, nrow = length(BM_phyc_gamma$gamma_percent), ncol = num_cols))
      BM_sync_gamma_df <- as.data.frame(matrix(NA, nrow = length(BM_sync_gamma$gamma_percent), ncol = num_cols))
      WN_phyc_gamma_df <- as.data.frame(matrix(NA, nrow = length(WN_phyc_gamma$gamma_percent), ncol = num_cols))
      WN_sync_gamma_df <- as.data.frame(matrix(NA, nrow = length(WN_sync_gamma$gamma_percent), ncol = num_cols))

      colnames(BM_phyc_gamma_df) <- colnames(BM_sync_gamma_df) <-
        colnames(WN_phyc_gamma_df) <- colnames(WN_sync_gamma_df) <- summary_col_names

      BM_phyc_gamma_df$phenotype_phylogenetic_signal <- BM_sync_gamma_df$phenotype_phylogenetic_signal <- "BM"
      WN_phyc_gamma_df$phenotype_phylogenetic_signal <- WN_sync_gamma_df$phenotype_phylogenetic_signal <- "WN"

      BM_phyc_gamma_df$test <- WN_phyc_gamma_df$test <- "phyc"
      BM_sync_gamma_df$test <- WN_sync_gamma_df$test <- "sync"

      BM_phyc_gamma_df$tree_id <- i
      BM_phyc_gamma_df$phenotype_id <- j
      BM_phyc_gamma_df$gamma_avg <- rep(BM_phyc_gamma$gamma_avg, nrow(BM_phyc_gamma_df))
      BM_phyc_gamma_df$gamma_percent <- BM_phyc_gamma$gamma_percent
      BM_phyc_gamma_df$gamma_count <- BM_phyc_gamma$gamma_count
      BM_phyc_gamma_df$num_hi_conf_edges <- BM_phyc_gamma$num_hi_conf_edges
      BM_phyc_gamma_df$pheno_beta <- BM_phyc_gamma$pheno_beta
      BM_phyc_gamma_df$geno_beta <- BM_phyc_gamma$geno_beta
      BM_phyc_gamma_df$epsilon <- BM_phyc_gamma$epsilon
      BM_phyc_gamma_df$genotype <- BM_phyc_gamma$genotype
      
      BM_sync_gamma_df$tree_id <- i
      BM_sync_gamma_df$phenotype_id <- j
      BM_sync_gamma_df$gamma_avg <- rep(BM_sync_gamma$gamma_avg, nrow(BM_sync_gamma_df))
      BM_sync_gamma_df$gamma_percent <- BM_sync_gamma$gamma_percent
      BM_sync_gamma_df$gamma_count <- BM_sync_gamma$gamma_count
      BM_sync_gamma_df$num_hi_conf_edges <- BM_sync_gamma$num_hi_conf_edges
      BM_sync_gamma_df$pheno_beta <- BM_sync_gamma$pheno_beta
      BM_sync_gamma_df$geno_beta <- BM_sync_gamma$geno_beta
      BM_sync_gamma_df$epsilon <- BM_sync_gamma$epsilon
      BM_sync_gamma_df$genotype <- BM_sync_gamma_df$genotype
      
      WN_phyc_gamma_df$tree_id <- i
      WN_phyc_gamma_df$phenotype_id <- j
      WN_phyc_gamma_df$gamma_avg <- rep(WN_phyc_gamma$gamma_avg, nrow(WN_phyc_gamma_df))
      WN_phyc_gamma_df$gamma_percent <- WN_phyc_gamma$gamma_percent
      WN_phyc_gamma_df$gamma_count <- WN_phyc_gamma$gamma_count
      WN_phyc_gamma_df$num_hi_conf_edges <- WN_phyc_gamma$num_hi_conf_edges
      WN_phyc_gamma_df$pheno_beta <- WN_phyc_gamma$pheno_beta
      WN_phyc_gamma_df$geno_beta <- WN_phyc_gamma$geno_beta
      WN_phyc_gamma_df$epsilon <- WN_phyc_gamma$epsilon
      WN_phyc_gamma_df$genotype <- WN_phyc_gamma$genotype
      
      WN_sync_gamma_df$tree_id <- i
      WN_sync_gamma_df$phenotype_id <- j
      WN_sync_gamma_df$gamma_avg <- rep(WN_sync_gamma$gamma_avg, nrow(WN_sync_gamma_df))
      WN_sync_gamma_df$gamma_percent <- WN_sync_gamma$gamma_percent
      WN_sync_gamma_df$gamma_count <- WN_sync_gamma$gamma_count
      WN_sync_gamma_df$num_hi_conf_edges <- WN_sync_gamma$num_hi_conf_edges
      WN_sync_gamma_df$pheno_beta <- WN_sync_gamma$pheno_beta
      WN_sync_gamma_df$geno_beta <- WN_sync_gamma$geno_beta
      WN_sync_gamma_df$epsilon <- WN_sync_gamma$epsilon
      WN_sync_gamma_df$genotype <- WN_sync_gamma$genotype
      
      summary_df <- rbind(summary_df, BM_phyc_gamma_df, BM_sync_gamma_df, WN_phyc_gamma_df, WN_sync_gamma_df)
    }
  }

  summary_df <- summary_df[2:nrow(summary_df), , drop = FALSE] # remove the NA row
  save(summary_df, file = "../data/simulated_gamma_summary.RData")
  write.table(summary_df,
              sep = "\t",
              file = "../data/simulated_gamma_summary.tsv", 
              row.names = FALSE)
}
