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

      BM_geno_mat <- genotype_AR_mat_list[[i]][, BM_keepers, drop = FALSE]
      WN_geno_mat <- genotype_AR_mat_list[[i]][, WN_keepers, drop = FALSE]

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
                                "discrete_pheno_BM_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))

      temp_disc_WN_trait <-
        WN_phenotype_AR_mat_list[[i]][1:num_tips, j, drop = FALSE]
      write.table(temp_disc_WN_trait,
                  sep = "\t",
                  file = paste0("../data/",
                                "discrete_pheno_WN_tree_",
                                i,
                                "_pheno_",
                                j,
                                ".tsv"))

      # Save gamma results ----
      BM_phyc_gamma <- BM_phyc_gamma_list[[i]][[j]]
      save(BM_phyc_gamma, file = paste0("../data/",
                                        "discrete_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_phyc_gamma.RData"))
      BM_sync_gamma <- BM_sync_gamma_list[[i]][[j]]
      save(BM_sync_gamma, file = paste0("../data/",
                                        "discrete_pheno_BM_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_sync_gamma.RData"))

      WN_phyc_gamma <- WN_phyc_gamma_list[[i]][[j]]
      save(WN_phyc_gamma, file = paste0("../data/",
                                        "discrete_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_phyc_gamma.RData"))
      WN_sync_gamma <- WN_sync_gamma_list[[i]][[j]]
      save(WN_sync_gamma, file = paste0("../data/",
                                        "discrete_pheno_WN_tree_",
                                        i,
                                        "_pheno_",
                                        j,
                                        "_sync_gamma.RData"))
    }
  }
}
