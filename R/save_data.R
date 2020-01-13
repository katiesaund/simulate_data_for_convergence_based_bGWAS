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
                      WN_sync_gamma_list) {
  
  num_trees <- length(tree_list) 
  for (i in 1:num_trees) {
    current_tree <- tree_list[[i]]
    num_tips <- ape::Ntip(current_tree)
    num_pheno <- ncol(BM_phenotype_AR_mat_list[[i]])
    
    write.table(genotype_AR_mat_list[[i]], 
                sep = "\t", 
                row.names = TRUE, 
                file = paste0("../data/",
                              "simulated_genotype_",
                              i, 
                              ".tsv"))
    
    write.tree(current_tree, 
               file = paste0("../data/", 
                             "simulated_tree_",
                             i, 
                             ".tree"))
    
    for (j in 1:num_pheno) {
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
      write.table(temp_disc_BM_trait, 
                  sep = "\t", 
                  file = paste0("../data/", 
                                "discrete_pheno_WN_tree_", 
                                i, 
                                "_pheno_", 
                                j, 
                                ".tsv"))
    }
  }
}