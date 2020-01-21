save_geno_key_cont <- function(dir, prefix, num_tree, threshold, num_pheno){
  for (j in 1:num_pheno) { 
    for (i in 1:num_tree) {
      cont_fname <- paste0(dir, "hogwash_continuous_", prefix, i, "_pheno_", j, ".rda") 

      if (file.exists(cont_fname)){
        hog_cont_results <- local(get(load(cont_fname)))  
        cont_df <- as.data.frame(hog_cont_results$gamma$gamma_count /
                                   hog_cont_results$gamma$pheno_beta)
        colnames(cont_df) <- "gamma_over_beta"
        cont_df <- cont_df %>% 
          mutate("association_with_phenotype" = gamma_over_beta)
        cont_df$association_with_phenotype[cont_df$association_with_phenotype < threshold] <- "null"
        cont_df$association_with_phenotype[cont_df$association_with_phenotype != "null"] <- "associated"
        cont_df <- cont_df %>% 
          mutate("genotype" = row.names(hog_cont_results$hit_pvals))
        
        
        write.table(cont_df,
                    sep = "\t", 
                    file = paste0(dir, 
                                  "key_continuous_geno_", 
                                  prefix, 
                                  i, 
                                  "_pheno_", 
                                  j,
                                  ".tsv"), 
                    row.names = FALSE) 
      } else {
        print("File does not exist:")
        print(cont_fname)
      }
    }
  }
}

save_geno_key_disc <- function(dir, prefix, num_tree, threshold, num_pheno){
  for (j in 1:num_pheno) {
    for (i in 1:num_tree) {
      phyc_fname <- paste0(dir, "hogwash_phyc_", prefix, i, "_pheno_", j, ".rda") 
      sync_fname <- paste0(dir, "hogwash_synchronous_", prefix, i, "_pheno_", j, ".rda") 
      
      if (file.exists(phyc_fname) & file.exists(sync_fname)){
        hog_phyc_results <- local(get(load(phyc_fname)))  
        hog_sync_results <- local(get(load(sync_fname)))  
        phyc_df <- as.data.frame(hog_phyc_results$gamma$gamma_count / 
                                   hog_phyc_results$gamma$pheno_beta)
        colnames(phyc_df) <- "gamma_over_beta"
        phyc_df <- phyc_df %>% 
          mutate("association_with_phenotype" = gamma_over_beta)
        phyc_df$association_with_phenotype[phyc_df$association_with_phenotype < threshold] <- "null"
        phyc_df$association_with_phenotype[phyc_df$association_with_phenotype != "null"] <- "associated"
        phyc_df <- phyc_df %>% 
          mutate("genotype" = names(hog_phyc_results$contingency_table))
        
        sync_df <- as.data.frame(hog_sync_results$gamma$gamma_count / 
                                   hog_sync_results$gamma$pheno_beta)
        colnames(sync_df) <- "gamma_over_beta"
        sync_df <- sync_df %>% 
          mutate("association_with_phenotype" = gamma_over_beta)
        sync_df$association_with_phenotype[sync_df$association_with_phenotype < threshold] <- "null"
        sync_df$association_with_phenotype[sync_df$association_with_phenotype != "null"] <- "associated"
        sync_df <- sync_df %>% 
          mutate("genotype" = names(hog_sync_results$contingency_table))
        
        write.table(phyc_df,
                    sep = "\t", 
                    file = paste0(dir, 
                                  "key_phyc_geno_", 
                                  prefix, 
                                  i, 
                                  "_pheno_", 
                                  j, 
                                  ".tsv"), 
                    row.names = FALSE) 
        write.table(sync_df,
                    sep = "\t", 
                    file = paste0(dir, 
                                  "key_synchronous_geno_", 
                                  prefix, 
                                  i, 
                                  "_pheno_",
                                  j, 
                                  ".tsv"), 
                    row.names = FALSE) 
      } else {
        print("Files do not exist:")
        print(phyc_fname)
        print(sync_fname)
      }
    }
  }
}
