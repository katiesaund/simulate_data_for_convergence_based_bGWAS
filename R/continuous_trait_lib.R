make_continuous_phenotypes <- function(tree_list, num_pheno){
  num_trees <- length(tree_list)
  set.seed(1)
  for (j in 1:num_pheno) {
    for (i in 1:num_trees) {
      continuous_BM_pheno <- phytools::fastBM(tree = tree_list[[i]])
      
      # When lambda is close to zero, the phylogenetic signal is low
      lambda_has_high_signal <- TRUE
      while (lambda_has_high_signal) {
        jumbled_pheno <- sample(unname(continuous_BM_pheno), 
                                size = ape::Ntip(tree_list[[i]]), 
                                replace = FALSE)
        names(jumbled_pheno) <- tree_list[[i]]$tip.label
        jumbled_lambda <- phytools::phylosig(tree = tree_list[[i]],
                                             x = jumbled_pheno,
                                             method = "lambda")
        lambda_has_high_signal <- jumbled_lambda$lambda > 0.1
      }
      cont_low_signal_pheno <- jumbled_pheno
      
      continuous_BM_pheno <- 
        matrix(c(names(continuous_BM_pheno), continuous_BM_pheno), ncol = 2)  
      row.names(continuous_BM_pheno) <- continuous_BM_pheno[, 1]
      continuous_BM_pheno <- continuous_BM_pheno[, 2, drop = FALSE]
      
      cont_low_signal_pheno <- matrix(c(names(cont_low_signal_pheno), cont_low_signal_pheno), ncol = 2) 
      row.names(cont_low_signal_pheno) <- cont_low_signal_pheno[, 1]
      cont_low_signal_pheno <- cont_low_signal_pheno[, 2, drop = FALSE]
      colnames(continuous_BM_pheno) <- colnames(cont_low_signal_pheno) <- "pheno"
      
      write.table(continuous_BM_pheno, 
                  sep = "\t", 
                  file = paste0("../data/", 
                                "continuous_pheno_BM_tree_", 
                                i, 
                                "_pheno_", 
                                j, 
                                ".tsv"))
      write.table(cont_low_signal_pheno, 
                  sep = "\t", 
                  file = paste0("../data/", 
                                "continuous_pheno_WN_tree_", 
                                i, 
                                "_pheno_",
                                j,
                                ".tsv"))
    }
  }
}

                  