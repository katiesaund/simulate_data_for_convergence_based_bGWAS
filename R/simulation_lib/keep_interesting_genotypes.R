keep_good_genotypes <- function(gamma_list, bin_size){
  keepers <- list()
  num_trees <- length(gamma_list)
  for (i in 1:num_trees) {
    num_pheno <- length(gamma_list[[i]])
    temp_keepers <- list()
    for (j in 1:num_pheno) {
      current_index <- NULL
      num_geno <- length(gamma_list[[i]][[j]]$geno_beta)
      unique_geno_beta <- unique(gamma_list[[i]][[j]]$geno_beta)
      unique_pheno_beta <- unique(gamma_list[[i]][[j]]$pheno_beta)
      for (k in 1:length(unique_geno_beta)) {
        for (l in 1:length(unique_pheno_beta)) {
          temp_index <- which(gamma_list[[i]][[j]]$geno_beta == unique_geno_beta[k] &
                                gamma_list[[i]][[j]]$pheno_beta == unique_pheno_beta[l])
          if (length(temp_index) > bin_size) {
            temp_index <- temp_index[1:bin_size]
          }
          current_index <- c(current_index, temp_index)
          current_index <- unique(sort(current_index))
        }
      }
    temp_keepers[[j]] <- current_index
    }
  keepers[[i]] <- temp_keepers
  }
  return(keepers)
}
