#' Keep a subset of all of the genotypes that have been made
#' 
#' @description Rather than keep all of the simulated data, which might be
#'   highly similar, we reduce the amount of redundant data. We create a "bin"
#'   size which is the max number of genotypes to keep that have any particular
#'   beta(genotype) & beta(phenotype) value pair. The hope is to reduce the size
#'   of the dataset while maintaining the diversity of the dataset.
#'
#' @param convergence_list Convergence values
#' @param bin_size Integer
#'
#' @return
#' @noRd
keep_good_genotypes <- function(convergence_list, bin_size){

  # Initialize list of genotypes to keep 
  keepers <- list()
  num_trees <- length(convergence_list)
  
  for (i in 1:num_trees) {
    num_pheno <- length(convergence_list[[i]])
    temp_keepers <- list()
    
    for (j in 1:num_pheno) {
      current_index <- NULL
      num_geno <- length(convergence_list[[i]][[j]]$geno_beta)
      unique_geno_beta <- unique(convergence_list[[i]][[j]]$geno_beta)
      unique_pheno_beta <- unique(convergence_list[[i]][[j]]$pheno_beta)
      
      for (k in 1:length(unique_geno_beta)) {
        
        for (l in 1:length(unique_pheno_beta)) {
          
          # Find the indices of the genotypes inside this beta value bin
          temp_index <- 
            which(convergence_list[[i]][[j]]$geno_beta == unique_geno_beta[k] &
                    convergence_list[[i]][[j]]$pheno_beta == unique_pheno_beta[l])
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
