# keep_good_genotypes <- function(BM_phyc_gamma_list, 
#                                 BM_sync_gamma_list,
#                                 WN_phyc_gamma_list,
#                                 WN_sync_gamma_list){
#   num_trees <- length(BM_phyc_gamma_list)
#   for (i in 1:num_trees) {
#     num_pheno <- length(BM_phyc_gamma_list[[i]])
#     for (j in 1:num_pheno) {
#       num_geno <- length(BM_phyc_gamma_list[[i]][[j]]$geno_beta)
#       unique_geno_beta <- unique(BM_phyc_gamma_list[[i]][[j]]$geno_beta)
#       unique_pheno_beta <- unique(BM_phyc_gamma_list[[i]][[j]]$pheno_beta)
#       
#     }
#   }
# }