# TODO Change this script into a function
# TODO add generate_association_key functions here?
# run from /lib/
source("generate_association_key.R")
# F1 Score stuff
library(tidyverse)
assoc_threshold <- 0.50 # TODO Change to user defined input

data_dir <- "../data/"
simulated_data_dir <- "../data/" # TODO change to user defined input
phenotype_types <- "discrete" #c("discrete", "continuous")
data_types <- c("BM", "WN")
tree_index <- c(1:3)
num_tree <- length(tree_index)
pheno_index <- c(1:2)
num_pheno <- length(pheno_index)
phenotype <- "pheno"
genotype <- "geno"
discrete_test_types <- c("phyc", "synchronous")
continuous_test_types <- "continuous"
alpha <- -log(0.10)

# Generate association key
# save_geno_key_cont(data_dir, "continuous_pheno_BM_tree_", num_tree, assoc_threshold, num_pheno)
# save_geno_key_cont(data_dir, "continuous_pheno_WN_tree_", num_tree, assoc_threshold, num_pheno)
save_geno_key_disc(data_dir, "discrete_pheno_BM_tree_", num_tree, assoc_threshold, num_pheno)
save_geno_key_disc(data_dir, "discrete_pheno_WN_tree_", num_tree, assoc_threshold, num_pheno)

# Back to F1
discrete_num_rows <- length(pheno_index) * length(discrete_test_types) * length(data_types) * length(tree_index)
# continuous_num_rows <- length(pheno_index) * length(data_types) * length(tree_index)
num_col <- 21

discrete_f1_mat <- matrix(NA, nrow = discrete_num_rows, ncol = num_col)
# continuous_f1_mat <- matrix(NA, nrow = continuous_num_rows, ncol = num_col)

# colnames(continuous_f1_mat) <-
  colnames(discrete_f1_mat) <- 
  c("phenotype_type", 
    "phenotype_phylogenetic_signal",
    "tree_id", 
    "test",
    "F1_score", 
    "true_positive", 
    "true_negative", 
    "false_positive",
    "false_negative", 
    "num_tip", 
    "num_geno", 
    "observed_gamma_value", 
    "recall",  # aka sensitivity
    "precision", 
    "false_positive_rate", 
    "positive_predictive_value", 
    "phenotype_id", 
    "observed_delta_epsilon", 
    "observed_pheno_beta", 
    "observed_mean_geno_beta", 
    "phenotype_state_balance")

discrete_f1_mat[, 1] <- "discrete"
# continuous_f1_mat[, 1] <- "continuous"
current_disc_row <- current_cont_row <- 0 

# DISCRETE
for (h in 1:length(pheno_index)) {
    for (j in 1:length(data_types)) {
      for (k in 1:length(tree_index)) {
          for (l in 1:length(discrete_test_types)) {
            current_file_name <- paste0(data_dir, 
                                        "hogwash_",
                                        discrete_test_types[l],
                                        "_discrete_", 
                                        phenotype,
                                        "_", 
                                        data_types[j],
                                        "_tree_", 
                                        tree_index[k],
                                        "_pheno_", 
                                        pheno_index[h],
                                        ".rda")
            current_key_name <- paste0(data_dir,
                                       "key_", 
                                       discrete_test_types[l],
                                       "_",
                                       genotype,
                                       "_discrete_",
                                       phenotype, 
                                       "_",
                                       data_types[j], 
                                       "_tree_", 
                                       tree_index[k],
                                       "_pheno_", 
                                       pheno_index[h],
                                       ".tsv")
            
            if (file.exists(current_file_name) & file.exists(current_key_name)) {
              hogwash_results <- local(get(load(current_file_name)))
              pvals <- rownames_to_column(hogwash_results$hit_pvals,
                                          var = "genotype") 
              key <- read.table(current_key_name, 
                                header = 1, 
                                stringsAsFactors = FALSE)
              
              pvals_and_key <- left_join(pvals, key, by = "genotype")
              
              true_positive <- pvals_and_key %>% filter(association_with_phenotype == "associated" & fdr_corrected_pvals > alpha) %>% nrow()
              true_negative <- pvals_and_key %>% filter(association_with_phenotype == "null" & fdr_corrected_pvals < alpha) %>% nrow()
              false_positive <- pvals_and_key %>% filter(association_with_phenotype == "null" & fdr_corrected_pvals > alpha) %>% nrow()
              false_negative <- pvals_and_key %>% filter(association_with_phenotype == "associated" & fdr_corrected_pvals < alpha) %>% nrow()
              
              precision <- true_positive / sum(true_positive, false_positive)
              recall <- true_positive / sum(false_negative, true_positive)
              false_positive_rate <- false_positive / sum(false_positive, true_negative)
              positive_predictive_value <- true_positive / sum(true_positive, false_positive)
              
              f1_score <- (2 * precision * recall) / (precision + recall)
              
              gamma_value <- hogwash_results$gamma$gamma_avg
              
              # calculate delta epsilon: 
              null_epsilon_dist <- hogwash_results$gamma$epsilon[key$association_with_phenotype == "null"]
              assoc_epsilon_dist <- hogwash_results$gamma$epsilon[key$association_with_phenotype == "associated"]
              delta_epsilon <- mean(assoc_epsilon_dist) - mean(null_epsilon_dist)
              beta_pheno <- hogwash_results$gamma$pheno_beta
              mean_beta_geno <- mean(hogwash_results$gamma$geno_beta)
              
              pheno_presence <- pheno_absence <- NULL 
              for (p in 1:length(hogwash_results$contingency_table)) {
                pheno_presence <- c(pheno_presence, sum(hogwash_results$contingency_table[[p]][, 1]))
                pheno_absence <- c(pheno_absence, sum(hogwash_results$contingency_table[[p]][, 2]))
              }
              avg_pheno_state_balance <- mean(pheno_presence) / sum(mean(pheno_presence), mean(pheno_absence))
              
              hogwash_results <- NULL

              current_disc_row <- current_disc_row + 1
              discrete_f1_mat[current_disc_row, 2] <- as.character(data_types[j])
              discrete_f1_mat[current_disc_row, 3] <- tree_index[k]
              discrete_f1_mat[current_disc_row, 4] <- discrete_test_types[l]
              discrete_f1_mat[current_disc_row, 5] <- f1_score
              discrete_f1_mat[current_disc_row, 6] <- true_positive
              discrete_f1_mat[current_disc_row, 7] <- true_negative
              discrete_f1_mat[current_disc_row, 8] <- false_positive
              discrete_f1_mat[current_disc_row, 9] <- false_negative
              discrete_f1_mat[current_disc_row, 10] <- current_num_tip <- NA
              discrete_f1_mat[current_disc_row, 11] <- current_num_geno <- NA
              discrete_f1_mat[current_disc_row, 12] <- gamma_value
              discrete_f1_mat[current_disc_row, 13] <- recall
              discrete_f1_mat[current_disc_row, 14] <- precision
              discrete_f1_mat[current_disc_row, 15] <- false_positive_rate
              discrete_f1_mat[current_disc_row, 16] <- positive_predictive_value
              discrete_f1_mat[current_disc_row, 17] <- pheno_index[h]
              discrete_f1_mat[current_disc_row, 18] <- delta_epsilon
              discrete_f1_mat[current_disc_row, 19] <- beta_pheno
              discrete_f1_mat[current_disc_row, 20] <- mean_beta_geno
              discrete_f1_mat[current_disc_row, 21] <- avg_pheno_state_balance
            } else {
              print("File does not exist")
              print(current_file_name)
            }
      }
    }
  }
}

# CONTINUOUS
# for (h in 1:length(pheno_index)) {
  for (j in 1:length(data_types)) {
    for (k in 1:length(tree_index)) {
      current_file_name <- paste0(data_dir, 
                                  "hogwash_continuous_continuous_pheno_",
                                  data_types[j],
                                  "_tree_", 
                                  tree_index[k],
                                  "_pheno_", 
                                  pheno_index[h],
                                  ".rda")
      current_key_name <- paste0(data_dir,
                                 "key_continuous_geno_continuous_pheno_",
                                 data_types[j], 
                                 "_tree_", 
                                 tree_index[k],
                                 "_pheno_", 
                                 pheno_index[h],
                                 ".tsv")
      
      if (file.exists(current_file_name) & file.exists(current_key_name)) {
        hogwash_results <- local(get(load(current_file_name)))
        pvals <- rownames_to_column(hogwash_results$hit_pvals, var = "genotype") 
        key <- read.table(current_key_name, header = 1, stringsAsFactors = FALSE)
        
        pvals_and_key <- left_join(pvals, key, by = "genotype")
        
        true_positive <- pvals_and_key %>% filter(association_with_phenotype == "associated" & fdr_corrected_pvals > alpha) %>% nrow()
        true_negative <- pvals_and_key %>% filter(association_with_phenotype == "null" & fdr_corrected_pvals < alpha) %>% nrow()
        false_positive <- pvals_and_key %>% filter(association_with_phenotype == "null" & fdr_corrected_pvals > alpha) %>% nrow()
        false_negative <- pvals_and_key %>% filter(association_with_phenotype == "associated" & fdr_corrected_pvals < alpha) %>% nrow()
        
        precision <- true_positive / sum(true_positive, false_positive)
        recall <- true_positive / sum(false_negative, true_positive)
        false_positive_rate <- false_positive / sum(false_positive, true_negative)
        positive_predictive_value <- true_positive / sum(true_positive, false_positive)
        
        f1_score <- (2 * precision * recall) / (precision + recall)
        
        hogwash_results <- local(get(load(current_file_name)))
        gamma_value <- hogwash_results$gamma$gamma_avg
        # calculate delta epsilon: 
        null_epsilon_dist <- hogwash_results$gamma$epsilon[key$association_with_phenotype == "null"]
        assoc_epsilon_dist <- hogwash_results$gamma$epsilon[key$association_with_phenotype == "associated"]
        delta_epsilon <- mean(assoc_epsilon_dist) - mean(null_epsilon_dist)
        beta_pheno <- hogwash_results$gamma$pheno_beta
        mean_beta_geno <- mean(hogwash_results$gamma$geno_beta)
        
        pheno_presence <- pheno_absence <- NULL 
        for (p in 1:length(hogwash_results$contingency_table)) {
          pheno_presence <- c(pheno_presence, sum(hogwash_results$contingency_table[[p]][, 1]))
          pheno_absence <- c(pheno_absence, sum(hogwash_results$contingency_table[[p]][, 2]))
        }
        avg_pheno_state_balance <- mean(pheno_presence) / sum(mean(pheno_presence), mean(pheno_absence))
        
        
        hogwash_results <- NULL
        
        current_cont_row <- current_cont_row + 1
        continuous_f1_mat[current_cont_row, 2] <- as.character(data_types[j])
        continuous_f1_mat[current_cont_row, 3] <- tree_index[k]
        continuous_f1_mat[current_cont_row, 4] <- "continuous"
        continuous_f1_mat[current_cont_row, 5] <- f1_score
        continuous_f1_mat[current_cont_row, 6] <- true_positive
        continuous_f1_mat[current_cont_row, 7] <- true_negative
        continuous_f1_mat[current_cont_row, 8] <- false_positive
        continuous_f1_mat[current_cont_row, 9] <- false_negative
        continuous_f1_mat[current_cont_row, 10] <- current_num_tip <- NA
        continuous_f1_mat[current_cont_row, 11] <- current_num_geno <- NA
        continuous_f1_mat[current_cont_row, 12] <- gamma_value
        continuous_f1_mat[current_cont_row, 13] <- recall
        continuous_f1_mat[current_cont_row, 14] <- precision
        continuous_f1_mat[current_cont_row, 15] <- false_positive_rate
        continuous_f1_mat[current_cont_row, 16] <- positive_predictive_value
        continuous_f1_mat[current_cont_row, 17] <- pheno_index[h]
        continuous_f1_mat[current_cont_row, 18] <- delta_epsilon
        continuous_f1_mat[current_cont_row, 19] <- beta_pheno
        continuous_f1_mat[current_cont_row, 20] <- mean_beta_geno
        continuous_f1_mat[current_cont_row, 21] <- avg_pheno_state_balance
        
      } else {
        print("file does not exist:")
        print(current_file_name)
      }
        
    }
  }
}
      
# f1_mat <- rbind(continuous_f1_mat, discrete_f1_mat)
f1_mat <- discrete_f1_mat
f1_mat[f1_mat[, 4] == "synchronous", 4] <- "sync"

write.table(f1_mat, 
            file = "../data/F1_scores.tsv", 
            sep = "\t", 
            col.names = TRUE,
            row.names = FALSE)
