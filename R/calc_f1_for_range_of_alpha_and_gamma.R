# TODO Change this script into a function
# run from /lib/

library(tidyverse)
data_dir <- "../data/"
simulated_data_dir <- "../data/" # TODO change to user defined input
phenotype_types <- "discrete" #c("discrete", "continuous")
data_types <- c("BM", "WN")
tree_index <- c(1:1)
num_tree <- length(tree_index)
pheno_index <- c(1:3)
num_pheno <- length(pheno_index)
phenotype <- "pheno"
genotype <- "geno"
discrete_test_types <- c("phyc", "synchronous")
# continuous_test_types <- "continuous"
alphas <- c(-log(0.05), -log(0.01), -log(0.005), -log(0.001), -log(0.0005), -log(0.0001))
epsilons <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)

discrete_num_rows <- 
  length(pheno_index) * 
  length(discrete_test_types) * 
  length(data_types) * 
  length(tree_index) * 
  length(alphas) * 
  length(epsilons)

column_names <- c("phenotype_type", 
                  "phenotype_phylogenetic_signal",
                  "tree_id", 
                  "test",
                  "phenotype_id", 
                  "num_tip", 
                  "num_geno", 
                  "alpha_threshold",
                  "epsilon_threshold",
                  "observed_gamma_value", 
                  "observed_delta_epsilon", 
                  "observed_pheno_beta", 
                  "observed_mean_geno_beta", 
                  "phenotype_state_balance",
                  "true_positive", 
                  "true_negative", 
                  "false_positive",
                  "false_negative", 
                  "recall",  # aka sensitivity
                  "precision", 
                  "F1_score", 
                  "false_positive_rate", 
                  "positive_predictive_value")

num_col <- length(column_names)
discrete_f1_mat <- matrix(NA, nrow = discrete_num_rows, ncol = num_col)
colnames(discrete_f1_mat) <- column_names

discrete_f1_tb <- as_tibble(discrete_f1_mat)
discrete_f1_tb$phenotype_type <- "discrete"

current_disc_row <- 0 

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
        
        if (file.exists(current_file_name)) {
          hogwash_results <- local(get(load(current_file_name)))
          pvals <- rownames_to_column(hogwash_results$hit_pvals, 
                                      var = "genotype") 
          pvals_epsilon_tb <- cbind(pvals, hogwash_results$gamma$epsilon)
          colnames(pvals_epsilon_tb) <- 
            c("genotype", "fdr_corrected_pvals", "epsilon")
          
          gamma_value <- hogwash_results$gamma$gamma_avg
          beta_pheno <- hogwash_results$gamma$pheno_beta
          mean_beta_geno <- mean(hogwash_results$gamma$geno_beta)
          pheno_presence <- pheno_absence <- NULL 
          for (p in 1:length(hogwash_results$contingency_table)) {
            pheno_presence <- c(pheno_presence, sum(hogwash_results$contingency_table[[p]][, 1]))
            pheno_absence <- c(pheno_absence, sum(hogwash_results$contingency_table[[p]][, 2]))
          }
          avg_pheno_state_balance <- mean(pheno_presence) / sum(mean(pheno_presence), mean(pheno_absence))
          
          hogwash_results <- NULL
          
          
          for (p in 1:length(alphas)) {
            for (q in 1:length(epsilons)) {
              true_positive <- pvals_epsilon_tb %>% filter(epsilon >= epsilons[q] & fdr_corrected_pvals >= alphas[p]) %>% nrow()
              true_negative <- pvals_epsilon_tb %>% filter(epsilon < epsilons[q] & fdr_corrected_pvals < alphas[p]) %>% nrow()
              false_positive <- pvals_epsilon_tb %>% filter(epsilon < epsilons[q] & fdr_corrected_pvals >= alphas[p]) %>% nrow()
              false_negative <- pvals_epsilon_tb %>% filter(epsilon >= epsilons[q] & fdr_corrected_pvals < alphas[p]) %>% nrow()
              
              precision <- true_positive / sum(true_positive, false_positive)
              recall <- true_positive / sum(false_negative, true_positive)
              false_positive_rate <- false_positive / sum(false_positive, true_negative)
              positive_predictive_value <- true_positive / sum(true_positive, false_positive)
              
              f1_score <- (2 * precision * recall) / (precision + recall)
              
              # calculate delta epsilon: 
              null_epsilon_dist <- pvals_epsilon_tb$epsilon[pvals_epsilon_tb$epsilon < epsilons[q]]
              assoc_epsilon_dist <- pvals_epsilon_tb$epsilon[pvals_epsilon_tb$epsilon >= epsilons[q]]
              delta_epsilon <- mean(assoc_epsilon_dist) - mean(null_epsilon_dist)
              
              current_disc_row <- current_disc_row + 1
              discrete_f1_tb$phenotype_phylogenetic_signal[current_disc_row] <- as.character(data_types[j])
              discrete_f1_tb$tree_id[current_disc_row]  <- tree_index[k]
              discrete_f1_tb$test[current_disc_row] <- discrete_test_types[l]
              discrete_f1_tb$phenotype_id[current_disc_row] <- pheno_index[h]
              discrete_f1_tb$F1_score[current_disc_row] <- f1_score
              discrete_f1_tb$true_positive[current_disc_row] <- true_positive
              discrete_f1_tb$true_negative[current_disc_row] <- true_negative
              discrete_f1_tb$false_negative[current_disc_row] <- false_negative
              discrete_f1_tb$false_positive[current_disc_row] <- false_positive
              discrete_f1_tb$num_geno[current_disc_row] <- nrow(pvals_epsilon_tb)
              discrete_f1_tb$observed_gamma_value[current_disc_row] <- gamma_value
              discrete_f1_tb$recall[current_disc_row] <- recall
              discrete_f1_tb$precision[current_disc_row] <- precision
              discrete_f1_tb$false_positive_rate[current_disc_row] <- false_positive_rate
              discrete_f1_tb$positive_predictive_value[current_disc_row] <- positive_predictive_value
              discrete_f1_tb$observed_delta_epsilon[current_disc_row] <- delta_epsilon
              discrete_f1_tb$observed_pheno_beta[current_disc_row] <- beta_pheno
              discrete_f1_tb$observed_mean_geno_beta[current_disc_row] <- mean_beta_geno
              discrete_f1_tb$phenotype_state_balance[current_disc_row] <- avg_pheno_state_balance
              discrete_f1_tb$alpha_threshold[current_disc_row] <- alphas[p]
              discrete_f1_tb$epsilon_threshold[current_disc_row] <- epsilons[q]
            }
          }
        } else {
          print("File does not exist")
          print(current_file_name)
        }
      }
    }
  }
}


discrete_f1_tb$test[discrete_f1_tb$test == "synchronous"] <- "sync"

write_tsv(discrete_f1_tb,
          path = "../data/F1_scores_range_of_alpha_gamma.tsv",
          col_names = TRUE)

