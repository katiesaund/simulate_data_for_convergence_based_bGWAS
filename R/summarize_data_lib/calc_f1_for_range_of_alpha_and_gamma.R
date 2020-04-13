# TODO Change this script into a function
# run from /lib/
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- "../data/"
phenotype_types <- c("discrete", "continuous")
data_types <- c("BM", "WN")
tree_index <- c(1:as.numeric(args[1]))
num_tree <- length(tree_index)
pheno_index <- c(1:as.numeric(args[2]))
num_pheno <- length(pheno_index)
phenotype <- "pheno"
genotype <- "geno"
discrete_test_types <- c("phyc", "synchronous")
continuous_test_types <- "continuous"
alphas <- c(-log(0.05), -log(0.01), -log(0.005), -log(0.001), -log(0.0005), -log(0.0001), -log(0.00005), -log(0.00001))
epsilons <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)

discrete_num_rows <- 
  length(pheno_index) * 
  length(discrete_test_types) * 
  length(data_types) * 
  length(tree_index) * 
  length(alphas) * 
  length(epsilons)

continuous_num_rows <- 
  length(pheno_index) * 
  length(continuous_test_types) * 
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
                  "observed_mean_pheno_beta", 
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
continuous_f1_mat <- matrix(NA, nrow = continuous_num_rows, ncol = num_col)
colnames(discrete_f1_mat) <- colnames(continuous_f1_mat) <- column_names

discrete_f1_tb <- as_tibble(discrete_f1_mat)
continuous_f1_tb <- as_tibble(continuous_f1_mat)

col_names_2 <- c(column_names,"genotype", "fdr_corrected_pvals", "epsilon")
discrete_f1_tb_and_geno <- matrix(NA, nrow = 1, ncol = length(col_names_2))
colnames(discrete_f1_tb_and_geno) <- col_names_2
discrete_f1_tb_and_geno <- as.data.frame(discrete_f1_tb_and_geno)

continuous_f1_tb_and_geno <- matrix(NA, nrow = 1, ncol = length(col_names_2))
colnames(continuous_f1_tb_and_geno) <- col_names_2
continuous_f1_tb_and_geno <- as.data.frame(continuous_f1_tb_and_geno)

current_cont_row <- current_disc_row <- 0 
# DISCRETE
for (h in 1:length(pheno_index)) {
  for (j in 1:length(data_types)) {
    for (k in 1:length(tree_index)) {
      for (l in 1:length(discrete_test_types)) {
        current_file_name <- paste0(data_dir, 
                                    "hogwash_",
                                    discrete_test_types[l],
                                    "_",
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
          mean_beta_pheno <- mean(hogwash_results$gamma$pheno_beta)
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
              
              precision <- positive_predictive_value <- true_positive / sum(true_positive, false_positive)
              recall <- true_positive / sum(false_negative, true_positive)
              false_positive_rate <- false_positive / sum(false_positive, true_negative)

              f1_score <- (2 * precision * recall) / (precision + recall)
              
              # calculate delta epsilon: 
              null_epsilon_dist <- pvals_epsilon_tb$epsilon[pvals_epsilon_tb$epsilon < epsilons[q]]
              assoc_epsilon_dist <- pvals_epsilon_tb$epsilon[pvals_epsilon_tb$epsilon >= epsilons[q]]
              delta_epsilon <- mean(assoc_epsilon_dist) - mean(null_epsilon_dist)
              
              # TODO: add genotype data
              temp_tb <- as_tibble(matrix(NA, nrow = 1, ncol = num_col))
              colnames(temp_tb) <- column_names
              big_tb <- cbind(pvals_epsilon_tb, temp_tb)
              
              
              current_disc_row <- current_disc_row + 1
              discrete_f1_tb$phenotype_type[current_disc_row] <- big_tb$phenotype_type <- "discrete"
              discrete_f1_tb$phenotype_phylogenetic_signal[current_disc_row] <- big_tb$phenotype_phylogenetic_signal <- as.character(data_types[j])
              discrete_f1_tb$tree_id[current_disc_row] <- big_tb$tree_id <- tree_index[k]
              discrete_f1_tb$test[current_disc_row] <- big_tb$test <- discrete_test_types[l]
              discrete_f1_tb$phenotype_id[current_disc_row] <- big_tb$phenotype_id <- pheno_index[h]
              discrete_f1_tb$F1_score[current_disc_row] <- big_tb$F1_score <- f1_score
              discrete_f1_tb$true_positive[current_disc_row] <- big_tb$true_positive <- true_positive
              discrete_f1_tb$true_negative[current_disc_row] <- big_tb$true_negative <- true_negative
              discrete_f1_tb$false_negative[current_disc_row] <- big_tb$false_negative <- false_negative
              discrete_f1_tb$false_positive[current_disc_row] <- big_tb$false_positive <- false_positive
              discrete_f1_tb$num_geno[current_disc_row] <- big_tb$num_geno <- nrow(pvals_epsilon_tb)
              discrete_f1_tb$observed_gamma_value[current_disc_row] <- big_tb$observed_gamma_value <- gamma_value
              discrete_f1_tb$recall[current_disc_row] <- big_tb$recall <- recall
              discrete_f1_tb$precision[current_disc_row] <- big_tb$precision <- precision
              discrete_f1_tb$false_positive_rate[current_disc_row] <- big_tb$false_positive_rate <- false_positive_rate
              discrete_f1_tb$positive_predictive_value[current_disc_row] <- big_tb$positive_predictive_value <- positive_predictive_value
              discrete_f1_tb$observed_delta_epsilon[current_disc_row] <- big_tb$observed_delta_epsilon <- delta_epsilon
              discrete_f1_tb$observed_mean_pheno_beta[current_disc_row] <- big_tb$observed_mean_pheno_beta <- mean_beta_pheno
              discrete_f1_tb$observed_mean_geno_beta[current_disc_row] <- big_tb$observed_mean_geno_beta <- mean_beta_geno
              discrete_f1_tb$phenotype_state_balance[current_disc_row] <- big_tb$phenotype_state_balance <- avg_pheno_state_balance
              discrete_f1_tb$alpha_threshold[current_disc_row] <- big_tb$alpha_threshold <- alphas[p]
              discrete_f1_tb$epsilon_threshold[current_disc_row] <- big_tb$epsilon_threshold <- epsilons[q]
              
              discrete_f1_tb_and_geno <- rbind(discrete_f1_tb_and_geno, big_tb)
              
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


# CONTINUOUS
for (h in 1:length(pheno_index)) {
  for (j in 1:length(data_types)) {
    for (k in 1:length(tree_index)) {
      for (l in 1:length(continuous_test_types)) {
        current_file_name <- paste0(data_dir, 
                                    "hogwash_continuous_continuous_", 
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
          mean_beta_pheno <- mean(hogwash_results$gamma$pheno_beta)
          mean_beta_geno <- mean(hogwash_results$gamma$geno_beta)
          hogwash_results <- NULL
          for (p in 1:length(alphas)) {
            for (q in 1:length(epsilons)) {
              true_positive <- pvals_epsilon_tb %>% filter(epsilon >= epsilons[q] & fdr_corrected_pvals >= alphas[p]) %>% nrow()
              true_negative <- pvals_epsilon_tb %>% filter(epsilon < epsilons[q] & fdr_corrected_pvals < alphas[p]) %>% nrow()
              false_positive <- pvals_epsilon_tb %>% filter(epsilon < epsilons[q] & fdr_corrected_pvals >= alphas[p]) %>% nrow()
              false_negative <- pvals_epsilon_tb %>% filter(epsilon >= epsilons[q] & fdr_corrected_pvals < alphas[p]) %>% nrow()
              
              precision <- positive_predictive_value <- true_positive / sum(true_positive, false_positive)
              recall <- true_positive / sum(false_negative, true_positive)
              false_positive_rate <- false_positive / sum(false_positive, true_negative)
              
              f1_score <- (2 * precision * recall) / (precision + recall)
              
              # calculate delta epsilon: 
              null_epsilon_dist <- pvals_epsilon_tb$epsilon[pvals_epsilon_tb$epsilon < epsilons[q]]
              assoc_epsilon_dist <- pvals_epsilon_tb$epsilon[pvals_epsilon_tb$epsilon >= epsilons[q]]
              delta_epsilon <- mean(assoc_epsilon_dist) - mean(null_epsilon_dist)
              
              # TODO: add genotype data
              temp_tb <- as_tibble(matrix(NA, nrow = 1, ncol = num_col))
              colnames(temp_tb) <- column_names
              big_tb <- cbind(pvals_epsilon_tb, temp_tb)
              
              
              current_cont_row <- current_cont_row + 1
              continuous_f1_tb$phenotype_type[current_cont_row] <- big_tb$phenotype_type <- "continuous"
              continuous_f1_tb$phenotype_phylogenetic_signal[current_cont_row] <- big_tb$phenotype_phylogenetic_signal <- as.character(data_types[j])
              continuous_f1_tb$tree_id[current_cont_row] <- big_tb$tree_id <- tree_index[k]
              continuous_f1_tb$test[current_cont_row] <- big_tb$test <- continuous_test_types[l]
              continuous_f1_tb$phenotype_id[current_cont_row] <- big_tb$phenotype_id <- pheno_index[h]
              continuous_f1_tb$F1_score[current_cont_row] <- big_tb$F1_score <- f1_score
              continuous_f1_tb$true_positive[current_cont_row] <- big_tb$true_positive <- true_positive
              continuous_f1_tb$true_negative[current_cont_row] <- big_tb$true_negative <- true_negative
              continuous_f1_tb$false_negative[current_cont_row] <- big_tb$false_negative <- false_negative
              continuous_f1_tb$false_positive[current_cont_row] <- big_tb$false_positive <- false_positive
              continuous_f1_tb$num_geno[current_cont_row] <- big_tb$num_geno <- nrow(pvals_epsilon_tb)
              continuous_f1_tb$observed_gamma_value[current_cont_row] <- big_tb$observed_gamma_value <- gamma_value
              continuous_f1_tb$recall[current_cont_row] <- big_tb$recall <- recall
              continuous_f1_tb$precision[current_cont_row] <- big_tb$precision <- precision
              continuous_f1_tb$false_positive_rate[current_cont_row] <- big_tb$false_positive_rate <- false_positive_rate
              continuous_f1_tb$positive_predictive_value[current_cont_row] <- big_tb$positive_predictive_value <- positive_predictive_value
              continuous_f1_tb$observed_delta_epsilon[current_cont_row] <- big_tb$observed_delta_epsilon <- delta_epsilon
              continuous_f1_tb$observed_mean_pheno_beta[current_cont_row] <- big_tb$observed_mean_pheno_beta <- mean_beta_pheno
              continuous_f1_tb$observed_mean_geno_beta[current_cont_row] <- big_tb$observed_mean_geno_beta <- mean_beta_geno
              continuous_f1_tb$alpha_threshold[current_cont_row] <- big_tb$alpha_threshold <- alphas[p]
              continuous_f1_tb$epsilon_threshold[current_cont_row] <- big_tb$epsilon_threshold <- epsilons[q]
              
              continuous_f1_tb_and_geno <- rbind(continuous_f1_tb_and_geno, big_tb)
              
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

continuous_f1_tb_and_geno <- continuous_f1_tb_and_geno %>% filter(!is.na(phenotype_type))

discrete_f1_tb$test[discrete_f1_tb$test == "synchronous"] <- "sync"
discrete_f1_tb_and_geno$test[discrete_f1_tb_and_geno$test == "synchronous"] <- "sync"
discrete_f1_tb_and_geno <- discrete_f1_tb_and_geno %>% filter(!is.na(phenotype_type))

write_tsv(discrete_f1_tb,
          path = "../data/F1_scores_range_of_alpha_gamma_discrete.tsv",
          col_names = TRUE)

write_tsv(discrete_f1_tb_and_geno,
          path = "../data/F1_scores_by_genotype_range_of_alpha_gamma_discrete.tsv",
          col_names = TRUE)

write_tsv(continuous_f1_tb,
          path = "../data/F1_scores_range_of_alpha_gamma_continuous.tsv",
          col_names = TRUE)

write_tsv(continuous_f1_tb_and_geno,
          path = "../data/F1_scores_by_genotype_range_of_alpha_gamma_continuous.tsv",
          col_names = TRUE)

combined_f1_tb_and_geno <- rbind(discrete_f1_tb_and_geno, continuous_f1_tb_and_geno)
combined_f1_tb <- rbind(discrete_f1_tb, continuous_f1_tb)

write_tsv(combined_f1_tb,
          path = "../data/F1_scores_range_of_alpha_gamma_combined.tsv",
          col_names = TRUE)

write_tsv(combined_f1_tb_and_geno,
          path = "../data/F1_scores_by_genotype_range_of_alpha_gamma_combined.tsv",
          col_names = TRUE)

