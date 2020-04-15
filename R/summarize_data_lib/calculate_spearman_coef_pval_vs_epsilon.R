# This script calculates Spearman's rank correlation coefficient for the
# -log(P-value) vs epsilon for each run of hogwash. It saves the individual 
# rho values as well as summary values. 

# Inputs ----
library(tidyverse)

df <-
  read_tsv(file = "../data/aggregated_hogwash_data_by_genotype_range_of_alpha_gamma_combined.tsv", 
           col_names = TRUE)

df$phenotype_id <- paste0("pheno_", df$phenotype_id)
df$tree_id <- paste0("tree_", df$tree_id)

num_test <- length(unique(df$test))
num_tree <- length(unique(df$tree_id))
num_signal <- length(unique(df$phenotype_phylogenetic_signal))
num_pheno <- length(unique(df$phenotype_id))

alpha_to_filter_to <- df$alpha_threshold[1]
epsilon_to_filter_to <- df$epsilon_threshold[1]

num_row = num_test * num_tree * num_signal * num_pheno

# Initialize storage data.frame ----
spearman_df <- data.frame(matrix(NA, nrow = num_row, ncol = 8))
colnames(spearman_df) <- c("test",
                           "tree_id", 
                           "phenotype_id",
                           "signal", 
                           "spearman_pvalue", 
                           "spearman_rho", 
                           "nrow",
                           "report_num_geno")

row_index <- 1
for (i in 1:num_test) {
  for (j in 1:num_tree) {
    for (k in 1:num_pheno) {
      for (l in 1:num_signal) {
        
        current_test <- unique(df$test)[i]
        current_tree <- unique(df$tree_id)[j]
        current_pheno <- unique(df$phenotype_id)[k]
        current_signal <- unique(df$phenotype_phylogenetic_signal)[l]
        
        temp_data <- df %>% 
          # Subset dataframe to one specific test-tree-pheno-signal hogwash run
          filter(test == current_test, 
                 tree_id == current_tree, 
                 phenotype_id == current_pheno, 
                 phenotype_phylogenetic_signal == current_signal, 
                 alpha_threshold == alpha_to_filter_to, 
                 epsilon_threshold == epsilon_to_filter_to)
        
        # Run spearman's rank correlation coefficient
        spearman_output <- stats::cor.test(x = temp_data$fdr_corrected_pvals, 
                                           y = temp_data$epsilon, 
                                           method = "spearman", # Rank 
                                           alternative = "greater") 
        # Greater means it's testing for a positive correlation
        # NA is returned for p-value and rho if either x or y is flat (like,
        #   if all of the fdr corrected pvalues are 1)
        spearman_df[row_index, ] <- c(current_test, 
                                 current_tree, 
                                 current_pheno, 
                                 current_signal, 
                                 spearman_output$p.value, 
                                 spearman_output$estimate,
                                 nrow(temp_data), 
                                 temp_data$num_geno[1])
        row_index <- row_index + 1
      }
    }
  }
}

storage.mode(spearman_df$spearman_pvalue) <- "numeric"
storage.mode(spearman_df$spearman_rho) <- "numeric"
storage.mode(spearman_df$nrow) <- "numeric"
storage.mode(spearman_df$report_num_geno) <- "numeric"

spearman_df$spearman_pvalue[which(is.nan(spearman_df$spearman_pvalue))] <- NA

# Plot rho values ----
spearman_df %>% 
  ggplot(aes(x = test, y = spearman_rho, fill = signal)) + 
  geom_boxplot(alpha = 0.5) + 
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               position = position_dodge(0.75), 
               dotsize = 0.4) + 
  theme_bw() + 
  ylab("Spearman Rho: -log(p) ~ epsilon") + 
  xlab("") + 
  ggsave(filename = "../figures/pval_vs_epsilon_spearman_rho_boxplot.pdf")

# Save output -----
write_tsv(spearman_df,
          path = "../data/pval_vs_epsilon_spearman.tsv", 
          col_names = TRUE)

# Initialize summary data.frame ----
rho_medians_df <- data.frame(matrix(NA, nrow = num_test * num_signal, ncol = 7))
colnames(rho_medians_df) <- c("test", 
                              "signal", 
                              "median_spearman_rho",
                              "mean_spearman_rho", 
                              "max_spearman_rho",
                              "min_spearman_rho", 
                              "num_no_spearman_rho")

row_id <- 1
for (i in 1:num_test) {
  for (j in 1:num_signal) {
    current_test <- unique(df$test)[i]
    current_signal <- unique(df$phenotype_phylogenetic_signal)[j]
    
    temp_data <- spearman_df %>% 
      filter(test == current_test, 
             signal == current_signal)
    
    temp_data$spearman_rho[which(is.nan(temp_data$spearman_rho))] <- NA
    
    if (nrow(temp_data) != (num_tree * num_pheno)) { 
      stop("Dim mismatch") 
    }
    if (sum(temp_data$nrow != temp_data$report_num_geno) > 0) { 
      stop("Wrong number of genotypes")
    }

    current_median <- median(temp_data$spearman_rho, na.rm = TRUE)
    current_mean <- mean(temp_data$spearman_rho, na.rm = TRUE)
    current_max <- max(temp_data$spearman_rho, na.rm = TRUE)
    current_min <- min(temp_data$spearman_rho, na.rm = TRUE)
    current_NA <- sum(is.na(temp_data$spearman_rho))
    
    rho_medians_df[row_id, ] <- c(current_test, 
                                  current_signal, 
                                  current_median,
                                  current_mean,
                                  current_max,
                                  current_min,
                                  current_NA)
    row_id <- row_id + 1
  }
}

storage.mode(rho_medians_df$median_spearman_rho) <- "numeric"
storage.mode(rho_medians_df$mean_spearman_rho ) <- "numeric"
storage.mode(rho_medians_df$max_spearman_rho) <- "numeric"
storage.mode(rho_medians_df$min_spearman_rho) <- "numeric"
storage.mode(rho_medians_df$num_no_spearman_rho) <- "numeric"

# Save summary output ----
write_tsv(rho_medians_df,
          path = "../data/pval_vs_epsilon_spearman_rho_summaries.tsv", 
          col_names = TRUE)