library(tidyverse)

df <-
  read_tsv(file = "../data/F1_scores_by_genotype_range_of_alpha_gamma_combined.tsv", 
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
rsq_df <- data.frame(matrix(NA, nrow = num_row, ncol = 7))
colnames(rsq_df) <- c("test", "tree_id", "phenotype_id", "signal", "rsq", "nrow", "report_num_geno")

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
          filter(test == current_test, 
                 tree_id == current_tree, 
                 phenotype_id == current_pheno, 
                 phenotype_phylogenetic_signal == current_signal, 
                 alpha_threshold == alpha_to_filter_to, 
                 epsilon_threshold == epsilon_to_filter_to)
        
        model <- lm(formula = fdr_corrected_pvals ~ epsilon,
                    data = temp_data)
        
        r2 <- round(summary(model)$r.squared, 2)
        rsq_df[row_index, ] <- c(current_test, 
                                 current_tree, 
                                 current_pheno, 
                                 current_signal, 
                                 r2, 
                                 nrow(temp_data), 
                                 temp_data$num_geno[1])
        row_index <- row_index + 1
      }
    }
  }
}

rsq_df$rsq[which(is.nan(rsq_df$rsq))] <- NA

rsq_df %>% 
  ggplot(aes(x = test, y = rsq, fill = signal)) + 
  geom_boxplot(alpha = 0.5) + 
  #geom_point(alpha = 0.35, aes(color = signal, fill = signal)) + 
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               position = position_dodge(0.75), 
               dotsize = 0.4) + 
  theme_bw() + 
  ylab("R Squared: -log(p) ~ epsilon") + 
  xlab("") + 
  ggsave(filename = "../figures/pval_vs_epsilon_rsquared_boxplot.pdf")

write_tsv(rsq_df,
          path = "../data/pval_vs_epsilon_rsquared.tsv", 
          col_names = TRUE)

# get summary data
rsq_medians_df <- data.frame(matrix(NA, nrow = num_test * num_signal, ncol = 7))
colnames(rsq_medians_df) <- c("test", "signal", "median_rsq", "mean_rq", "max_rsq", "min_rsq", "num_no_rsq")

row_id <- 1
for (i in 1:num_test) {
  for (j in 1:num_signal) {
    current_test <- unique(df$test)[i]
    current_signal <- unique(df$phenotype_phylogenetic_signal)[j]
    
    temp_data <- rsq_df %>% 
      filter(test == current_test, 
             signal == current_signal)
    
    temp_data$rsq[which(is.nan(temp_data$rsq))] <- NA
    
    print(temp_data)
    if (nrow(temp_data) != (num_tree * num_pheno)) { 
      stop("Dim mismatch") 
    }
    if (sum(temp_data$nrow != temp_data$report_num_geno) > 0) { 
      stop("Wrong number of genotypes")
    }
    print("Rsq")
    print(as.numeric(temp_data$rsq))
    
    print("Median rsq plain then na.rm = TRUE")
    print(median(as.numeric(temp_data$rsq)))
    print(median(as.numeric(temp_data$rsq, na.rm = TRUE)))
    
    print("sum(is.na(temp_data$rsq), na.rm = FALSE)")
    print(sum(is.na(as.numeric(temp_data$rsq, na.rm = FALSE))))
    
    current_median <- median(as.numeric(temp_data$rsq, na.rm = TRUE))
    current_mean <- mean(as.numeric(temp_data$rsq, na.rm = TRUE))
    current_max <- max(as.numeric(temp_data$rsq, na.rm = TRUE))
    current_min <- min(as.numeric(temp_data$rsq, na.rm = TRUE))
    current_NA <- sum(is.na(as.numeric(temp_data$rsq)))
    
    print(c(current_test, 
            current_signal, 
            current_median,
            current_mean,
            current_max,
            current_min,
            current_NA))
    
    rsq_medians_df[row_id, ] <- c(current_test, 
                                  current_signal, 
                                  current_median,
                                  current_mean,
                                  current_max,
                                  current_min,
                                  current_NA)
    print( rsq_medians_df[row_id, ])
    row_id <- row_id + 1
  }
}

write_tsv(rsq_medians_df,
          path = "../data/pval_vs_epsilon_rsquared_summaries.tsv", 
          col_names = TRUE)