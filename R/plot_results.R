library(tidyverse)

df <- read_tsv(file = "../data/F1_scores_by_genotype_range_of_alpha_gamma.tsv", col_names = TRUE)

df$phenotype_id <- paste0("pheno_", df$phenotype_id)
df$tree_id <- paste0("tree_", df$tree_id)

df %>% 
  ggplot() + 
  geom_histogram(mapping = aes(x = epsilon, fill = epsilon > 0.50)) +
  facet_grid(tree_id + test ~ phenotype_id + phenotype_phylogenetic_signal) + 
  theme_bw()
ggsave("../figures/epsilon_histogram.pdf", height = 6, width = 8.5, units = "in")

df %>% 
  filter(!is.na(observed_delta_epsilon),alpha_threshold < 3, epsilon_threshold < 0.150) %>% 
  select(c(observed_delta_epsilon, tree_id, phenotype_id, test, phenotype_phylogenetic_signal)) %>% 
  unique() %>% 
  ggplot(mapping = aes(y = observed_delta_epsilon, x = test,col = phenotype_phylogenetic_signal)) + 
  # geom_boxplot() + 
  geom_jitter() + 
  facet_grid(tree_id ~ phenotype_id) + 
  theme_bw() 
ggsave("../figures/delta_epsilon_vs_test_tree_vs_phenotype.pdf", height = 6, width = 8.5, units = "in")


alpha_label <- c("a<0.01", "a<0.001", "a<0.0001")
names(alpha_label) <- sort(unique(df$alpha_threshold))
epsilon_label <- c("e>0.1", "e>0.2","e>0.3","e>0.4","e>0.5","e>0.6","e>0.7","e>0.8","e>0.9")
names(epsilon_label) <- sort(unique(df$epsilon_threshold))

df %>% 
  ggplot() + 
  geom_point(mapping = aes(x = observed_delta_epsilon, 
                           y = F1_score, 
                           col = phenotype_phylogenetic_signal, 
                           shape = test)) + 
  theme_bw() + 
  facet_grid(alpha_threshold ~ epsilon_threshold, labeller = labeller(alpha_threshold = alpha_label, epsilon_threshold = epsilon_label))
ggsave("../figures/F1_score_vs_delta_epsilon_alpha_vs_epsilon.pdf", height = 6, width = 8.5, units = "in")

df %>% 
  select(F1_score, phenotype_phylogenetic_signal, tree_id, test, phenotype_id, alpha_threshold, epsilon_threshold) %>%
  unique() %>% 
  filter(!is.na(F1_score)) %>% 
  ggplot(mapping = aes(x = test, y = F1_score)) + 
  geom_jitter(aes(col = phenotype_phylogenetic_signal, shape = test)) + 
  theme_bw() + 
  facet_grid(alpha_threshold ~ epsilon_threshold, labeller = labeller(alpha_threshold = alpha_label, epsilon_threshold = epsilon_label)) + 
  labs(title = "F1 score: FDR P-value (a) by Epsilon (e)", col = "Phylo Signal", shape = "Test") + 
  xlab("Test type") + 
  ylab("F1 Score")
ggsave("../figures/F1_score_vs_test_alpha_vs_epsilon.pdf", height = 6, width = 8.5, units = "in")

# PPV
df %>% 
  select(positive_predictive_value, phenotype_phylogenetic_signal, tree_id, test, phenotype_id, alpha_threshold, epsilon_threshold) %>% 
  unique() %>% 
  filter(!is.na(positive_predictive_value)) %>% 
  unique() %>% 
  ggplot() + 
  geom_jitter(mapping = aes(x = test, y = positive_predictive_value, col = phenotype_phylogenetic_signal, shape = test)) + 
  theme_bw() + 
  labs(title = "PPV vs. test", col = "Phylo Signal", shape = "Test") + 
  xlab("Test") + 
  ylab("PPV") + 
  facet_grid(alpha_threshold ~ epsilon_threshold, labeller = labeller(alpha_threshold = alpha_label, epsilon_threshold = epsilon_label))
ggsave("../figures/PPV_vs_test_alpha_vs_epsilon.pdf", height = 6, width = 8.5, units = "in")

df %>% 
  select(test, phenotype_phylogenetic_signal, phenotype_id, tree_id, phenotype_state_balance) %>% 
  unique() %>% 
  ggplot() + 
  geom_histogram(mapping = aes(x = phenotype_state_balance)) + 
  theme_bw() + 
  labs(title = "Phenotype state balance") + 
  xlab("Phenotype state balance") + 
  ylab("Count")
ggsave("../figures/Phenotype_state_balance_histrogram.pdf", height = 6, width = 8.5, units = "in")
