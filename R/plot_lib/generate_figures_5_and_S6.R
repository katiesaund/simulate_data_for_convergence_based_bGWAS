# Generate Figure 5 and S6. 
# Additionally, generate other figures to exploratore the hogwash results.

library(tidyverse)
library(gridExtra)

df <- 
  read_tsv(file = "../data/aggregated_hogwash_data_by_genotype_range_of_alpha_epsilon_combined.tsv",
           col_names = TRUE)

df$phenotype_id <- paste0("pheno_", df$phenotype_id)
df$tree_id <- paste0("tree_", df$tree_id)

### FOR PAPER FIGURE 5 ----
# Figure 5A
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "phyc") %>% 
  ggplot() +
  geom_point(aes(x = epsilon, y = fdr_corrected_pvals), size = 0.5, 
             alpha = 0.1) + 
  xlab("") + 
  theme_bw() + 
  ylab("") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_phylogenetic_signal + test) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0)) + 
  ggsave("../figures/Fig_5A_pval_vs_epsilon_dot_plot_only_phyc.pdf", height = 16, width = 3.5, units = "in")

df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "phyc") %>% 
  ggplot() +
  geom_point(aes(x = epsilon, y = fdr_corrected_pvals), size = 0.5, 
             alpha = 0.1) + 
  xlab("") + 
  theme_bw() + 
  ylab("") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_phylogenetic_signal + test) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0)) + 
  ggsave("../figures/Fig_5A_pval_vs_epsilon_dot_plot_only_phyc.png", height = 16, width = 3.5, units = "in")

# Figure 5B
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "sync") %>% 
  ggplot() +
  geom_point(aes(x = epsilon, y = fdr_corrected_pvals), size = 0.5, 
             alpha = 0.1) + 
  xlab("") + 
  theme_bw() + 
  ylab("") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_phylogenetic_signal + test) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0)) + 
  ggsave("../figures/Fig_5B_pval_vs_epsilon_dot_plot_only_sync.pdf", height = 16, width = 3.5, units = "in")

df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "sync") %>% 
  ggplot() +
  geom_point(aes(x = epsilon, y = fdr_corrected_pvals), size = 0.5, 
             alpha = 0.1) + 
  xlab("") + 
  theme_bw() + 
  ylab("") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_phylogenetic_signal + test) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0)) + 
  ggsave("../figures/Fig_5B_pval_vs_epsilon_dot_plot_only_sync.png", height = 16, width = 3.5, units = "in")

# Figure 5C
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "continuous") %>% 
  ggplot() +
  geom_point(aes(x = epsilon, y = fdr_corrected_pvals), 
             size = 0.5, 
             alpha = 0.1) + 
  xlab("") + 
  theme_bw() + 
  ylab("") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_phylogenetic_signal + test) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0)) +
  ggsave("../figures/Fig_5C_pval_vs_epsilon_dot_plot_only_continuous.pdf", height = 16, width = 3.5, units = "in")

df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "continuous") %>% 
  ggplot() +
  geom_point(aes(x = epsilon, y = fdr_corrected_pvals), 
             size = 0.5, 
             alpha = 0.1) + 
  xlab("") + 
  theme_bw() + 
  ylab("") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_phylogenetic_signal + test) + 
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = .5, vjust = .5),
        axis.text.y = element_text(color = "black", size = 10, angle = 0)) +
  ggsave("../figures/Fig_5C_pval_vs_epsilon_dot_plot_only_continuous.png", height = 16, width = 3.5, units = "in")


### END FIGURE 5

# Figure S6: Resource Usage ----
usage_df <- read_csv("../data/hogwash_resource_usage.csv")

memory_plot <- usage_df %>% 
  ggplot(aes(x = test, y = `Memory (GB)`)) + 
  geom_boxplot(aes(fill = test)) + 
  geom_jitter() + 
  theme_bw() + 
  xlab("") + 
  theme(legend.position = "none") + 
  ylim(0, max(usage_df$`Memory (GB)`)) + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_fill_manual(values = c("dodgerblue", "orchid", "orange"))

time_plot <- usage_df %>% 
  ggplot(aes(x = test, y = `Time (Hours)`)) + 
  geom_boxplot(aes(fill = test)) + 
  geom_jitter() + 
  theme_bw() + 
  xlab("") + 
  ylim(0, max(usage_df$`Time (Hours)`)) + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_fill_manual(values = c("dodgerblue", "orchid", "orange"))

pdf("../figures/Figure_S3_hogwash_resource_usage.pdf", width = 8, height = 4)
grid.arrange(memory_plot, time_plot, ncol = 2)
dev.off()
# End Figure S6 ----

# The rest of the figures may be of interest to users, but were not used in the 
# hogwash methods paper ----
df %>% 
  ggplot() + 
  geom_histogram(mapping = aes(x = epsilon, fill = epsilon > 0.50)) +
  facet_grid(tree_id + test ~ phenotype_id + phenotype_phylogenetic_signal) + 
  theme_bw()
ggsave("../figures/epsilon_histogram.pdf", height = 6, width = 8.5, units = "in")

df %>% 
  filter(!is.na(observed_delta_epsilon), 
         alpha_threshold < 3, 
         epsilon_threshold < 0.150) %>% 
  select(c(observed_delta_epsilon, 
           tree_id, 
           phenotype_id,
           test, 
           phenotype_phylogenetic_signal)) %>% 
  unique() %>% 
  ggplot(mapping = aes(y = observed_delta_epsilon, 
                       x = test, 
                       col = phenotype_phylogenetic_signal)) + 
  # geom_boxplot() + 
  geom_jitter() + 
  facet_grid(tree_id ~ phenotype_id) + 
  theme_bw() 
ggsave("../figures/delta_epsilon_vs_test_tree_vs_phenotype.pdf", 
       height = 6, width = 8.5, units = "in")


alpha_label <- c("a<0.05", "a<0.01", "a<0.005", "a<0.001", 
                 "a<0.0005", "a<0.0001", "a<0.00005", "a<0.00001")
names(alpha_label) <- sort(unique(df$alpha_threshold))
epsilon_label <- c("e>0.01", "e>0.05", "e>0.10", "e>0.20", "e>0.30", "e>0.40",
                   "e>0.50", "e>0.60", "e>0.70", "e>0.80", "e>0.90", "e>0.95", 
                   "e>0.99")
names(epsilon_label) <- sort(unique(df$epsilon_threshold))

df %>% 
  ggplot() + 
  geom_point(mapping = aes(x = observed_delta_epsilon, 
                           y = F1_score, 
                           col = phenotype_phylogenetic_signal, 
                           shape = test)) + 
  theme_bw() + 
  facet_grid(alpha_threshold ~ epsilon_threshold, 
             labeller = labeller(alpha_threshold = alpha_label,
                                 epsilon_threshold = epsilon_label))
ggsave("../figures/F1_score_vs_delta_epsilon_alpha_vs_epsilon.pdf", 
       height = 6, width = 8.5, units = "in")

df %>% 
  select(F1_score, 
         phenotype_phylogenetic_signal,
         tree_id, 
         test,
         phenotype_id, 
         alpha_threshold, 
         epsilon_threshold) %>%
  unique() %>% 
  filter(!is.na(F1_score)) %>% 
  ggplot(mapping = aes(x = test, y = F1_score)) + 
  geom_jitter(aes(col = phenotype_phylogenetic_signal, shape = test)) + 
  theme_bw() + 
  facet_grid(alpha_threshold ~ epsilon_threshold, 
             labeller = labeller(alpha_threshold = alpha_label, 
                                 epsilon_threshold = epsilon_label)) + 
  labs(title = "F1 score: FDR P-value (a) by Epsilon (e)", 
       col = "Phylo Signal", 
       shape = "Test") + 
  xlab("Test type") + 
  ylab("F1 Score")
ggsave("../figures/F1_score_vs_test_alpha_vs_epsilon.pdf",
       height = 6, width = 8.5, units = "in")

# Save as line plots
df %>% 
  filter(alpha_threshold < 3) %>% 
  select(F1_score, phenotype_phylogenetic_signal, tree_id, test, 
         phenotype_id, alpha_threshold, epsilon_threshold) %>%
  unique() %>% 
  filter(!is.na(F1_score)) %>% 
  ggplot(mapping = aes(x = epsilon_threshold,
                       y = F1_score)) +  
  geom_smooth(aes(linetype = test, 
                  color = phenotype_phylogenetic_signal),
              se = FALSE) + 
  theme_bw() + 
  labs(title = "F1 score vs Epsilon (Loess line of best fit)") + 
  xlab("Epsilon Threshold") + 
  ylab("F1 Score")
ggsave("../figures/F1_score_vs_epsilon_line_plot.pdf", 
       height = 6, width = 8.5, units = "in")


# PPV
df %>% 
  select(positive_predictive_value, phenotype_phylogenetic_signal, tree_id, 
         test, phenotype_id, alpha_threshold, epsilon_threshold) %>% 
  unique() %>% 
  filter(!is.na(positive_predictive_value)) %>% 
  unique() %>% 
  ggplot() + 
  geom_jitter(mapping = aes(x = test, y = positive_predictive_value,
                            col = phenotype_phylogenetic_signal, 
                            shape = test)) + 
  theme_bw() + 
  labs(title = "PPV vs. test", col = "Phylo Signal", shape = "Test") + 
  xlab("Test") + 
  ylab("PPV") + 
  facet_grid(alpha_threshold ~ epsilon_threshold, 
             labeller = labeller(alpha_threshold = alpha_label, 
                                 epsilon_threshold = epsilon_label))
ggsave("../figures/PPV_vs_test_alpha_vs_epsilon.pdf",
       height = 6, width = 8.5, units = "in")

# Save as line plots
df %>% 
  filter(alpha_threshold < 3) %>% 
  select(positive_predictive_value, phenotype_phylogenetic_signal, tree_id,
         test, phenotype_id, alpha_threshold, epsilon_threshold) %>%
  unique() %>% 
  filter(!is.na(positive_predictive_value)) %>% 
  ggplot(mapping = aes(x = epsilon_threshold,
                       y = positive_predictive_value)) +  
  geom_smooth(aes(linetype = test, 
                  color = phenotype_phylogenetic_signal), 
              se = FALSE) + 
  theme_bw() + 
  labs(title = "PPV vs Epsilon (Loess line of best fit)") + 
  xlab("Epsilon Threshold") + 
  ylab("Positive Predictive Value")
ggsave("../figures/PPV_vs_epsilon_line_plot.pdf",
       height = 6, width = 8.5, units = "in")


df %>% 
  select(test, phenotype_phylogenetic_signal, phenotype_id,
         tree_id, phenotype_state_balance) %>% 
  unique() %>% 
  ggplot() + 
  geom_histogram(mapping = aes(x = phenotype_state_balance)) + 
  theme_bw() + 
  labs(title = "Phenotype state balance") + 
  xlab("Phenotype state balance") + 
  ylab("Count")
ggsave("../figures/Phenotype_state_balance_histrogram.pdf",
       height = 6, width = 8.5, units = "in")


# false_positive_rate
# Save as line plots
df %>% 
  filter(alpha_threshold < 3) %>% 
  select(false_positive_rate, phenotype_phylogenetic_signal, tree_id, 
         test, phenotype_id, alpha_threshold, epsilon_threshold) %>%
  unique() %>% 
  filter(!is.na(false_positive_rate)) %>% 
  ggplot(mapping = aes(x = epsilon_threshold,
                       y = false_positive_rate)) +  
  geom_smooth(aes(linetype = test, 
                  color = phenotype_phylogenetic_signal), 
              se = FALSE) + 
  theme_bw() + 
  labs(title = "FPR vs Epsilon (Loess line of best fit)") + 
  xlab("Epsilon Threshold") + 
  ylab("False Positive Rate")
ggsave("../figures/FPR_vs_epsilon_line_plot.pdf",
       height = 6, width = 8.5, units = "in")

# false_negative_rate
# Save as line plots
df %>% 
  mutate(false_negative_rate = false_negative / sum(false_negative, true_positive)) %>% 
  filter(alpha_threshold < 3) %>% 
  filter(!is.na(false_negative_rate)) %>% 
  select(false_negative_rate, phenotype_phylogenetic_signal, tree_id, test, 
         phenotype_id, epsilon_threshold) %>%
  unique() %>% 
  ggplot(mapping = aes(x = epsilon_threshold,
                       y = false_negative_rate)) +  
  geom_smooth(aes(linetype = test, 
                  color = phenotype_phylogenetic_signal), 
              se = FALSE) + 
  theme_bw() + 
  labs(title = "FNR vs Epsilon (Loess line of best fit)") + 
  xlab("Epsilon Threshold") + 
  ylab("False Negative Rate")
ggsave("../figures/FNR_vs_epsilon_line_plot.pdf", 
       height = 6, width = 8.5, units = "in")

# Save as line plots
df %>% 
  filter(alpha_threshold < 3) %>% 
  mutate(tp_over_all = 100 * true_positive / sum(true_negative, true_negative, false_positive, false_negative), 
         tn_over_all = 100 * true_negative / sum(true_negative, true_negative, false_positive, false_negative),
         fp_over_all = 100 * false_positive / sum(true_negative, true_negative, false_positive, false_negative),
         fn_over_all = 100 * false_negative / sum(true_negative, true_negative, false_positive, false_negative)) %>% 
  select(epsilon_threshold, test, phenotype_phylogenetic_signal,
         tp_over_all, tn_over_all, fp_over_all, fn_over_all) %>% 
  reshape2::melt(id.vars = c("epsilon_threshold",
                             "test", 
                             "phenotype_phylogenetic_signal"), 
                 variable.name = "metric", 
                 value.name = "percent") %>% 
  ggplot(aes(x = epsilon_threshold, y = percent, col = metric)) +  
  geom_smooth() + 
  theme_bw() + 
  labs(title = "TP/TN/FP/FN vs epsilon threshold by group") + 
  xlab("Epsilon Threshold") + 
  ylab("TP/TN/FP/FN (%)") + 
  facet_grid(test ~ phenotype_phylogenetic_signal)
ggsave("../figures/tp_tn_fp_fn_vs_epsilon_line_plot.pdf", 
       height = 6, width = 8.5, units = "in")

# Make faceted bar plot, which catalogs what percent of each epsilon category 
#  get called significant
df <- df %>% mutate(epsilon_category = epsilon)
df$epsilon_category[df$epsilon_category >= 0.66] <- "high"
df$epsilon_category[df$epsilon_category < 0.66 & df$epsilon_category >= 0.33] <- "medium"
df$epsilon_category[df$epsilon_category < 0.33] <- "low"

df$epsilon_category <- factor(df$epsilon_category, 
                              levels = c("low", "medium", "high"))
alpha <- 0.0005

df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  group_by(phenotype_type, phenotype_phylogenetic_signal, 
           tree_id, phenotype_id, epsilon_category, test) %>% 
  mutate(percent_signficant = 100 * sum(fdr_corrected_pvals > -log(alpha)) / n(), 
         percent_data = 100 * n() / num_geno) %>% 
  ggplot() +
  geom_bar(aes(x = epsilon_category,
               y = percent_signficant, 
               fill = epsilon_category), 
           stat = "identity", 
           position = "dodge") + 
  geom_point(aes(x = epsilon_category, y = percent_data), color = "blue") + 
  xlab("Genotype Epsilon") + 
  theme_bw() + 
  ylab("Percent") + 
  scale_fill_manual(values = c("pink", "red", "darkred")) + 
  ggtitle("Relative significance by genotype epsilon; dot=percent of data") +
  facet_grid(phenotype_type + phenotype_phylogenetic_signal + test ~ tree_id + phenotype_id)

ggsave("../figures/percent_signficance_by_genotype_epsilon_category_bar_plot.pdf", 
       height = 6, width = 8.5, units = "in")

# signficance vs epsilon 
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  ggplot() +
  geom_jitter(aes(x = epsilon, y = fdr_corrected_pvals)) + 
  xlab("Epsilon") + 
  theme_bw() + 
  ylab("-log(P-value)") + 
  ggtitle("P-value vs. Epsilon") +
  facet_grid(phenotype_type + phenotype_phylogenetic_signal + test ~ tree_id + phenotype_id) + 
  geom_hline(aes(yintercept = -log(0.005)), color = "red")
ggsave("../figures/pval_vs_epsilon_dot_plot_bin_and_cont.pdf",
       height = 6, width = 12, units = "in")

# signficance vs epsilon - Just continuous
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "continuous") %>% 
  ggplot() +
  geom_jitter(aes(x = epsilon, y = fdr_corrected_pvals)) + 
  xlab("Epsilon") + 
  theme_bw() + 
  ylab("-log(P-value)") + 
  ggtitle("P-value vs. Epsilon") +
  facet_grid(phenotype_type + phenotype_phylogenetic_signal + test ~ tree_id + phenotype_id) + 
  geom_hline(aes(yintercept = -log(0.05)), color = "red")
ggsave("../figures/pval_vs_epsilon_dot_plot_only_continuous.pdf", 
       height = 6, width = 12, units = "in")

# signficance vs epsilon - Just binary
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test != "continuous") %>% 
  ggplot() +
  geom_jitter(aes(x = epsilon, y = fdr_corrected_pvals)) + 
  xlab("Epsilon") + 
  theme_bw() + 
  ylab("-log(P-value)") + 
  ggtitle("P-value vs. Epsilon") +
  facet_grid(phenotype_type + phenotype_phylogenetic_signal + test ~ tree_id + phenotype_id) + 
  geom_hline(aes(yintercept = -log(0.005)), color = "red")
ggsave("../figures/pval_vs_epsilon_dot_plot_only_binary.pdf", 
       height = 6, width = 12, units = "in")

# signficance vs epsilon - Just continuous BM
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "continuous", phenotype_phylogenetic_signal == "BM") %>% 
  ggplot() +
  geom_jitter(aes(x = epsilon, y = fdr_corrected_pvals)) + 
  xlab("Epsilon") + 
  theme_bw() + 
  ylab("-log(P-value)") + 
  ggtitle("P-value vs. Epsilon") +
  facet_grid(phenotype_type + phenotype_phylogenetic_signal + test ~ tree_id + phenotype_id) + 
ggsave("../figures/pval_vs_epsilon_dot_plot_only_continuous_BM.pdf", 
       height = 6, width = 12, units = "in")


# P-value histogram
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  ggplot() +
  geom_histogram(aes(x = fdr_corrected_pvals)) + 
  xlab("-log(P-value); red = -log(0.05)") + 
  theme_bw() + 
  ylab("Count") + 
  ggtitle("P-value") +
  geom_vline(aes(xintercept = -log(0.05)), color = "red") + 
  facet_grid(tree_id + phenotype_id ~ phenotype_type + phenotype_phylogenetic_signal + test) + 
  ggsave("../figures/pval_histogram.pdf", height = 6, width = 12, units = "in")

# P-value histogram
df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test == "continuous", phenotype_phylogenetic_signal == "BM") %>% 
  ggplot() +
  geom_histogram(aes(x = fdr_corrected_pvals)) + 
  xlab("-log(P-value)") + 
  theme_bw() + 
  ylab("Count") + 
  ggtitle("P-value") +
  facet_grid(tree_id + phenotype_id ~ phenotype_type + test + phenotype_phylogenetic_signal) + 
  ggsave("../figures/pval_histogram_just_continuous_BM.pdf", 
         height = 6, width = 12, units = "in")

df %>% 
  filter(epsilon_threshold < 0.04, 
         alpha_threshold == -log(0.05)) %>% 
  filter(test != "continuous") %>% 
  ggplot() +
  geom_histogram(aes(x = fdr_corrected_pvals)) + 
  xlab("-log(P-value)") + 
  theme_bw() + 
  ylab("Count") + 
  ggtitle("P-value") +
  facet_grid(tree_id + phenotype_id ~ phenotype_type + test + phenotype_phylogenetic_signal) + 
  ggsave("../figures/pval_histogram_just_binary.pdf", 
         height = 6, width = 12, units = "in")