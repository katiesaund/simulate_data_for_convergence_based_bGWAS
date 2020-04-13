library(tidyverse)

df <- 
  read_tsv(file = "../data/aggregated_hogwash_data_by_genotype_range_of_alpha_gamma_combined.tsv",
           col_names = TRUE)

df$phenotype_id <- paste0("pheno_", df$phenotype_id)
df$tree_id <- paste0("tree_", df$tree_id)

df %>% 
  group_by(phenotype_type, phenotype_phylogenetic_signal, tree_id, test, phenotype_id) %>% 
  select(num_geno) %>%
  unique() %>% 
  write_tsv("summary_of_num_genotypes_per_test.tsv")
  