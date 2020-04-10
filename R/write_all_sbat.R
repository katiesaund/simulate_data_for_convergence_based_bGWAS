# This bash script writes the 4 sbatch files I need to run on the cluster: 

# Read in user inputs: 
library(tidyverse)
source("../../simulate_data_for_convergence_based_bGWAS/R/write_hogwash_sbat_for_sim_data.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/write_sim_sbat.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/write_plot_sbat.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/write_f1_sbat.R")
source("../../simulate_data_for_convergence_based_bGWAS/R/write_spearman_sbat.R")

user_inputs <- read_tsv("simulation_input_values.tsv", col_names = FALSE)

num_tree <- user_inputs$X1[1]
num_pheno <- user_inputs$X1[2]
num_tip <- user_inputs$X1[3]
num_geno <- user_inputs$X1[4]

write_simulate_sbat(num_tree, num_pheno, num_tip, num_geno)
write_hoghwash_sbat_for_sim(num_tree, num_pheno)
write_f1_sbat(num_tree, num_pheno)
write_plot_sbat()
write_spearman_sbat()