# This bash script writes the 6 sbatch files needed to run data simulation, 
#   hogwash, summarize the output data, and plot all of the data.  
# To use this script users must be in their data analysis folder (where they
# will perform data simulation and hogwash), provide the relative path from that
# directory to this simulate_data_for_convergence_based_bGWAS repo, and have
# written a small file named: simulation_input_values.tsv. The entries for
# simulation_input_values.tsv should each be on a new line. They are as follows:
# 1. Number of trees.
# 2. Number of phenotypes per tree. 
# 3. Number of tips per tree. 
# 4. Number of genotypes desired. This number is approximate.   
# TODO add a better description of the relationship between the input number of genotypes and the ACTUAL number of genotypes.

# Read in user inputs: 
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
relative_path <- args[1]

source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/write_hogwash_sbat_for_sim_data.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/write_sim_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/write_f1_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/write_plot_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/write_spearman_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/write_resource_usage_sbat.R"))

user_inputs <- read_tsv("simulation_input_values.tsv", col_names = FALSE)

num_tree <- user_inputs$X1[1]
num_pheno <- user_inputs$X1[2]
num_tip <- user_inputs$X1[3]
num_geno <- user_inputs$X1[4]

#1
write_simulate_sbat(num_tree, num_pheno, num_tip, num_geno)
#2
write_hoghwash_sbat_for_sim(num_tree, num_pheno)
#3
write_f1_sbat(num_tree, num_pheno)
#4
write_spearman_sbat()
#5
write_resource_usage_sbat()
#6
write_plot_sbat()

# Users should run the jobs in order, not submitting one until the previous have
# completed.