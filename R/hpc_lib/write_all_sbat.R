# This bash script writes the 6 sets of sbatch files needed to run data 
#  simulation, hogwash, summarize the output data, and plot all of the data.  
# To use this script users must be in their data analysis folder (where they
# will perform data simulation and hogwash), provide the relative path from that
# directory to this simulate_data_for_convergence_based_bGWAS repo, and have
# written a small file named: simulation_input_values.tsv. The entries for
# simulation_input_values.tsv should each be on a new line. They are as follows:
# 1. Number of trees.
# 2. Number of phenotypes per tree. Actual number will be twice this because we
#    generate this many Brownian motion phenotypes and this many white noise
#    phenotypes.
# 3. Number of tips per tree. 
# 4. Number of genotypes desired. The actual number will be smaller because
#    steps are taken during the simulation to remove redundant genotypes and
#    gentoypes with extreme phylogenetic signal.

# Read in user inputs: 
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
relative_path <- args[1]

source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/hpc_lib/write_hogwash_sbat_for_sim_data.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/hpc_lib/write_sim_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/hpc_lib/write_aggregate_data_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/hpc_lib/write_plot_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/hpc_lib/write_spearman_sbat.R"))
source(paste0(relative_path, "/", "simulate_data_for_convergence_based_bGWAS/R/hpc_lib/write_resource_usage_sbat.R"))

if (!file.exists("simulation_input_values.tsv")) {
  stop("This directory must contain a file called simulation_input_values.tsv")
}
user_inputs <- read_tsv("simulation_input_values.tsv", col_names = FALSE)

print(length(user_inputs$X1))

num_tree <- user_inputs$X1[1]
num_pheno <- user_inputs$X1[2]
num_tip <- user_inputs$X1[3]
num_geno <- user_inputs$X1[4]

# 1 Write .sbat files to simulate both continuous and discrete data
write_simulate_sbat(num_tree, num_pheno, num_tip, num_geno, relative_path)
# 2 Write .sbat files to run hogwash for all three tests
write_hogwash_sbat_for_sim(num_tree, num_pheno, relative_path)
# 3 Summarize hogwash output data
write_aggregate_data_sbat(num_tree, num_pheno, relative_path)
# 4 Calculate Spearman's rank correlation coefficient for -log(P) vs epsilon
write_spearman_sbat(relative_path)
# 5 Record the amount of memory and time required to run hogwash on these data
write_resource_usage_sbat(relative_path)
# 6 Make summary plots of hogwash output data
write_plot_sbat(relative_path)

# Users should run the jobs in order, not submitting one until the previous have
# completed.