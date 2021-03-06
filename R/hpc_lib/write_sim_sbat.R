# This function writes the two .sbat files necessary to simulate the trees, 
# phenotypes, and genotypes. 
write_simulate_sbat <- function(num_tree, num_pheno, num_tip, num_geno, path) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=sim_disc_data"),
               paste0("#SBATCH --output=simulate_discrete_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=ACCOUNT_NAME",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=40G --time=24:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/simulate_discrete_data.R"), 
                     num_tree, 
                     num_pheno,
                     num_tip,
                     num_geno,
                     path,
                     sep = " ")),
             paste0(getwd(), "/", "1A_simulate_discrete_data.sbat"),
             sep = "\n")
  
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=sim_cont_data"),
               paste0("#SBATCH --output=simulate_continuous_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=ACCOUNT_NAME",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=22G --time=24:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/simulation_lib/simulate_continuous_data.R"),
                     num_tree,
                     num_pheno, 
                     num_tip, 
                     num_geno,
                     path, 
                     sep = " ")),
             paste0(getwd(), "/", "1B_simulate_continuous_data.sbat"),
             sep = "\n")
}
