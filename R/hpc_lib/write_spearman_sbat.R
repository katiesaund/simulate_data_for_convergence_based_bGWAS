# This function writes the .sbat file needed to calculate the Spearman's rank 
# correlation coefficient for -log(p) vs. epsilon for each hogwash run.
write_spearman_sbat <- function(path) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=spearman_sim_data"),
               paste0("#SBATCH --output=spearman_sim_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=ACCOUNT_NAME",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G --time=00:30:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/summarize_data_lib/calculate_spearman_coef_pval_vs_epsilon.R")),
             paste0(getwd(), "/", "4_calculate_spearman.sbat"),
             sep = "\n")
}
