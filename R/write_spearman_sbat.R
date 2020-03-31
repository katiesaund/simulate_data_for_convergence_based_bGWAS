write_rsq_sbat <- function() {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=rsq_sim_data"),
               paste0("#SBATCH --output=rsq_sim_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=esnitkin1",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=10G --time=1:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               "Rscript ../../simulate_data_for_convergence_based_bGWAS/R/calculate_r_squared.R"),
             paste0(getwd(), "/", "calculate_r_squared.sbat"),
             sep = "\n")
}
