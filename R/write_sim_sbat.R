write_simulate_sbat <- function(num_tree, num_pheno, num_tip, num_geno) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=simulate_data"),
               paste0("#SBATCH --output=simulate_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=esnitkin1",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=5G --time=240:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste("Rscript ../../simulate_data_for_convergence_based_bGWAS/R/simulate_data.R", num_tree, num_pheno, num_tip, num_geno, sep = " ")),
             paste0(getwd(), "/", "simulate_data.sbat"),
             sep = "\n")
}
