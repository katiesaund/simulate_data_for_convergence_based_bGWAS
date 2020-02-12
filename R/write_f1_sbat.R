write_f1_sbat <- function(num_tree, num_pheno) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=calc_f1"),
               paste0("#SBATCH --output=calc_f1.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=esnitkin1",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=10G --time=10:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste("Rscript ../../simulate_data_for_convergence_based_bGWAS/R/calc_f1_for_range_of_alpha_and_gamma.R", num_tree, num_pheno, sep = " ")),
             paste0(getwd(), "/", "calc_f1.sbat"),
             sep = "\n")
}
