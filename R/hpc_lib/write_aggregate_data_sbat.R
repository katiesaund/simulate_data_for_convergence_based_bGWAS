# Write .sbat file to collect all of the hogwash output data from all of the 
# various runs together into one file to be used by other scripts to calculate
# summary statistics and for plotting.
write_aggregate_data_sbat <- function(num_tree, num_pheno, path) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=aggregate_hogwash_output"),
               paste0("#SBATCH --output=aggregate_hogwash_output.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=ACCOUNT_NAME",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=10G --time=7:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/summarize_data_lib/count_num_geno_per_simulation.R"), 
                     num_tree,
                     sep = " "), 
               paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/summarize_data_lib/aggregate_data_for_range_of_alpha_and_epsilon.R"), 
                     num_tree,
                     num_pheno,
                     sep = " ")),
             paste0(getwd(), "/", "3_aggregate_hogwash_output.sbat"),
             sep = "\n")
}
