# This function writes the .sbat file necessary to capture the resource usage
# for each hogwash run. 
write_resource_usage_sbat <- function(path) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=resource_usage"),
               paste0("#SBATCH --output=resource_usage.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=ACCOUNT_NAME",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=500M --time=00:10:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste(paste0("Rscript ", 
                            path, 
                            "/simulate_data_for_convergence_based_bGWAS/R/summarize_data_lib/capture_hogwash_resource_usage.R"), 
                     num_tree, 
                     num_pheno, 
                     path, sep = " ")),
             paste0(getwd(), "/", "5_record_resource_usage.sbat"),
             sep = "\n")
}
