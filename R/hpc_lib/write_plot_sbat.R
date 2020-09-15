# This function writes the .sbat file that is necessary to plot the paper's 
# Figure 5 and S6. 
write_plot_sbat <- function(path) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=plot_sim_data"),
               paste0("#SBATCH --output=plot_sim_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=ACCOUNT_NAME",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=26G --time=4:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/plot_lib/generate_figures_5_and_S6.R")),
             paste0(getwd(), "/", "6_plot_sim_data.sbat"),
             sep = "\n")
}
