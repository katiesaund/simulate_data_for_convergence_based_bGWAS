write_plot_sbat <- function(path) {
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=plot_sim_data"),
               paste0("#SBATCH --output=plot_sim_data.out"),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=esnitkin1",
               "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=30G --time=24:00:00",
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/summarize_data_lib/plot_results.R")),
             paste0(getwd(), "/", "6_plot_sim_data.sbat"),
             sep = "\n")
}
