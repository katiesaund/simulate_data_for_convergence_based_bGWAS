# Katie Saund
# 
# Write sbat scripts for hogwash on simulated data. 
write_hogwash_sbat_for_sim <- function(num_tree, num_pheno, path) {
  data_dir <- "../data/"
  perm <- 100000
  binary_fdr <- 0.0005
  continuous_fdr <- 0.05
  bootstrap <- 0.70
  memory <- "20G"
  time <- "200:00:00"
  
  for (i in 1:num_tree) {
    for (j in 1:num_pheno) {
      # Continuous BM ----
      temp_tree <- paste0(data_dir, "simulated_continuous_tree_", i, ".tree")
      temp_geno <- paste0(data_dir, "simulated_genotype_for_continuous_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_continuous_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("continuous_pheno_BM_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "both"
      command <- paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/hpc_lib/run_hogwash_sbatch.R "),
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       continuous_fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "2A_hogwash_", temp_name, ".sbat")
      writeLines(c("#!/bin/sh",
                   paste0("#SBATCH --job-name=", temp_name),
                   paste0("#SBATCH --output=", temp_name, ".out"),
                   "#SBATCH --mail-user=katiephd@umich.edu",  
                   "#SBATCH --mail-type=END",
                   "#SBATCH --export=ALL",
                   "#SBATCH --partition=standard",
                   "#SBATCH --account=esnitkin1",
                   paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                   "cd $SLURM_SUBMIT_DIR",
                   "echo $SLURM_SUBMIT_DIR",
                   "echo $SLURM_JOB_ID",
                   command),
                 fname,
                 sep = "\n")
      
      # Continuous WN
      temp_geno <- paste0(data_dir, "simulated_genotype_for_continuous_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_continuous_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("continuous_pheno_WN_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "both"
      command <- paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/hpc_lib/run_hogwash_sbatch.R "),
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       continuous_fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "2B_hogwash_", temp_name, ".sbat")
      writeLines(c("#!/bin/sh",
                   paste0("#SBATCH --job-name=", temp_name),
                   paste0("#SBATCH --output=", temp_name, ".out"),
                   "#SBATCH --mail-user=katiephd@umich.edu",  
                   "#SBATCH --mail-type=END",
                   "#SBATCH --export=ALL",
                   "#SBATCH --partition=standard",
                   "#SBATCH --account=esnitkin1",
                   paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                   "cd $SLURM_SUBMIT_DIR",
                   "echo $SLURM_SUBMIT_DIR",
                   "echo $SLURM_JOB_ID",
                   command),
                 fname,
                 sep = "\n")
      
      # Discrete BM ---- PHYC
      temp_tree <- paste0(data_dir, "simulated_discrete_tree_", i, ".tree")
      temp_geno <- paste0(data_dir, "simulated_genotype_for_discrete_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_discrete_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("phyc_discrete_pheno_BM_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "phyc"
      command <- paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/hpc_lib/run_hogwash_sbatch.R "),
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       binary_fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "2C_hogwash_", temp_name, ".sbat")
      writeLines(c("#!/bin/sh",
                   paste0("#SBATCH --job-name=", temp_name),
                   paste0("#SBATCH --output=", temp_name, ".out"),
                   "#SBATCH --mail-user=katiephd@umich.edu",  
                   "#SBATCH --mail-type=END",
                   "#SBATCH --export=ALL",
                   "#SBATCH --partition=standard",
                   "#SBATCH --account=esnitkin1",
                   paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                   "cd $SLURM_SUBMIT_DIR",
                   "echo $SLURM_SUBMIT_DIR",
                   "echo $SLURM_JOB_ID",
                   command),
                 fname,
                 sep = "\n")
      
      # Discrete BM ---- SYNC
      temp_tree <- paste0(data_dir, "simulated_discrete_tree_", i, ".tree")
      temp_geno <- paste0(data_dir, "simulated_genotype_for_discrete_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_discrete_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("synchronous_discrete_pheno_BM_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "synchronous"
      command <- paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/hpc_lib/run_hogwash_sbatch.R "),
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       binary_fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "2D_hogwash_", temp_name, ".sbat")
      writeLines(c("#!/bin/sh",
                   paste0("#SBATCH --job-name=", temp_name),
                   paste0("#SBATCH --output=", temp_name, ".out"),
                   "#SBATCH --mail-user=katiephd@umich.edu",  
                   "#SBATCH --mail-type=END",
                   "#SBATCH --export=ALL",
                   "#SBATCH --partition=standard",
                   "#SBATCH --account=esnitkin1",
                   paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                   "cd $SLURM_SUBMIT_DIR",
                   "echo $SLURM_SUBMIT_DIR",
                   "echo $SLURM_JOB_ID",
                   command),
                 fname,
                 sep = "\n")
      
      # Discrete WN -- PHYC
      temp_geno <- paste0(data_dir, "simulated_genotype_for_discrete_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_discrete_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("phyc_discrete_pheno_WN_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "phyc"
      command <- paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/hpc_lib/run_hogwash_sbatch.R "),
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       binary_fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "2E_hogwash_", temp_name, ".sbat")
      writeLines(c("#!/bin/sh",
                   paste0("#SBATCH --job-name=", temp_name),
                   paste0("#SBATCH --output=", temp_name, ".out"),
                   "#SBATCH --mail-user=katiephd@umich.edu",  
                   "#SBATCH --mail-type=END",
                   "#SBATCH --export=ALL",
                   "#SBATCH --partition=standard",
                   "#SBATCH --account=esnitkin1",
                   paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                   "cd $SLURM_SUBMIT_DIR",
                   "echo $SLURM_SUBMIT_DIR",
                   "echo $SLURM_JOB_ID",
                   command),
                 fname,
                 sep = "\n")
      # Discrete WN -- SYNC
      temp_geno <- paste0(data_dir, "simulated_genotype_for_discrete_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_discrete_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("synchronous_discrete_pheno_WN_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "synchronous"
      command <- paste(paste0("Rscript ", path, "/simulate_data_for_convergence_based_bGWAS/R/hpc_lib/run_hogwash_sbatch.R "),
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       binary_fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "2F_hogwash_", temp_name, ".sbat")
      writeLines(c("#!/bin/sh",
                   paste0("#SBATCH --job-name=", temp_name),
                   paste0("#SBATCH --output=", temp_name, ".out"),
                   "#SBATCH --mail-user=katiephd@umich.edu",  
                   "#SBATCH --mail-type=END",
                   "#SBATCH --export=ALL",
                   "#SBATCH --partition=standard",
                   "#SBATCH --account=esnitkin1",
                   paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                   "cd $SLURM_SUBMIT_DIR",
                   "echo $SLURM_SUBMIT_DIR",
                   "echo $SLURM_JOB_ID",
                   command),
                 fname,
                 sep = "\n")
    }
  }
}
