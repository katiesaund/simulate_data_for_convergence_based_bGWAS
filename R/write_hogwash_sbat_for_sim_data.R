# Katie Saund
# 
# Write sbat scripts for hogwash on simulated data. 
write_hoghwash_sbat_for_sim <- function(num_tree, num_pheno) {
  data_dir <- "../data/"
  perm <- 50000
  fdr <- 0.0005
  bootstrap <- 0.70
  memory <- "10G"
  time <- "200:00:00"
  
  for (i in 1:num_tree) {
    for (j in 1:num_pheno) {
      # Continuous BM ----
      temp_tree <- paste0(data_dir, "simulated_continuous_tree_", i, ".tree")
      temp_geno <- paste0(data_dir, "simulated_genotype_for_continuous_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_pheno <- paste0(data_dir, "simulated_continuous_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
      temp_name <- paste0("pheno_BM_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "both"
      command <- paste("Rscript /nfs/esnitkin/Project_Cdiff/Analysis/hogwash_methods/simulate_data_for_convergence_based_bGWAS/R/run_hogwash_sbatch.R ",
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "hogwash_continuous_", temp_name, ".sbat")
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
      temp_name <- paste0("pheno_WN_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "both"
      command <- paste("Rscript /nfs/esnitkin/Project_Cdiff/Analysis/hogwash_methods/simulate_data_for_convergence_based_bGWAS/R/run_hogwash_sbatch.R ",
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "hogwash_continuous_", temp_name, ".sbat")
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
      temp_name <- paste0("discrete_pheno_BM_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "phyc"
      command <- paste("Rscript /nfs/esnitkin/Project_Cdiff/Analysis/hogwash_methods/simulate_data_for_convergence_based_bGWAS/R/run_hogwash_sbatch.R ",
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "hogwash_phyc_", temp_name, ".sbat")
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
      temp_name <- paste0("discrete_pheno_BM_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "synchronous"
      command <- paste("Rscript /nfs/esnitkin/Project_Cdiff/Analysis/hogwash_methods/simulate_data_for_convergence_based_bGWAS/R/run_hogwash_sbatch.R ",
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "hogwash_synchronous_", temp_name, ".sbat")
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
      temp_name <- paste0("discrete_pheno_WN_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "phyc"
      command <- paste("Rscript /nfs/esnitkin/Project_Cdiff/Analysis/hogwash_methods/simulate_data_for_convergence_based_bGWAS/R/run_hogwash_sbatch.R ",
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "hogwash_phyc_", temp_name, ".sbat")
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
      temp_name <- paste0("discrete_pheno_WN_tree_", i, "_pheno_", j)
      temp_key <- NULL
      temp_test <- "synchronous"
      command <- paste("Rscript /nfs/esnitkin/Project_Cdiff/Analysis/hogwash_methods/simulate_data_for_convergence_based_bGWAS/R/run_hogwash_sbatch.R ",
                       temp_pheno,
                       temp_geno, 
                       temp_tree,
                       temp_name,
                       getwd(),
                       perm, 
                       fdr,
                       bootstrap, 
                       temp_test, 
                       temp_key,
                       sep = " ")
      fname <- paste0(getwd(), "/", "hogwash_synchronous_", temp_name, ".sbat")
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
