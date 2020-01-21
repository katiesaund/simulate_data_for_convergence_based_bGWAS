# Katie Saund
# 
# Write sbat scripts for hogwash on simulated data. 

num_tree <- 10 # TODO change to user input
num_pheno <- 5 # TODO change to user input
data_dir <- "../data/"
date <- ""
perm <- 10000
fdr <- 0.10
bootstrap <- 0.70
memory <- "750m"  # TODO change to user input
time <- "10:00:00"

for (i in 1:num_tree) {
  for (j in 1:num_pheno) {
    # Continuous BM -----
    # temp_geno <- paste0(data_dir, date, "simulated_genotype_", i, ".tsv")
    # temp_tree <- paste0(data_dir, date, "simulated_tree_", i, ".tree")
    # temp_pheno <- paste0(data_dir, date, "continuous_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
    # temp_name <- paste0("continuous_pheno_BM_tree_", i, "_pheno_", j)
    # command <- paste("Rscript /nfs/esnitkin/bin_group/pipeline/Github/gwas/gwas_support/run_hogwash_from_pbs.R",
    #                 temp_pheno,
    #                 temp_geno, 
    #                 temp_tree,
    #                 temp_name,
    #                 getwd(),
    #                 perm, 
    #                 fdr,
    #                 bootstrap, 
    #                 sep = " ")
    # fname <- paste0(getwd(), "/", "hogwash_", temp_name, ".sbat")
    # writeLines(c("#!/bin/sh",
    #            paste0("#SBATCH --job-name=", temp_name),
    #            paste0("#SBATCH --output=", temp_name, ".out"),
    #            "#SBATCH --mail-user=katiephd@umich.edu",  
    #            "#SBATCH --mail-type=NONE",
    #            "#SBATCH --export=ALL",
    #            "#SBATCH --partition=standard",
    #            "#SBATCH --account=esnitkin1",
    #            "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=5g --time=01:00:00",
    #            "cd $SLURM_SUBMIT_DIR",
    #            "echo $SLURM_SUBMIT_DIR",
    #            command),
    #          fname,
    #          sep = "\n")
    #
    # Continuous WN
    # temp_pheno <- paste0(data_dir, date, "continuous_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
    # temp_name <- paste0("continuous_pheno_WN_tree_", i, "_pheno_", j)
    # command <- paste("Rscript /nfs/esnitkin/bin_group/pipeline/Github/gwas/gwas_support/run_hogwash_from_pbs.R",
    #                 temp_pheno,
    #                 temp_geno, 
    #                 temp_tree,
    #                 temp_name,
    #                 getwd(),
    #                 perm, 
    #                 fdr,
    #                 bootstrap, 
    #                 sep = " ")
    # fname <- paste0(getwd(), "/", "hogwash_", temp_name, ".sbat")
    # writeLines(c("#!/bin/sh",
    #             paste0("#SBATCH --job-name=", temp_name),
    #             paste0("#SBATCH --output=", temp_name, ".out"),
    #             "#SBATCH --mail-user=katiephd@umich.edu",  
    #             "#SBATCH --mail-type=NONE",
    #             "#SBATCH --export=ALL",
    #             "#SBATCH --partition=standard",
    #             "#SBATCH --account=esnitkin1",
    #             "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=5g --time=01:00:00",
    #             "cd $SLURM_SUBMIT_DIR",
    #             "echo $SLURM_SUBMIT_DIR",
    #             command),
    #           fname,
    #           sep = "\n")
    #

    # Discrete BM ----
    temp_tree <- paste0(data_dir, "simulated_tree_", i, ".tree")
    temp_geno <- paste0(data_dir, "simulated_genotype_for_discrete_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
    temp_pheno <- paste0(data_dir, "simulated_discrete_pheno_BM_tree_", i, "_pheno_", j, ".tsv")
    temp_name <- paste0("discrete_pheno_BM_tree_", i, "_pheno_", j)
    command <- paste("Rscript /nfs/esnitkin/bin_group/pipeline/Github/gwas/gwas_support/run_hogwash_from_pbs.R",
                     temp_pheno,
                     temp_geno, 
                     temp_tree,
                     temp_name,
                     getwd(),
                     perm, 
                     fdr,
                     bootstrap, 
                     sep = " ")
    fname <- paste0(getwd(), "/", "hogwash_", temp_name, ".sbat")
    writeLines(c("#!/bin/sh",
                 paste0("#SBATCH --job-name=", temp_name),
                 paste0("#SBATCH --output=", temp_name, ".out"),
                 "#SBATCH --mail-user=katiephd@umich.edu",  
                 "#SBATCH --mail-type=NONE",
                 "#SBATCH --export=ALL",
                 "#SBATCH --partition=standard",
                 "#SBATCH --account=esnitkin1",
                 paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                 "cd $SLURM_SUBMIT_DIR",
                 "echo $SLURM_SUBMIT_DIR",
                 command),
               fname,
               sep = "\n")
    
    # Discrete WN
    temp_geno <- paste0(data_dir, "simulated_genotype_for_discrete_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
    temp_pheno <- paste0(data_dir, "simulated_discrete_pheno_WN_tree_", i, "_pheno_", j, ".tsv")
    temp_name <- paste0("discrete_pheno_WN_tree_", i, "_pheno_", j)
    command <- paste("Rscript /nfs/esnitkin/bin_group/pipeline/Github/gwas/gwas_support/run_hogwash_from_pbs.R",
                     temp_pheno,
                     temp_geno, 
                     temp_tree,
                     temp_name,
                     getwd(),
                     perm, 
                     fdr,
                     bootstrap, 
                     sep = " ")
    fname <- paste0(getwd(), "/", "hogwash_", temp_name, ".sbat")
    writeLines(c("#!/bin/sh",
                 paste0("#SBATCH --job-name=", temp_name),
                 paste0("#SBATCH --output=", temp_name, ".out"),
                 "#SBATCH --mail-user=katiephd@umich.edu",  
                 "#SBATCH --mail-type=NONE",
                 "#SBATCH --export=ALL",
                 "#SBATCH --partition=standard",
                 "#SBATCH --account=esnitkin1",
                 paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=", memory, " --time=", time),
                 "cd $SLURM_SUBMIT_DIR",
                 "echo $SLURM_SUBMIT_DIR",
                 command),
               fname,
               sep = "\n")
  }
}
