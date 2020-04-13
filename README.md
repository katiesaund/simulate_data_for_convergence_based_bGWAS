# Simulate data for convergence based bacterial GWAS method benchmarking paper

# Repo contents
This repository contains the R code necessary to generate simulated tree, genotype, and phenotype data for use in bGWAS. In addition to simulating data, this code also reports the amount of genotype convergence, phenotype convergence, and their intersection on a phylogenetic tree. 

These simulated data are being created to specifically benchmark the performance of [hogwash](https://github.com/katiesaund/hogwash). The results of the benchmarking will be reported in a methods paper (check back later in April for a manuscript posted to bioRxiv).

## How to use this repository: 
- Fork or clone the repo to your computer. 
- Simulate data:  
  - Data simulation is slow so I ran it on University of Michigan's high performance computer cluster (HPC), which uses the slurm scheduling system. I've included all of the scripts to create the batch files. You can use these with minimal changes if you use slurm, but a larger overhaul of the scheduling scripts will be required if you use a different scheduler. 
  - The two scripts that simulate data are simulate_continuous_data.R and simulate_discrete_data.R. 
  - Before you can run these two scripts you'll need to: 
    - Create a data/ directory in the location where you'll the data simulation
    - Write a file in that data/ dir called simulation_input_values.tsv. 
      -In this file include the number of trees, number of phenotypes per tree, number of tips per tree, and number of genotypes to generate. Each value should be in that specified order and each on its own line.
    - Create the batch scheduler job files that call these two scripts. In your data/ dir run: `$ Rscript /relative/path/to/simulation_data_for_convergence_based_bGWAS/R/run write_all_sbat.R "/relative/path/to/"`
  - Submit the 1A_simulate_discrete_data.sbat and 1B_simulate_continuous_data.sbat files (ex: `$ sbatch 1B_simulate_continuous_data.sbat`)
    - Note: all of the .sbat files are written with the order of submission as the file prefix. 
  - When the two simulation jobs are finished running you'll create several kinds of data files. For each tree-phenotype pair you have you'll generate:
    - 2 Trees: 
      - simulated_continuous_tree_NUMBER.tree
      - simulated_discrete_tree_NUMBER.tree
    - 4 Genotype matrices:
      - simulated_genotype_for_continuous_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv
      - simulated_genotype_for_continuous_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv
      - simulated_genotype_for_discrete_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv
      - simulated_genotype_for_discrete_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv
    - 4 Phenotype matrices: 
      - simulated_continuous_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv
      - simulated_continuous_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv
      - simulated_discrete_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv
      - simulated_discrete_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv


