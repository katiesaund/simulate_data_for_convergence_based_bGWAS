# Simulate data for convergence based bacterial GWAS method benchmarking paper

### Manuscript Title
hogwash: Three Methods for Genome-Wide Association Studies in Bacteria 

### Mansucript Authors
[Katie Saund](https://orcid.org/0000-0002-6214-6713) and Evan Snitkin

# Repo contents
This repository contains the R code necessary to generate simulated tree, genotype, and phenotype data for use in bGWAS. In addition to simulating data, this code also reports the amount of genotype convergence, phenotype convergence, and their intersection on a phylogenetic tree. 

These simulated data are being created to specifically benchmark the performance of [hogwash](https://github.com/katiesaund/hogwash). The results of the benchmarking will be reported in a methods paper (check back later in April for a manuscript posted to bioRxiv).

## How to use this repository: 
- Fork or clone the repo to your computer. 
- Simulate data:  
  - Data simulation is slow so I ran it on University of Michigan's high performance computer cluster (HPC), which uses the slurm scheduling system. I've included all of the scripts to create the batch files. You can use these with minimal changes if you use slurm, but a larger overhaul of the scheduling scripts will be required if you use a different scheduler. 
  - The two scripts that simulate data are simulate_continuous_data.R and simulate_discrete_data.R. 
  - Before you can run these two scripts you'll need to: 
    - Create a `data/` directory in the location where you'll the data simulation
    - Write a file called `data/simulation_input_values.tsv` 
      - In this file include: the number of trees, number of phenotypes per tree, number of tips per tree, and number of genotypes to generate. The values should be in that specific order and each value on its own line.
    - Create the batch scheduler job files that call these two scripts. In your `data/` dir run: `$ Rscript /relative/path/to/simulation_data_for_convergence_based_bGWAS/R/run write_all_sbat.R "/relative/path/to/"`
  - Submit the 1A_simulate_discrete_data.sbat and 1B_simulate_continuous_data.sbat files to the scheduling system (ex: `$ sbatch 1B_simulate_continuous_data.sbat`)
    - Note: all of the .sbat files are written with the order of submission as the file prefix. 
  - When the two simulation jobs are finished running you'll create several kinds of data files. For each tree-phenotype pair you have you'll generate:
    - 2 Trees: 
      - `simulated_continuous_tree_NUMBER.tree`
      - `simulated_discrete_tree_NUMBER.tree`
    - 4 Genotype matrices:
      - `simulated_genotype_for_continuous_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv`
      - `simulated_genotype_for_continuous_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv`
      - `simulated_genotype_for_discrete_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv`
      - `simulated_genotype_for_discrete_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv`
    - 4 Phenotype matrices: 
      - `simulated_continuous_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv`
      - `simulated_continuous_pheno_BM_tree_NUMBER_pheno_NUMBER.tsv`
      - `simulated_discrete_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv`
      - `simulated_discrete_pheno_WN_tree_NUMBER_pheno_NUMBER.tsv`
- Run hogwash on the simulated data
  - Install the hogwash package to your computer. 
    ```
         $ R
         > devtools::install_github("katiesaund/hogwash", ref = "1.0.0")
         > library(hogwash) 
    ```
  - Submit the hogwash jobs to the scheduler. Ex: `$ for i in 2*sbat; do sbatch $i; done`
  - For each hogwash job you'll generate: 
    - 2 files: 
      - A pdf with the output plots
      - A .rda with the output data
- Aggregate and summarize the output data from hogwash
  - Submit the various summarization jobs to the scheduler: 
    - `$ sbatch 3_calculate_F1.sbat`
      - This will generate several summary files which are used in for `4_calculate_spearman.sbat` and `6_plot_sim_data.sbat`  
    - When `3_calculate_F1.sbat` finishes running, submit the next summary job: 
      - `$ sbatch 4_calculate_spearman.sbat`
        - This will generate two summary files, this data describes the Spearman's rank correlation coefficients reported in the hogwash methods paper. 
          - `data/pval_vs_epsilon_spearman_rho_summaries.tsv` and
          - `data/pval_vs_epsilon_spearman.tsv`
    - `sbatch 5_record_resource_usage.sbat`
      - This creates a file `data/hogwash_resource_usage.csv` that grabs the amount of time and memory it took for each hogwash run. This particular `.sbat` file is highly specific to our HPC set up and is unlikely to translate to your machine.
- Plot the figures that summarize hogwash output for all of the simulated data. Generate Figure 5 from the hogwash paper. 
  - Create a `figures/` directory at the same level as `data/`
  - Submit the plotting job to the scheduler: 
    - `$ sbatch 6_plot_sim_data.sbat`
  - The plots generated in `figures` summarize the hogwash output from the simulated data in many forms. Of particular interest are the three subfigures that make up Figure 5 from the hogwash methods paper. 
    - `Fig_5A_pval_vs_epsilon_dot_plot_only_phyc.pdf` 
    - `Fig_5B_pval_vs_epsilon_dot_plot_only_sync.pdf`   
    - `Fig_5C_pval_vs_epsilon_dot_plot_only_continuous.pdf`   
  
## Questions or bugs? 
Please see the hogwash methods preprint (to be released soon) for more details on the data simulation process. Feel free to contact me at katiephd@umich.edu or open an issue here on github. 
