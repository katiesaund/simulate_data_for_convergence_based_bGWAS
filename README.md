# Simulate data for convergence based bacterial GWAS method benchmarking paper

# Repo contents
This repository contains the R code necessary to generate simulated tree, genotype, and phenotype data for use in bGWAS. In addition to simulating data, this code also reports the amount of genotype convergence, phenotype convergence, and their intersection on a phylogenetic tree. 

These simulated data are being created to specifically benchmark the performance of [hogwash](https://github.com/katiesaund/hogwash). The results of the benchmarking will be reported in a methods paper (check back later in April for a manuscript posted to bioRxiv).

## How to use this repository: 
1. Fork or clone the repo to your computer. 
2. Simulate data: 
  a. Data simulation is slow so I ran it on University of Michigan's high performance computer cluster (HPC), which uses slurm as the scheduling system. I've included all of the scripts to create the scheduling files. You can use these with minimal changes if you use slurm, but a large overhaul of the scheduling scripts will be required if you use a different scheduler. 
  b. The two scripts that need to be run to simulate data are simulate_continuous_data.R and simulate_discrete_data.R. 
  c. Before you can run these two scripts you'll need to 
  c. To create batch scheduler job files to run these two scripts run write_all_sbat.R. 



