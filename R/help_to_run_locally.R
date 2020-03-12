suppressWarnings(library(ape))
suppressWarnings(library(caper))
suppressWarnings(library(phytools))
suppressWarnings(library(scales))

source("R/tree.R")
source("R/discrete_trait_lib.R")
source("R/continuous_trait_lib.R")
source("R/transition_edges.R")
source("R/gamma.R")
source("R/save_data.R")
source("R/keep_interesting_genotypes.R")
source("R/ancestral_reconstruction.R")
source("R/high_confidence.R")

# Initialize variables / read in user input
num_trees <- 2
num_phenos <- 2
num_tips <- 50
num_start_trait <- 50