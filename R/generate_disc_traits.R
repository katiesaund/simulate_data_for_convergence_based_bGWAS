# 1. Generate coalescent trees. ------------------------------------------------
simulated_trees <- rep(list(0), num_sim_trees)
set.seed(1)
for (i in 1:num_sim_trees) {
  # Use rcoal instead of rtree to get more realistic branch lengths. 
  simulated_trees[[i]] <- ape::rcoal(n = num_sim_tips, rooted = TRUE) 
  simulated_trees[[i]]$node.label <- rep(100, Nnode(simulated_trees[[i]]))
}

# 2. Generate genotypes; on average D-stat == 0 (BM) & save genotype/tree combos
make_genotype_mat(simulated_trees, num_sim_genotypes)
  
