generate_trees <- function(num_trees, num_tips){
  # Generate a series of trees with a set number of tips
  simulated_trees <- rep(list(0), num_trees)
  set.seed(1)
  for (i in 1:num_trees) {
    # Use rcoal instead of rtree to get realistic branch lengths. 
    simulated_trees[[i]] <- ape::rcoal(n = num_tips, rooted = TRUE) 
    simulated_trees[[i]]$node.label <- rep(100, Nnode(simulated_trees[[i]]))
  }
  return(simulated_trees)
}


