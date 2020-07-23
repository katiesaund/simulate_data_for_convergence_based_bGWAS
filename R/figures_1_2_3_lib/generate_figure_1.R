# This script creates Figure 1B and 1C. 
# 1B: 4 tip tree used to introduce node, tip, and edge nomenclature
# 1C: Two trees used to introduce the concept of convergence

# Source functions ----
source("plot_lib.R")

# Simulate data -----
set.seed(10)
tree <- ape::rcoal(n = 12)
genotype_w_convergence <- c(0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0)
genotype_no_convergence <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)

# Genotype ancestral reconstruction and transition edges ----
make_recon_by_edges <- function(genotype, tree, color1, color2) {
  discrete_genotype_recon_ER <- ape::ace(x = genotype, 
                                         phy = tree, 
                                         type = "discrete",
                                         method = "ML", 
                                         marginal = FALSE, 
                                         model = "ER")
  
  disc_geno_ML_anc_rec <-
    as.numeric(colnames(discrete_genotype_recon_ER$lik.anc)[apply(discrete_genotype_recon_ER$lik.anc,
                                                                  1,
                                                                  which.max)])
  names(disc_geno_ML_anc_rec) <- 
    c((ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
  disc_geno_tip_and_node_recon <- c(genotype, disc_geno_ML_anc_rec)
  names(disc_geno_tip_and_node_recon) <- 
    c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))
  
  
  # Prep for plots ---- 
  disc_geno_recon_by_edges <- 
    reorder_tip_and_node_to_edge(disc_geno_tip_and_node_recon, tree)
  
  disc_geno_recon_edge_color <- disc_geno_recon_by_edges
  disc_geno_recon_edge_color[disc_geno_recon_edge_color == 1] <- color1
  disc_geno_recon_edge_color[disc_geno_recon_edge_color == 0] <- color2
  return(disc_geno_recon_edge_color)
}

disc_geno_recon_edge_color_w_convergence <- 
  make_recon_by_edges(genotype_w_convergence, tree, "hot pink", "black")
disc_geno_recon_edge_color_no_convergence <- 
  make_recon_by_edges(genotype_no_convergence, tree, "hot pink", "black")


cex_value <- 1
edge_width <- 2.5
tip_label_log <- FALSE
star <- "*"
dot_size <- 1.5

# Make tree tip dots for no convergence tree: 
blue_no_con <- purple_no_con <- green_no_con <- c("t7", "t12", "t9", "t5") 
yellow_no_con <- c("t7", "t12")
orange_no_con <- c("t2", "t11", "t6")
teal_no_con <- c("t8", "t3")

# Tree tip dots for convergence tree
blue_con <- c("t5", "t3", "t4", "t2")
purple_con <- "t2"
green_con <- c("t1", "t10", "t5") 
yellow_con <- c("t8", "t3", "t4", "t10", "t11")
orange_con <-  c("t7", "t12", "t9", "t5") 
teal_con <- c("t2", "t11")


first_offset <- 0.2
second_offset <- 0.4
third_offset <- 0.6
fourth_offset <- 0.8
# Make small tree
small_tree <- drop.tip(tree, 1:8)
small_tree$tip.label <- paste0("tip ", 4:1)
small_tree$node.label <- paste0("node ", c(1, 3, 2))

# Make plot ----
pdf(file = "../../img/Figure_1C_convergence.pdf", width = 3.6, height = 2.5)
graphics::par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))

# Tree - no convergence
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_geno_recon_edge_color_no_convergence,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% blue_no_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "deepskyblue1",
#           offset = first_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% purple_no_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "darkmagenta",
#           offset = second_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% green_no_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "palegreen3",
#           offset = third_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% yellow_no_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "yellow",
#           offset = fourth_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% orange_no_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "orange",
#           offset = first_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% teal_no_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "snow4",
#           offset = first_offset)

# Tree - with convergence
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_geno_recon_edge_color_w_convergence,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% blue_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "deepskyblue1",
#           offset = third_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% purple_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "darkmagenta",
#           offset = fourth_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% green_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "palegreen3",
#           offset = second_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% yellow_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "yellow",
#           offset = first_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% orange_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "orange",
#           offset = first_offset)
# tiplabels(tip = c(1:Ntip(tree))[tree$tip.label %in% teal_con],
#           pch = star, 
#           cex = dot_size, 
#           col = "snow4",
#           offset = second_offset)
# 

dev.off()

pdf(file = "../../img/Figure_1B_nomeclature.pdf", width = 2.5, height = 2.5)
graphics::par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))

# Tree - no convergence
graphics::plot(small_tree, 
               show.node.label = FALSE,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = FALSE,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
dev.off()
