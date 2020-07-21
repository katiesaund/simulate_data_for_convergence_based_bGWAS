# This script creates Figure 1B and 1C. 
# 1B: 4 tip tree used to introduce node, tip, and edge nomenclature
# 1C: Two trees used to introduce the concept of convergence

# Source functions ----
source("plot_lib.R")

# Simulate data -----
set.seed(10)
tree <- ape::rcoal(n = 12)
genotype_w_convergence <- c(0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1)
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
dev.off()

pdf(file = "../../img/Figure_1B_nomeclature.pdf", width = 2.5, height = 2.5)
graphics::par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))

# Tree - no convergence
graphics::plot(small_tree, 
               show.node.label = TRUE,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = TRUE,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
dev.off()
