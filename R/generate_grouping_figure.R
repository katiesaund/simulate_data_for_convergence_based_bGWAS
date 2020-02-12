source("R/schema_plot_lib.R")

set.seed(10)
tree <- ape::rcoal(n = 12)
plot(tree)
ungrouped_phenotype   <- c("a",  "a",  "a", "b", "b", "b", "c", "c",  "b",  "b",  "d",  "d")
grouped_phenotype   <- c( 1,  1,  1, 0, 0, 0, 1, 1,  0,  0,  1,  1)

# genotype             <- c( 1,  1,  0, 1, 1, 1, 1, 0,  0,  0,  1,  1)


# ungrouped phenotype
ungrouped_phenotype_recon_ER <- ape::ace(x = ungrouped_phenotype, 
                                        phy = tree, 
                                        type = "discrete",
                                        method = "ML", 
                                        marginal = FALSE, 
                                        model = "ER")
ungrouped_phenotype_recon_ARD <- ape::ace(x = ungrouped_phenotype, 
                                         phy = tree, 
                                         type = "discrete",
                                         method = "ML", 
                                         marginal = FALSE, 
                                         model = "ARD")

AIC(ungrouped_phenotype_recon_ER) < AIC(ungrouped_phenotype_recon_ARD)

disc_geno_sync_cont_edge_color[disc_geno_sync_cont_edge_color == 1] <- "blue"
disc_geno_sync_cont_edge_color[disc_geno_sync_cont_edge_color == 0] <- "black"

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_pheno_recon_edge_color,
               main = expression(paste("Phenotype Reconstruction")),
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)