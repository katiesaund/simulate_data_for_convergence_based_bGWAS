source("R/figures_2_3_lib/fig_2_3_plot_lib.R")

set.seed(10)
tree <- ape::rcoal(n = 12)
ungrouped_genotype   <- c("a",  "a",  "a", "b", "b", "b", "c", "c",  "b",  "b",  "d",  "d")
grouped_genotype   <- c( 1,  1,  1, 0, 0, 0, 1, 1,  0,  0,  1,  1)

# genotype             <- c( 1,  1,  0, 1, 1, 1, 1, 0,  0,  0,  1,  1)

snp_geno_1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)
snp_geno_2 <- c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0)
snp_geno_3 <- c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
gene_geno <- snp_geno_1 + snp_geno_2 + snp_geno_3

# Ancestral reconstruction
snp_geno_1_recon <- ape::ace(x = snp_geno_1, 
                              phy = tree, 
                              type = "discrete",
                              method = "ML", 
                              marginal = FALSE, 
                              model = "ER")
snp_geno_2_recon <- ape::ace(x = snp_geno_2, 
                              phy = tree, 
                              type = "discrete",
                              method = "ML", 
                              marginal = FALSE, 
                              model = "ER")
snp_geno_3_recon <- ape::ace(x = snp_geno_3, 
                              phy = tree, 
                              type = "discrete",
                              method = "ML", 
                              marginal = FALSE, 
                              model = "ER")
gene_geno_recon <- ape::ace(x = gene_geno, 
                             phy = tree, 
                             type = "discrete",
                             method = "ML", 
                             marginal = FALSE, 
                             model = "ER")


# Find transition edges for each of the three genotypes
# SNP 1
snp_geno_1_ML_anc_rec <-
  as.numeric(colnames(snp_geno_1_recon$lik.anc)[apply(snp_geno_1_recon$lik.anc,
                                                                 1,
                                                                 which.max)])
names(snp_geno_1_ML_anc_rec) <- c( (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
snp_geno_1_tip_and_node_recon <- c(snp_geno_1, snp_geno_1_ML_anc_rec)
names(snp_geno_1_tip_and_node_recon) <- c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))

snp_geno_1_trans <- identify_transition_edges(tr = tree, 
                                               vec = snp_geno_1,
                                               node_recon = snp_geno_1_ML_anc_rec, 
                                               disc_cont = "discrete")
snp_geno_1_recon_by_edges <- reorder_tip_and_node_to_edge(snp_geno_1_tip_and_node_recon, tree)


# SNP 2
snp_geno_2_ML_anc_rec <-
  as.numeric(colnames(snp_geno_2_recon$lik.anc)[apply(snp_geno_2_recon$lik.anc,
                                                       1,
                                                       which.max)])
names(snp_geno_2_ML_anc_rec) <- c( (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
snp_geno_2_tip_and_node_recon <- c(snp_geno_2, snp_geno_2_ML_anc_rec)
names(snp_geno_2_tip_and_node_recon) <- c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))

snp_geno_2_trans <- identify_transition_edges(tr = tree, 
                                               vec = snp_geno_2,
                                               node_recon = snp_geno_2_ML_anc_rec, 
                                               disc_cont = "discrete")
snp_geno_2_recon_by_edges <- reorder_tip_and_node_to_edge(snp_geno_2_tip_and_node_recon, tree)

# SNP 3
snp_geno_3_ML_anc_rec <-
  as.numeric(colnames(snp_geno_3_recon$lik.anc)[apply(snp_geno_3_recon$lik.anc,
                                                       1,
                                                       which.max)])
names(snp_geno_3_ML_anc_rec) <- c( (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
snp_geno_3_tip_and_node_recon <- c(snp_geno_3, snp_geno_3_ML_anc_rec)
names(snp_geno_3_tip_and_node_recon) <- c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))

snp_geno_3_trans <- identify_transition_edges(tr = tree, 
                                               vec = snp_geno_3,
                                               node_recon = snp_geno_3_ML_anc_rec, 
                                               disc_cont = "discrete")
snp_geno_3_recon_by_edges <- reorder_tip_and_node_to_edge(snp_geno_3_tip_and_node_recon, tree)

# Gene geno
gene_geno_phyc_trans_index <- unique(c(which(snp_geno_1_trans$trans_dir == 1), 
                                       which(snp_geno_2_trans$trans_dir == 1),
                                       which(snp_geno_3_trans$trans_dir == 1)))

# gene_geno_ML_anc_rec <-
#   as.numeric(colnames(gene_geno_recon$lik.anc)[apply(gene_geno_recon$lik.anc,
#                                                       1,
#                                                       which.max)])
# names(gene_geno_ML_anc_rec) <- c( (ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
# gene_geno_tip_and_node_recon <- c(gene_geno, gene_geno_ML_anc_rec)
# names(gene_geno_tip_and_node_recon) <- c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))
# 
# gene_geno_trans <- identify_transition_edges(tr = tree, 
#                                               vec = gene_geno,
#                                               node_recon = gene_geno_ML_anc_rec, 
#                                               disc_cont = "discrete")
# gene_geno_recon_by_edges <- reorder_tip_and_node_to_edge(gene_geno_tip_and_node_recon, tree)

# Prep colors for plots
# Reconstruction
snp_geno_1_recon_edge_color <- snp_geno_1_recon_by_edges
snp_geno_1_recon_edge_color[snp_geno_1_recon_edge_color == 1] <- "lightpink"
snp_geno_1_recon_edge_color[snp_geno_1_recon_edge_color == 0] <- "black"

snp_geno_2_recon_edge_color <- snp_geno_2_recon_by_edges
snp_geno_2_recon_edge_color[snp_geno_2_recon_edge_color == 1] <- "skyblue1"
snp_geno_2_recon_edge_color[snp_geno_2_recon_edge_color == 0] <- "black"

snp_geno_3_recon_edge_color <- snp_geno_3_recon_by_edges
snp_geno_3_recon_edge_color[snp_geno_3_recon_edge_color == 1] <- "palegreen3"
snp_geno_3_recon_edge_color[snp_geno_3_recon_edge_color == 0] <- "black"

# Transitions
snp_geno_1_trans_edge_color <- snp_geno_1_trans$trans_dir
snp_geno_1_trans_edge_color[snp_geno_1_trans_edge_color == 1] <- "violetred2"
snp_geno_1_trans_edge_color[snp_geno_1_trans_edge_color %in% c(0, -1)] <- "black"

snp_geno_2_trans_edge_color <- snp_geno_2_trans$trans_dir
snp_geno_2_trans_edge_color[snp_geno_2_trans_edge_color == 1] <- "dodgerblue"
snp_geno_2_trans_edge_color[snp_geno_2_trans_edge_color %in% c(0, -1)] <- "black"

snp_geno_3_trans_edge_color <- snp_geno_3_trans$trans_dir
snp_geno_3_trans_edge_color[snp_geno_3_trans_edge_color == 1] <- "springgreen4"
snp_geno_3_trans_edge_color[snp_geno_3_trans_edge_color %in% c(0, -1)] <- "black"


# Gene
gene_geno_trans_phyc_edge_color <- rep("black", ape::Nedge(tree))
gene_geno_trans_phyc_edge_color[gene_geno_phyc_trans_index] <- "mediumpurple"

cex_value <- 1
edge_width <- 2.5
tip_label_log <- FALSE

pdf(file = "img/Figure_3_grouping_examples.pdf", width = 7, height = 5)
graphics::par(mfrow = c(3, 3), mar = c(1, 3, 1, 1))

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = snp_geno_1_recon_edge_color,
               #main = "SNP 1 Gene A Reconstruction",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = snp_geno_1_trans_edge_color,
               #main = "SNP 1 Gene A PhyC Transitions",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(0,
               type = 'n',
               axes = FALSE,
               ann = FALSE)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = snp_geno_2_recon_edge_color,
               #main = "SNP 2 Gene A Reconstruction",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = snp_geno_2_trans_edge_color,
               #main = "SNP 2 Gene A PhyC Transitions",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = gene_geno_trans_phyc_edge_color,
               #main = "Gene A PhyC Transitions",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = snp_geno_3_recon_edge_color,
               #main = "SNP 3 Gene A Reconstruction",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = snp_geno_3_trans_edge_color,
               #main = "SNP 3 Gene A PhyC Transitions",
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
dev.off()
