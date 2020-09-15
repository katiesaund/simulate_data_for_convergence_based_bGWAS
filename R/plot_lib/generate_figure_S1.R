# Generate Figure S1 B, C, E, and F.
# Import functions ----
source("plot_lib.R")

# Set up tree and genotypes ----
set.seed(1)
tree <- ape::rcoal(13)
tree$node.label <- rep("100", ape::Nnode(tree))
tree$edge.length[1] <- .3
tree$edge.length[2] <- .3
tree$edge.length[3] <- .3
tree$edge.length[9] <- .3
tree$edge.length[10] <- .3
tree$edge.length[11] <- .3

# Genotype
genotype <- matrix(0, nrow = ape::Ntip(tree), ncol = 10)
row.names(genotype) <- tree$tip.label[13:1]
colnames(genotype) <- paste0("SNP", 1:10)
genotype[1,1] <-
  genotype[2, 2] <-
  genotype[3, 3] <-
  genotype[4, 4] <-
  genotype[5, 5] <-
  genotype[6, 6] <-
  genotype[7, 7] <-
  genotype[10, 8] <-
  genotype[11, 9] <- 1
genotype <- genotype[match(tree$tip.label, row.names(genotype)), ]
gene_genotype <- as.matrix(rowSums(genotype), ncol = 1)


# Post-AR stuff ----
snp_recon <- matrix(0, nrow = 12, ncol = 10)
for (i in 1:ncol(genotype)) {
  if (colSums(genotype[, i, drop = FALSE]) != 0) {
    print(genotype[, i])
    snp_geno_recon <- ape::ace(x = genotype[, i, drop = TRUE], 
                               phy = tree, 
                               type = "discrete",
                               method = "ML", 
                               marginal = FALSE, 
                               model = "ER")
    snp_geno_ML_anc_rec <-
      as.numeric(colnames(snp_geno_recon$lik.anc)[apply(snp_geno_recon$lik.anc,
                                                        1,
                                                        which.max)])
    snp_recon[, i] <- snp_geno_ML_anc_rec
  }
}
row.names(snp_recon) <- 
            c((ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
colnames(snp_recon) <- colnames(genotype)

snp_geno_tip_and_node_recon <- rbind(genotype, snp_recon)

# get trans matrix
snp_trans <- matrix(0, nrow = Nedge(tree), ncol = 10)
for (i in 1:ncol(genotype)) {
  temp_trans <- id_transition_edges_from_vec(tr = tree, 
                                            vec = genotype[, i],
                                            node_recon = snp_recon[, i], 
                                            disc_cont = "discrete")
  snp_trans[, i] <- temp_trans$transition
}

# Get edge order (instead of tip and node order)
snp_geno_recon_by_edges <- matrix(NA, nrow = Nedge(tree), ncol = 10)
for (i in 1:ncol(genotype)) {
  snp_geno_recon_by_edges[, i] <- 
    reorder_tip_and_node_to_edge(snp_geno_tip_and_node_recon[, i], tree)
}


# Pre-AR stuff ----
gene_geno_recon <- ape::ace(x = gene_genotype, 
                           phy = tree, 
                           type = "discrete",
                           method = "ML", 
                           marginal = FALSE, 
                           model = "ER")
gene_geno_ML_anc_rec <-
  as.numeric(colnames(gene_geno_recon$lik.anc)[apply(gene_geno_recon$lik.anc,
                                                    1,
                                                    which.max)])
names(gene_geno_ML_anc_rec) <- 
  c((ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))

gene_geno_tip_and_node_recon <- c(gene_genotype, gene_geno_ML_anc_rec)
temp_trans <- id_transition_edges_from_vec(tr = tree, 
                                           vec = gene_genotype,
                                           node_recon = gene_geno_ML_anc_rec, 
                                           disc_cont = "discrete")
gene_trans <- temp_trans$transition
gene_geno_recon_by_edges <- 
  reorder_tip_and_node_to_edge(gene_geno_tip_and_node_recon, tree)

# Plotting stuff ----

cex_value <- 1
edge_width <- 2.5
tip_label_log <- TRUE

# Post-AR -----------
# Set up edge colors
post_ar_pre_grouping_colors <- snp_geno_recon_by_edges
post_ar_pre_grouping_colors <- rowSums(post_ar_pre_grouping_colors)
colors <- rainbow(sum(post_ar_pre_grouping_colors))
post_ar_pre_grouping_colors[post_ar_pre_grouping_colors == 0] <- "black"
post_ar_pre_grouping_colors[post_ar_pre_grouping_colors == 1] <- colors

# Pre-AR grouping --------
pre_ar_recon_colors <- gene_geno_recon_by_edges
pre_ar_recon_colors[pre_ar_recon_colors == 0] <- "black"
pre_ar_recon_colors[pre_ar_recon_colors == 1] <- colors[1]

pre_ar_trans_colors <- gene_trans
pre_ar_trans_colors[pre_ar_trans_colors == 0] <- "black"
pre_ar_trans_colors[pre_ar_trans_colors == 1] <- colors[2:4]

pdf(file = "../../img/Pre_vs_post_ar_grouping_examples.pdf", width = 8, height = 8)
graphics::par(mfrow = c(2, 2), mar = c(1, 3, 1, 1))
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = pre_ar_recon_colors,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               main = "Pre-AR Grouping: Ancestral Reconstruction",
               cex = cex_value)

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = pre_ar_trans_colors,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               main = "Pre-AR Grouping: Transitions",
               cex = cex_value)

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = post_ar_pre_grouping_colors,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               main = "Post-AR Grouping: Ancestral Reconstruction",
               cex = cex_value)

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = post_ar_pre_grouping_colors,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               main = "Post-AR Grouping: Transitions",
               cex = cex_value)
dev.off()