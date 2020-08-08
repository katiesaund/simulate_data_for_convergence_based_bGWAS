# This script creates Figure 2 which is a schematic of how beta-genotype and 
# beta-phenotype are defined as well as their intersection. 

# Source functions ----
source("plot_lib.R")

# Simulate data -----
set.seed(10)
tree <- ape::rcoal(n = 12)
discrete_phenotype <- c( 1,  0,  0, 1, 1, 1, 0, 1,  0,  0,  1,  1)
genotype <- c( 1,  1,  0, 1, 1, 1, 1, 0,  0,  0,  1,  1)
continuous_phenotype <- c(11, 11, 14, 7, 7, 6, 3, 2, 13, 12, 11, 10)

# Discrete phenotype ancestral reconstruction and transition edges ----
discrete_phenotype_recon_ER <- ape::ace(x = discrete_phenotype, 
                                     phy = tree, 
                                     type = "discrete",
                                     method = "ML", 
                                     marginal = FALSE, 
                                     model = "ER")
discrete_phenotype_recon_ARD <- ape::ace(x = discrete_phenotype, 
                                        phy = tree, 
                                        type = "discrete",
                                        method = "ML", 
                                        marginal = FALSE, 
                                        model = "ARD")

AIC(discrete_phenotype_recon_ER) < AIC(discrete_phenotype_recon_ARD)
# use ER
disc_pheno_ML_anc_rec <-
  as.numeric(colnames(discrete_phenotype_recon_ER$lik.anc)[apply(discrete_phenotype_recon_ER$lik.anc,
                                                    1,
                                                    which.max)])
names(disc_pheno_ML_anc_rec) <- 
  c((ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
disc_pheno_tip_and_node_recon <- c(discrete_phenotype, disc_pheno_ML_anc_rec)
names(disc_pheno_tip_and_node_recon) <-
  c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))

disc_pheno_trans <-
  id_transition_edges_from_vec(tr = tree, 
                               vec = discrete_phenotype,
                               node_recon = disc_pheno_ML_anc_rec, 
                               disc_cont = "discrete")
# Continuous phenotype ancestral reconstruction and transition edges ----
cont_phenotype_recon <- ape::ace(continuous_phenotype,
                                 tree,
                                 model = "BM",
                                 type = "continuous",
                                 method = "ML",
                                 marginal = FALSE)
cont_pheno_ML_anc_rec <- cont_phenotype_recon$ace
cont_pheno_tip_and_node_recon <- c(continuous_phenotype, cont_pheno_ML_anc_rec)
names(cont_pheno_tip_and_node_recon) <- 
  c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))

named_cont_pheno <- continuous_phenotype
names(named_cont_pheno) <- tree$tip.label
cont_pheno_trans <- 
  id_transition_edges_from_vec(tr = tree, 
                               vec = continuous_phenotype,
                               node_recon = cont_pheno_ML_anc_rec, 
                               disc_cont = "continuous")

# Genotype ancestral reconstruction and transition edges ----
discrete_genotype_recon_ER <- ape::ace(x = genotype, 
                                        phy = tree, 
                                        type = "discrete",
                                        method = "ML", 
                                        marginal = FALSE, 
                                        model = "ER")
discrete_genotype_recon_ARD <- ape::ace(x = genotype, 
                                         phy = tree, 
                                         type = "discrete",
                                         method = "ML", 
                                         marginal = FALSE, 
                                         model = "ARD")

AIC(discrete_genotype_recon_ER) < AIC(discrete_genotype_recon_ARD)
# use ER
disc_geno_ML_anc_rec <-
  as.numeric(colnames(discrete_genotype_recon_ER$lik.anc)[apply(discrete_genotype_recon_ER$lik.anc,
                                                                 1,
                                                                 which.max)])
names(disc_geno_ML_anc_rec) <- 
  c((ape::Ntip(tree) + 1):(ape::Ntip(tree) + ape::Nnode(tree)))
disc_geno_tip_and_node_recon <- c(genotype, disc_geno_ML_anc_rec)
names(disc_geno_tip_and_node_recon) <- 
  c(1:sum(ape::Ntip(tree), ape::Nnode(tree)))

disc_geno_trans <- 
  id_transition_edges_from_vec(tr = tree, 
                               vec = genotype,
                               node_recon = disc_geno_ML_anc_rec, 
                               disc_cont = "discrete")
disc_geno_trans_phyc <- disc_geno_trans$trans_dir 
disc_geno_trans_phyc[disc_geno_trans_phyc == -1] <- 0
disc_geno_trans_sync_and_cont <- disc_geno_trans$transition

# Calculate delta pheno per edge -----
cont_pheno_edge_recon_mat <- 
  convert_to_edge_mat(tree, cont_pheno_tip_and_node_recon)
geno_trans_index <- which(disc_geno_trans_sync_and_cont == 1)
geno_non_trans_index <- which(disc_geno_trans_sync_and_cont == 0)

geno_trans_delta_pheno <- 
  calculate_phenotype_change_on_edge(edge_list = geno_trans_index,
                                     phenotype_by_edges = cont_pheno_edge_recon_mat)
geno_non_trans_delta_pheno <- 
  calculate_phenotype_change_on_edge(edge_list = geno_non_trans_index,
                                     phenotype_by_edges = cont_pheno_edge_recon_mat)

all_delta_pheno <- 
  calculate_phenotype_change_on_edge(edge_list = 1:ape::Nedge(tree), 
                                     phenotype_by_edges = cont_pheno_edge_recon_mat)
# Which genotype transition edges have median(delta pheno) > median(delta phenotype on non-transition edges)?
geno_non_trans_median_delta <- median(geno_non_trans_delta_pheno)

edges_with_high_delta <- which(all_delta_pheno > geno_non_trans_median_delta)
edges_with_geno_trans_and_high_delta <-
  geno_trans_index[geno_trans_delta_pheno > geno_non_trans_median_delta]

# Prep for plots ---- 
disc_pheno_recon_by_edges <- 
  reorder_tip_and_node_to_edge(disc_pheno_tip_and_node_recon, tree)

disc_pheno_recon_edge_color <- disc_pheno_recon_by_edges
disc_pheno_recon_edge_color[disc_pheno_recon_edge_color == 1] <- "hot pink"
disc_pheno_recon_edge_color[disc_pheno_recon_edge_color == 0] <- "black"

disc_pheno_recon_by_edges_phyc_beta <- disc_pheno_recon_edge_color
disc_pheno_recon_by_edges_phyc_beta[disc_pheno_recon_by_edges_phyc_beta == "hot pink"] <- "red"

disc_pheno_trans_edge_color <- disc_pheno_trans$transition
disc_pheno_trans_edge_color[disc_pheno_trans_edge_color == 1] <- "red"
disc_pheno_trans_edge_color[disc_pheno_trans_edge_color == 0] <- "black"

cont_pheno_high_delta_color <- rep("black", ape::Nedge(tree))
cont_pheno_high_delta_color[edges_with_high_delta] <- "red"

disc_geno_recon_by_edges <- 
  reorder_tip_and_node_to_edge(disc_geno_tip_and_node_recon, tree)

disc_geno_recon_edge_color <- disc_geno_recon_by_edges
disc_geno_recon_edge_color[disc_geno_recon_edge_color == 1] <- "cyan"
disc_geno_recon_edge_color[disc_geno_recon_edge_color == 0] <- "black"

disc_geno_phyc_trans_edge_color <- disc_geno_trans_phyc
disc_geno_phyc_trans_edge_color[disc_geno_phyc_trans_edge_color == 1] <- "blue"
disc_geno_phyc_trans_edge_color[disc_geno_phyc_trans_edge_color == 0] <- "black"

disc_geno_sync_cont_edge_color <- disc_geno_trans_sync_and_cont 
disc_geno_sync_cont_edge_color[disc_geno_sync_cont_edge_color == 1] <- "blue"
disc_geno_sync_cont_edge_color[disc_geno_sync_cont_edge_color == 0] <- "black"

# Calculate intersection -----
phyc_overlap_color <- 
  as.numeric(disc_pheno_recon_by_edges + disc_geno_trans_phyc == 2)
phyc_overlap_color[phyc_overlap_color == 1] <- "purple"
phyc_overlap_color[phyc_overlap_color == 0] <- "black"

sync_overlap_color <- 
  as.numeric(disc_pheno_trans$transition + disc_geno_trans_sync_and_cont  == 2)
sync_overlap_color[sync_overlap_color == 1] <- "purple"
sync_overlap_color[sync_overlap_color == 0] <- "black"

cont_overlap_color <- rep("black", ape::Nedge(tree))
cont_overlap_color[edges_with_geno_trans_and_high_delta] <- "purple"

cex_value <- 1
edge_width <- 2.5
tip_label_log <- FALSE

# Make plot ----
pdf(file = "../../img/Figure_2_hogwash_test_schema.pdf", width = 9, height = 7.5)
graphics::par(mfrow = c(3, 5), mar = c(1, 1, 1, 1))

 # phyc
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_pheno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_pheno_recon_by_edges_phyc_beta,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_geno_phyc_trans_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width,
               show.tip.label = tip_label_log,
               edge.color = phyc_overlap_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)


# sync
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_pheno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_pheno_trans_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)

graphics::plot(tree,
               font = 1,
               edge.width = edge_width,
               show.tip.label = tip_label_log,
               edge.color = disc_geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = disc_geno_sync_cont_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = sync_overlap_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)

# continuous phenotype
plot_p_recon <- phytools::contMap(tree = tree, type = "p",
                                  x = named_cont_pheno,
                                  method = "user",
                                  anc.states = cont_pheno_ML_anc_rec,
                                  plot = FALSE)

plot(plot_p_recon,
     add = TRUE,
     ylim = c(-1 / 25 * ape::Ntip(tree), ape::Ntip(tree)),
     colors = plot_p_recon$cols,
     lwd = 3,
     outline = FALSE,
     ftype = "off",
     offset = 1.7)


graphics::par(mar = c(1, 1, 1, 1))

graphics::plot(tree,
               font = 1,
               edge.width = edge_width, 
               show.tip.label = tip_label_log,
               edge.color = cont_pheno_high_delta_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)

graphics::plot(tree,
               font = 1,
               edge.width = edge_width,
               show.tip.label = tip_label_log,
               edge.color = disc_geno_recon_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               font = 1,
               show.tip.label = tip_label_log,
               edge.width = edge_width, 
               edge.color = disc_geno_sync_cont_edge_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
graphics::plot(tree,
               show.tip.label = tip_label_log,
               font = 1,
               edge.width = edge_width, 
               edge.color = cont_overlap_color,
               use.edge.length = FALSE,
               label.offset = 0.25,
               adj = 0,
               cex = cex_value)
dev.off()