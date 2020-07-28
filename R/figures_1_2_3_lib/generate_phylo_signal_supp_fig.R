# Goal: generate a multipanel figure for the supplement that provides an 
#       intuitive understanding of lamda, D, and phylogenetic signal. 

# TODO: fix the remaining two continuous phenotypes and then migrate to illustrator
 
# Source functions ----
source("plot_lib.R")
library(hogwash)
library(phytools)

# Simulate data -----
set.seed(10)
tree <- ape::rcoal(n = 12)

# Discrete -----
pheno_disc_clumpy <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1), ncol = 1)
pheno_disc_BM <- matrix(c(1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0), ncol = 1)
pheno_disc_intermediate <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0), ncol = 1)
pheno_disc_WN <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1), ncol = 1)
pheno_disc_overdispersed <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1), ncol = 1)

row.names(pheno_disc_clumpy) <- row.names(pheno_disc_BM) <- 
  row.names(pheno_disc_intermediate) <- row.names(pheno_disc_WN) <- 
  row.names(pheno_disc_overdispersed) <- tree$tip.label

# Calculated D for each disc phenotype 
# Should be Negative
hogwash::report_phylogenetic_signal(pheno = pheno_disc_clumpy, tree = tree) # -1.13

# Should be zero
hogwash::report_phylogenetic_signal(pheno = pheno_disc_BM, tree = tree) # 0.043

# should be 0.5
hogwash::report_phylogenetic_signal(pheno = pheno_disc_intermediate, tree = tree) # 0.379

# Should be 1
hogwash::report_phylogenetic_signal(pheno = pheno_disc_WN, tree = tree) # 0.971

# should be close to 2
hogwash::report_phylogenetic_signal(pheno = pheno_disc_overdispersed, tree = tree) # 1.782

# Continuous -----
set.seed(1)
pheno_cont_clumpy <- matrix(rTraitCont(tree, model = "BM"), ncol = 1)
set.seed(1)
pheno_cont_BM <- matrix(rTraitCont(tree, model = "BM"), ncol = 1)
pheno_cont_intermediate <- matrix(sample(rTraitCont(tree, model = "BM"), replace = FALSE), ncol = 1)
set.seed(2)
pheno_cont_WN <- matrix(sample(rTraitCont(tree, model = "BM"), replace = FALSE), ncol = 1)
set.seed(3)
pheno_cont_overdispersed <- matrix(sample(rTraitCont(tree, model = "BM"), replace = FALSE), ncol = 1)

row.names(pheno_cont_clumpy) <- row.names(pheno_cont_BM) <- 
  row.names(pheno_cont_intermediate) <- row.names(pheno_cont_WN) <- 
  row.names(pheno_cont_overdispersed) <- tree$tip.label

# Calculated lambda for each continuous phenotype 
# Should be closer to 2 # FIX!
hogwash::report_phylogenetic_signal(pheno = pheno_cont_clumpy, tree = tree) # 1.00

# Should be 1
hogwash::report_phylogenetic_signal(pheno = pheno_cont_BM, tree = tree) # 1.00

# should be 0.5
hogwash::report_phylogenetic_signal(pheno = pheno_cont_intermediate, tree = tree) # 0.328

# Should be 0
hogwash::report_phylogenetic_signal(pheno = pheno_cont_WN, tree = tree) # 0

# should be negative ## FIX
hogwash::report_phylogenetic_signal(pheno = pheno_cont_overdispersed, tree = tree) # 0


# Now plot them all together
graphics::par(mfrow = c(2, 5))

phytools::dotTree(tree, 
                  pheno_disc_clumpy,
                  legend = FALSE, 
                  ftype = "i")
phytools::dotTree(tree,
                  pheno_disc_BM,
                  legend = FALSE, 
                  ftype = "i")
phytools::dotTree(tree,
                  pheno_disc_intermediate,
                  legend = FALSE, 
                  ftype = "i")
phytools::dotTree(tree,
                 pheno_disc_WN,
                 legend = FALSE, 
                 ftype = "i")
phytools::dotTree(tree,
                  pheno_disc_overdispersed,
                  legend = FALSE,
                  ftype = "i")

# continuous
plot_p_recon <- phytools::contMap(tree,
                                  pheno_cont_clumpy[, 1, drop = TRUE], 
                                  plot = FALSE)
graphics::plot(plot_p_recon,
               add = TRUE,
               colors = plot_p_recon$cols,
               lwd = 4,
               ftype = "i",
               offset = 1.7)


plot_p_recon <- phytools::contMap(tree,
                                  pheno_cont_BM[, 1, drop = TRUE], 
                                  plot = FALSE)
graphics::plot(plot_p_recon,
               add = TRUE,
               colors = plot_p_recon$cols,
               lwd = 4,
               ftype = "i",
               offset = 1.7)
plot_p_recon <- phytools::contMap(tree,
                                  pheno_cont_intermediate[, 1, drop = TRUE], 
                                  plot = FALSE)
graphics::plot(plot_p_recon,
               add = TRUE,
               colors = plot_p_recon$cols,
               lwd = 4,
               ftype = "i",
               offset = 1.7)
plot_p_recon <- phytools::contMap(tree,
                                  pheno_cont_WN[, 1, drop = TRUE], 
                                  plot = FALSE)
graphics::plot(plot_p_recon,
               add = TRUE,
               colors = plot_p_recon$cols,
               lwd = 4,
               ftype = "i",
               offset = 1.7)
plot_p_recon <- phytools::contMap(tree,
                                  pheno_cont_overdispersed[, 1, drop = TRUE], 
                                  plot = FALSE)
graphics::plot(plot_p_recon,
               add = TRUE,
               colors = plot_p_recon$cols,
               lwd = 4,
               ftype = "i",
               offset = 1.7)
