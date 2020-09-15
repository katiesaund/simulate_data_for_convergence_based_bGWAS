# Generate Figure S2 A-F to provide an intuitive understanding of lamda, D, and 
# phylogenetic signal. 

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
phylo_signal_df <- as.data.frame(matrix(NA, nrow = 2, ncol = 4))
# colnames(phylo_signal_df) <- c("signal", "clumpy", "BM", "intermediate", "WN", "over_dispersed")
colnames(phylo_signal_df) <- c("signal", "BM", "intermediate", "WN")

phylo_signal_df$signal <- c("D", "lambda")
# Should be Negative
# phylo_signal_df$clumpy[1] <- gsub(pattern = ".*D = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_disc_clumpy, tree = tree)) # -1.13

# Should be zero
phylo_signal_df$BM[1] <- gsub(pattern = ".*D = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_disc_BM, tree = tree)) # 0.043

# should be 0.5
phylo_signal_df$intermediate[1] <- gsub(pattern = ".*D = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_disc_intermediate, tree = tree)) # 0.379

# Should be 1
phylo_signal_df$WN[1] <- gsub(pattern = ".*D = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_disc_WN, tree = tree)) # 0.971

# should be close to 2
# phylo_signal_df$over_dispersed[1] <- gsub(pattern = ".*D = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_disc_overdispersed, tree = tree)) # 1.782

# Continuous -----
set.seed(1)
# pheno_cont_clumpy <- matrix(rTraitCont(tree, model = "BM"), ncol = 1)
nums <- c(rep(1, 3), rep(0, 5), rep(1, 4))
nums <- nums * rnorm(n = 12, mean = 1, sd = 0.01)
pheno_cont_clumpy <- matrix(nums, ncol = 1)


set.seed(1)
pheno_cont_BM <- matrix(rTraitCont(tree, model = "BM"), ncol = 1)
pheno_cont_intermediate <- matrix(sample(rTraitCont(tree, model = "BM"), replace = FALSE), ncol = 1)
set.seed(2)
pheno_cont_WN <- matrix(sample(rTraitCont(tree, model = "BM"), replace = FALSE), ncol = 1)

set.seed(3)
#pheno_cont_overdispersed <- matrix(sample(rTraitCont(tree, model = "BM"), replace = FALSE), ncol = 1)
low_nums <- rnorm(6, mean = 0, sd = 0.01)
hi_nums <- rnorm(6, mean = 1, sd = 0.1)
mixed_nums <- c(low_nums[1], hi_nums[1], low_nums[2], hi_nums[2], 
                low_nums[3], hi_nums[3], low_nums[4], hi_nums[4], 
                low_nums[5], hi_nums[5], low_nums[6], hi_nums[6])
pheno_cont_overdispersed <- matrix(mixed_nums, ncol = 1)

row.names(pheno_cont_clumpy) <- row.names(pheno_cont_BM) <- 
  row.names(pheno_cont_intermediate) <- row.names(pheno_cont_WN) <- 
  row.names(pheno_cont_overdispersed) <- tree$tip.label

# Made pheno_cont_clumpy more clumpy
# pheno_cont_clumpy[10, 1] <- pheno_cont_clumpy[12, 1] * 1.01
# pheno_cont_clumpy[9, 1] <- pheno_cont_clumpy[12, 1] * 0.99
# pheno_cont_clumpy[8, 1] <- pheno_cont_clumpy[7, 1] * 1.01
# pheno_cont_clumpy[5, 1] <- pheno_cont_clumpy[4, 1] * 1.01
# pheno_cont_clumpy[6, 1] <- pheno_cont_clumpy[4, 1] * 0.99
# pheno_cont_clumpy[1, 1] <- pheno_cont_clumpy[2, 1] * 0.99


# Calculate lambda for each continuous phenotype 
# Should be closer to 2 # FIX!
# phylo_signal_df$clumpy[2] <- gsub(pattern = ".*lambda = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_cont_clumpy, tree = tree)) 
# hogwash::report_phylogenetic_signal(pheno = pheno_cont_clumpy, tree = tree) # 1.00

# Should be 1
phylo_signal_df$BM[2] <- gsub(pattern = ".*lambda = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_cont_BM, tree = tree)) # 1.00

# should be 0.5
phylo_signal_df$intermediate[2] <- gsub(pattern = ".*lambda = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_cont_intermediate, tree = tree)) # 0.328

# Should be 0
phylo_signal_df$WN[2] <- gsub(pattern = ".*lambda = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_cont_WN, tree = tree)) # 0

# should be negative ## FIX
# phylo_signal_df$over_dispersed[2] <- gsub(pattern = ".*lambda = ", "", hogwash::report_phylogenetic_signal(pheno = pheno_cont_overdispersed, tree = tree)) # 0


# Now plot them all together
pdf(file = "../../img/supplementary_phylogenetic_signal.pdf", 
    width = 10, height = 6)
# graphics::par(mfrow = c(2, 5), mar = c(2, 2, 2, 2))
graphics::par(mfrow = c(2, 3), mar = c(2, 2, 2, 2))

# phytools::dotTree(tree, 
#                   colors = "red",
#                   pheno_disc_clumpy,
#                   legend = TRUE, 
#                   ftype = "i")
phytools::dotTree(tree,
                  colors = "red",
                  pheno_disc_BM,
                  legend = TRUE, 
                  ftype = "i")
phytools::dotTree(tree,
                  colors = "red",
                  pheno_disc_intermediate,
                  legend = TRUE, 
                  ftype = "i")
phytools::dotTree(tree, 
                  colors = "red",
                 pheno_disc_WN,
                 legend = TRUE, 
                 ftype = "i")
# phytools::dotTree(tree,
#                   colors = "red",
#                   pheno_disc_overdispersed,
#                   legend = TRUE,
#                   ftype = "i")

# continuous
# phytools::dotTree(tree,
#                   pheno_cont_clumpy,
#                   legend = TRUE,
#                   ftype = "i")
phytools::dotTree(tree,
                  pheno_cont_BM,
                  legend = TRUE,
                  ftype = "i")
phytools::dotTree(tree,
                  pheno_cont_intermediate,
                  legend = TRUE,
                  ftype = "i")
phytools::dotTree(tree,
                  pheno_cont_WN,
                  legend = TRUE,
                  ftype = "i")
# phytools::dotTree(tree,
#                   pheno_cont_overdispersed,
#                   legend = TRUE,
#                   ftype = "i")
dev.off()

# Also save the phylogenetic signal values
write.csv(phylo_signal_df, "../../data/phylo_signal_supplementary_trees.csv")
phylo_signal_df
